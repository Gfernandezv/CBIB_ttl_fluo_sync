"""
Extracción de recortes de imagen centrados en ROIs de ImageJ.

Este módulo es independiente del pipeline de análisis (analisis_ttl.py,
roi_status_selector.py) y no modifica su lógica. Se usa para generar,
a partir de:

- una imagen/stack .tif (la misma sobre la que se dibujaron los ROIs en ImageJ)
- un RoiSet.zip del ROI Manager (mismas ROIs, mismo orden que el CSV exportado)

un recorte cuadrado de tamaño fijo centrado en el centroide de cada ROI
seleccionado, para inspección visual rápida.

Requiere: roifile, tifffile, numpy, imageio (o tifffile para guardar PNG/TIFF).

Supuesto clave (estándar de ImageJ ROI Manager -> Multi Measure):
    El orden de los ROIs dentro de RoiSet.zip coincide con el orden de las
    columnas IntDen1, IntDen2, ... IntDenN del CSV exportado, es decir:
        ROI1 -> primer ROI en el .zip
        ROI2 -> segundo ROI en el .zip
        ...
    Si en algún momento se borraron o reordenaron ROIs en el ROI Manager
    después de exportar el CSV, este supuesto deja de cumplirse y el
    matching por posición será incorrecto. Ver `verify_roi_count` más abajo
    para una verificación mínima antes de confiar en los recortes.
"""

from pathlib import Path

import numpy as np
import roifile
import tifffile


def _roi_centroid(roi):
    """
    Centroide de un ImagejRoi, en coordenadas de píxeles de la imagen.

    Usa el promedio de las coordenadas del contorno cuando están disponibles
    (más representativo para ROIs no simétricos, ej. polígonos irregulares).
    Si el ROI no expone coordenadas de contorno (por ejemplo un rectángulo
    simple sin puntos explícitos), cae al centro del bounding box.
    """
    try:
        coords = roi.coordinates()
        if coords is not None and len(coords) > 0:
            cx = float(np.mean(coords[:, 0]))
            cy = float(np.mean(coords[:, 1]))
            return cx, cy
    except Exception:
        pass

    # Fallback: centro del bounding box (top/left/bottom/right)
    cx = (roi.left + roi.right) / 2.0
    cy = (roi.top + roi.bottom) / 2.0
    return cx, cy


def load_rois_ordered(roi_zip_path):
    """
    Carga RoiSet.zip y devuelve la lista de ROIs en el mismo orden en que
    aparecen en el archivo (que debe coincidir con el orden de columnas
    IntDen1..N del CSV exportado por ImageJ).

    Returns
    -------
    list of roifile.ImagejRoi
    """
    roi_zip_path = Path(roi_zip_path)
    if not roi_zip_path.exists():
        raise FileNotFoundError(f"No existe el archivo de ROIs: {roi_zip_path}")

    rois = roifile.roiread(str(roi_zip_path))
    if isinstance(rois, roifile.ImagejRoi):
        rois = [rois]
    return list(rois)


def verify_roi_count(roi_zip_path, expected_n_rois):
    """
    Verificación mínima: chequea que la cantidad de ROIs en el .zip coincide
    con la cantidad de ROIs distintas presentes en los datos del pipeline
    (por ejemplo df_norm['ROI'].nunique() o processing_active['ROI'].nunique()).

    No garantiza que el orden sea correcto, solo detecta el caso más común
    de desalineación: ROIs agregadas/eliminadas después de exportar el CSV.
    """
    rois = load_rois_ordered(roi_zip_path)
    n_found = len(rois)
    if n_found != expected_n_rois:
        print(
            f"ADVERTENCIA: {roi_zip_path} tiene {n_found} ROIs, "
            f"pero el pipeline reporta {expected_n_rois} ROIs distintas. "
            "El matching por posición (ROI1 -> primer ROI del zip, etc.) "
            "puede estar desalineado. Revisar manualmente antes de confiar "
            "en los recortes."
        )
        return False
    return True


def _load_stack(tif_path):
    """
    Carga la imagen/stack completa tal cual está en disco.

    Devuelve siempre un array 3D (frames, height, width). Si el archivo es
    una imagen 2D simple, se le agrega una dimensión de frame único para
    mantener la misma interfaz.
    """
    tif_path = Path(tif_path)
    if not tif_path.exists():
        raise FileNotFoundError(f"No existe la imagen: {tif_path}")

    arr = tifffile.imread(str(tif_path))

    if arr.ndim == 2:
        return arr[np.newaxis, ...]
    if arr.ndim == 3:
        # Asume (frames, height, width), convención habitual de ImageJ
        return arr

    raise ValueError(
        f"Forma de imagen no soportada: {arr.shape}. "
        "Se esperaba 2D (height, width) o 3D (frames, height, width)."
    )


def find_sample_tif(sample_dir):
    """
    Busca el único archivo .tif/.tiff dentro de sample_dir.

    El nombre del archivo de imagen varía entre muestras (ej.
    'Oct172025_20_ROI.tif'), a diferencia de RoiSet.zip que es siempre el
    mismo nombre. Esta función evita tener que especificar el nombre exacto
    por muestra.

    Raises
    ------
    FileNotFoundError
        Si no hay ningún .tif/.tiff en la carpeta.
    ValueError
        Si hay más de un .tif/.tiff y no se puede elegir sin ambigüedad.
    """
    sample_dir = Path(sample_dir)
    candidates = sorted(sample_dir.glob("*.tif")) + sorted(sample_dir.glob("*.tiff"))

    if len(candidates) == 0:
        raise FileNotFoundError(f"No se encontró ningún .tif/.tiff en: {sample_dir}")

    if len(candidates) > 1:
        raise ValueError(
            f"Se encontró más de un .tif/.tiff en {sample_dir}: "
            f"{[c.name for c in candidates]}. "
            "Especificá tif_name explícitamente para desambiguar."
        )

    return candidates[0]


def format_crop_filename(genotype, sample, roi_name, extension="tif"):
    """
    Arma el nombre de archivo para un recorte según la convención del
    laboratorio: m<mutante>_i<imagen>_R<roi> (ej. 'm27_i11_R008.tif').

    Parameters
    ----------
    genotype : str
        Genotipo tal como aparece en el pipeline, ej. 'm65'. Se le quita
        la 'm' inicial para obtener el número de mutante.
    sample : str
        Identificador de muestra tal como aparece en el pipeline, ej.
        'sample_19'. Se toma el número final para el índice de imagen.
    roi_name : str
        Nombre de ROI en la convención del pipeline, ej. 'ROI8'. Se toma
        el número y se rellena a 3 dígitos (ROI8 -> '008').
    extension : str
        Extensión de archivo sin punto, ej. 'tif' o 'png'.

    Returns
    -------
    str
        Nombre de archivo, ej. 'm65_i19_R008.tif'.

    Raises
    ------
    ValueError
        Si genotype, sample o roi_name no siguen el formato esperado.
    """
    genotype_str = str(genotype).strip()
    mutante = genotype_str[1:] if genotype_str.lower().startswith("m") else genotype_str
    if not mutante:
        raise ValueError(f"No se pudo derivar el número de mutante de genotype={genotype!r}")

    sample_str = str(sample).strip()
    digits_sample = "".join(ch for ch in sample_str if ch.isdigit())
    if not digits_sample:
        raise ValueError(f"No se pudo derivar el número de imagen de sample={sample!r}")
    # Toma el último bloque de dígitos si hay varios (ej. 'sample_19' -> '19')
    imagen = digits_sample.lstrip("0") or "0"

    roi_str = str(roi_name).strip()
    digits_roi = "".join(ch for ch in roi_str if ch.isdigit())
    if not digits_roi:
        raise ValueError(f"No se pudo derivar el número de ROI de roi_name={roi_name!r}")
    roi_padded = digits_roi.zfill(3)

    return f"m{mutante}_i{imagen}_R{roi_padded}.{extension}"


def crop_selected_rois(
    sample_dir,
    selected_roi_names,
    output_dir,
    genotype,
    sample,
    tif_name=None,
    roi_zip_name="RoiSet.zip",
    box_size=20,
    save_format="tif",
):
    """
    Genera recortes cuadrados de `box_size` x `box_size` px centrados en el
    centroide de cada ROI seleccionado, para TODOS los frames del stack, y
    los guarda en `output_dir` con el nombre m<mutante>_i<imagen>_R<roi>.<ext>
    (ver format_crop_filename). Cada recorte es un stack multi-página con la
    misma cantidad de frames que la imagen original, mismo recorte espacial
    en todos ellos, para poder ver la evolución temporal de esa zona.

    Parameters
    ----------
    sample_dir : str or Path
        Carpeta de la muestra (ej. data/Proc_data/mut27_image11), donde se
        espera encontrar la imagen .tif/.tiff y roi_zip_name.
    selected_roi_names : list of str
        Nombres de ROI a recortar, en la convención del pipeline (ej.
        ["ROI3", "ROI17"]), tal como aparecen en la columna 'ROI' de
        roi_temp_summary_active o processing_active.
    output_dir : str or Path
        Carpeta donde se guardan los recortes (ej. .../ROIs). Se crea si
        no existe.
    genotype : str
        Genotipo de la muestra, ej. 'm65'. Usado para nombrar el archivo.
    sample : str
        Identificador de muestra, ej. 'sample_19'. Usado para nombrar el
        archivo.
    tif_name : str or None
        Nombre exacto del archivo de imagen dentro de sample_dir. Si es
        None (default), se busca automáticamente el único .tif/.tiff de
        la carpeta con find_sample_tif().
    roi_zip_name : str
        Nombre del archivo RoiSet.zip dentro de sample_dir (mismo nombre
        en todas las muestras).
    box_size : int
        Lado del recorte cuadrado, en píxeles.
    save_format : {"tif", "png"}
        Formato de salida de cada recorte. "png" solo tiene sentido si el
        stack original es de un solo frame (PNG no soporta multi-página).

    Returns
    -------
    dict
        {roi_name: output_path or None}. None si el ROI no se encontró o
        el recorte no pudo generarse (se imprime el motivo en cada caso).
    """
    sample_dir = Path(sample_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    tif_path = (sample_dir / tif_name) if tif_name is not None else find_sample_tif(sample_dir)
    roi_zip_path = sample_dir / roi_zip_name

    rois = load_rois_ordered(roi_zip_path)
    stack = _load_stack(tif_path)
    n_frames, height, width = stack.shape

    half = box_size / 2.0
    results = {}

    for roi_name in selected_roi_names:
        # "ROI3" -> índice 2 (0-indexed) en la lista ordenada de ROIs
        try:
            roi_number = int(str(roi_name).replace("ROI", "").strip())
        except ValueError:
            print(f"{roi_name}: no se pudo interpretar como 'ROI<numero>'. Se omite.")
            results[roi_name] = None
            continue

        roi_index = roi_number - 1
        if roi_index < 0 or roi_index >= len(rois):
            print(
                f"{roi_name}: índice {roi_number} fuera de rango "
                f"(el .zip tiene {len(rois)} ROIs). Se omite."
            )
            results[roi_name] = None
            continue

        roi = rois[roi_index]
        cx, cy = _roi_centroid(roi)

        x0 = int(round(cx - half))
        y0 = int(round(cy - half))

        # Clampeo para no salirse de los bordes de la imagen
        x0 = max(0, min(x0, width - box_size))
        y0 = max(0, min(y0, height - box_size))

        if width < box_size or height < box_size:
            print(
                f"{roi_name}: la imagen ({width}x{height}) es más chica que "
                f"box_size={box_size}. Se omite."
            )
            results[roi_name] = None
            continue

        crop = stack[:, y0:y0 + box_size, x0:x0 + box_size]

        filename = format_crop_filename(
            genotype=genotype,
            sample=sample,
            roi_name=roi_name,
            extension=save_format,
        )
        out_path = output_dir / filename
        if save_format == "tif":
            tifffile.imwrite(str(out_path), crop, imagej=True)
        elif save_format == "png":
            import imageio.v3 as iio
            iio.imwrite(str(out_path), crop)
        else:
            raise ValueError(f"save_format no soportado: {save_format}")

        results[roi_name] = out_path

    n_ok = sum(1 for v in results.values() if v is not None)
    print(f"Recortes generados: {n_ok} / {len(selected_roi_names)} en {output_dir}")

    return results


def crop_included_rois_from_notebook(
    df_phase,
    files,
    sample,
    genotype,
    status_col="ROI_status",
    included_status=1,
    roi_zip_name="RoiSet.zip",
    box_size=20,
    save_format="tif",
    crops_subfolder="ROIs",
):
    """
    Atajo pensado para usarse directamente en la celda de inclusión manual
    de 01_preprocessing.ipynb, justo después de que df_phase quede con el
    ROI_status final.

    Toma las ROI con status_col == included_status en df_phase, y genera
    sus recortes de box_size x box_size (todos los frames del stack) usando
    la imagen y el RoiSet.zip que ya están en files["output_dir"] (la misma
    carpeta de la muestra). Los recortes se guardan dentro de la carpeta
    ROIs/ existente (la misma donde vive RoiSet.zip), junto a los archivos
    originales de ImageJ.

    Parameters
    ----------
    df_phase : pandas.DataFrame
        Tabla con columnas 'ROI' y status_col, ya actualizada con las
        marcas de inclusión/exclusión manual.
    files : dict
        Dict devuelto por ttl.process_sample()["files"] o equivalente,
        debe contener la clave "output_dir" (carpeta de la muestra).
    sample : str
        Identificador de muestra, ej. 'sample_19'. Usado para nombrar
        los archivos de salida.
    genotype : str
        Genotipo de la muestra, ej. 'm65'. Usado para nombrar los
        archivos de salida.
    status_col : str
        Columna de estatus de inclusión en df_phase.
    included_status : int
        Valor de status_col que indica ROI incluida (por defecto 1).
    roi_zip_name : str
        Nombre del archivo RoiSet.zip dentro de la carpeta de la muestra.
    box_size : int
        Lado del recorte cuadrado, en píxeles.
    save_format : {"tif", "png"}
        Formato de salida de cada recorte.
    crops_subfolder : str
        Nombre de la subcarpeta dentro de la carpeta de la muestra donde
        se guardan los recortes (por defecto 'ROIs', la misma carpeta
        donde vive RoiSet.zip).

    Returns
    -------
    dict
        {roi_name: output_path or None}, igual que crop_selected_rois().
    """
    if "ROI" not in df_phase.columns:
        raise ValueError("df_phase debe contener una columna 'ROI'.")
    if status_col not in df_phase.columns:
        raise ValueError(f"df_phase no tiene la columna de estatus '{status_col}'.")

    included_rois = (
        df_phase.loc[df_phase[status_col] == included_status, "ROI"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    if not included_rois:
        print(
            f"No hay ROIs con {status_col}={included_status} en df_phase. "
            "No se generaron recortes."
        )
        return {}

    sample_dir = Path(files["output_dir"])
    output_dir = sample_dir / crops_subfolder

    return crop_selected_rois(
        sample_dir=sample_dir,
        selected_roi_names=included_rois,
        output_dir=output_dir,
        genotype=genotype,
        sample=sample,
        roi_zip_name=roi_zip_name,
        box_size=box_size,
        save_format=save_format,
    )
