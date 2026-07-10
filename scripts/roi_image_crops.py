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

import re
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


def _resolve_roi_index(roi_name, n_rois):
    """
    "ROI3" -> índice 2 (0-indexed) en la lista ordenada de ROIs, o None (con
    motivo impreso) si el nombre no es interpretable o queda fuera de rango.
    """
    try:
        roi_number = int(str(roi_name).replace("ROI", "").strip())
    except ValueError:
        print(f"{roi_name}: no se pudo interpretar como 'ROI<numero>'. Se omite.")
        return None

    roi_index = roi_number - 1
    if roi_index < 0 or roi_index >= n_rois:
        print(
            f"{roi_name}: índice {roi_number} fuera de rango "
            f"(el .zip tiene {n_rois} ROIs). Se omite."
        )
        return None

    return roi_index


def _crop_roi_stack(stack, roi, box_size):
    """
    Recorte cuadrado de box_size x box_size centrado en el centroide de roi,
    para todos los frames de stack (frames, height, width). Devuelve
    (crop, None) o (None, motivo) si la imagen es más chica que box_size.
    """
    n_frames, height, width = stack.shape
    if width < box_size or height < box_size:
        return None, f"la imagen ({width}x{height}) es más chica que box_size={box_size}"

    half = box_size / 2.0
    cx, cy = _roi_centroid(roi)

    x0 = int(round(cx - half))
    y0 = int(round(cy - half))

    # Clampeo para no salirse de los bordes de la imagen
    x0 = max(0, min(x0, width - box_size))
    y0 = max(0, min(y0, height - box_size))

    return stack[:, y0:y0 + box_size, x0:x0 + box_size], None


def sort_stack_by_temperature(stack, frame_temp_df, frame_col="frame", temp_col="temp_mean"):
    """
    Reordena los frames de un stack (frames, height, width) de menor a mayor
    temperatura, mezclando heating y cooling (se ordena solo por temp_col,
    sin distinguir fase).

    frame_temp_df trae una fila por frame del stack ORIGINAL (antes de
    cualquier recorte), con frame_col 1-indexado igual que la columna
    'frame' del pipeline (ver long_format en analisis_ttl.py) y temp_col la
    temperatura asociada a ese frame. Es la misma tabla para todas las ROI
    de una muestra, ya que la temperatura es una propiedad del frame, no
    del ROI: por ejemplo
    `preprocessed_all[preprocessed_all["source_folder"] == folder][["frame", "temp_mean"]].drop_duplicates()`.

    Frames de stack que no aparezcan en frame_temp_df (ej. descartados por
    no tener TTL, ver long_format) quedan excluidos del resultado.

    Parameters
    ----------
    stack : np.ndarray
        Array (frames, height, width), frame 1 en índice 0.
    frame_temp_df : pandas.DataFrame
        Con columnas frame_col y temp_col.
    frame_col : str
    temp_col : str

    Returns
    -------
    sorted_stack : np.ndarray
        Mismo recorte espacial, frames reordenados por temperatura.
    frame_order : np.ndarray
        Numero de frame (1-indexado) en el orden usado.
    temp_order : np.ndarray
        temp_mean en el mismo orden que frame_order.
    """
    n_frames = stack.shape[0]

    temp_lookup = (
        frame_temp_df[[frame_col, temp_col]]
        .dropna()
        .drop_duplicates(subset=frame_col)
        .set_index(frame_col)[temp_col]
    )

    valid_frames = [f for f in temp_lookup.index if 1 <= f <= n_frames]
    if not valid_frames:
        raise ValueError("Ningun frame de frame_temp_df cae dentro del rango del stack.")

    order = sorted(valid_frames, key=lambda f: temp_lookup[f])
    frame_order = np.array(order, dtype=int)
    temp_order = temp_lookup.loc[order].to_numpy()
    sorted_stack = stack[frame_order - 1]

    return sorted_stack, frame_order, temp_order


def select_indices_by_temp_bin(temps, bin_width=0.5):
    """
    Para un array de temperaturas (ej. ya ordenado ascendente, como
    temp_order de sort_stack_by_temperature), elige un indice representativo
    por bin de bin_width grados: el mas cercano al centro de cada bin.

    Pensado para reducir un stack ordenado por temperatura (con muchos
    frames pegados en las zonas donde la señal cambia poco) a un frame por
    bin_width °C, para que un collage/montage quede parejo en todo el rango
    de temperatura en vez de sobre-representar las zonas con mas frames.

    Parameters
    ----------
    temps : array-like
        Temperaturas, una por frame/indice.
    bin_width : float
        Ancho de bin en grados.

    Returns
    -------
    np.ndarray
        Indices (en el orden original de temps) elegidos, ordenados de
        menor a mayor temperatura, uno por bin no vacio.
    """
    temps = np.asarray(temps, dtype=float)
    bin_idx = np.floor(temps / bin_width).astype(int)
    bin_center = (bin_idx + 0.5) * bin_width
    dist = np.abs(temps - bin_center)

    order = np.lexsort((dist, bin_idx))
    _, first_pos = np.unique(bin_idx[order], return_index=True)
    picked = order[first_pos]

    return picked[np.argsort(temps[picked])]


def montage_from_stack(stack, ncols=None, padding=1, pad_value=None):
    """
    Arma un collage 2D de un stack (frames, height, width): cada frame queda
    como una celda de una grilla, para ver de un vistazo como cambia la
    intensidad a lo largo del stack sin depender de un visor de TIFF
    interactivo (Python no trae uno equivalente al de ImageJ/Fiji; esta es
    la alternativa mas simple para recortes chicos, via
    skimage.util.montage).

    Parameters
    ----------
    stack : np.ndarray
        (frames, height, width).
    ncols : int or None
        Columnas de la grilla. None (default) arma una grilla lo mas
        cuadrada posible.
    padding : int
        Pixeles de separacion entre celdas.
    pad_value : int/float or None
        Valor de relleno del padding y de las celdas vacias al final de la
        grilla (si n_frames no completa la grilla). None (default) usa el
        minimo del stack, para que el padding no se vea como senal.

    Returns
    -------
    np.ndarray
        Imagen 2D (grilla completa), mismo dtype que stack.
    """
    from skimage.util import montage as skimage_montage

    n_frames = stack.shape[0]
    if ncols is None:
        ncols = int(np.ceil(np.sqrt(n_frames)))
    nrows = int(np.ceil(n_frames / ncols))

    if pad_value is None:
        pad_value = stack.min()

    return skimage_montage(
        stack,
        grid_shape=(nrows, ncols),
        padding_width=padding,
        fill=pad_value,
    )


def save_montage_png(stack, out_path, ncols=None, padding=1):
    """
    Arma el collage de montage_from_stack y lo guarda como PNG de 8 bits,
    escalando linealmente min-max del stack completo a 0-255 (para que el
    contraste sea comparable entre celdas de un mismo collage).

    Returns
    -------
    Path
        out_path.
    """
    import imageio.v3 as iio

    grid = montage_from_stack(stack, ncols=ncols, padding=padding, pad_value=stack.min())

    lo, hi = float(stack.min()), float(stack.max())
    if hi > lo:
        scaled = ((grid.astype(np.float32) - lo) / (hi - lo) * 255.0).clip(0, 255).astype(np.uint8)
    else:
        scaled = np.zeros_like(grid, dtype=np.uint8)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    iio.imwrite(str(out_path), scaled)
    return out_path


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

    results = {}

    for roi_name in selected_roi_names:
        roi_index = _resolve_roi_index(roi_name, len(rois))
        if roi_index is None:
            results[roi_name] = None
            continue

        crop, error = _crop_roi_stack(stack, rois[roi_index], box_size)
        if crop is None:
            print(f"{roi_name}: {error}. Se omite.")
            results[roi_name] = None
            continue

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


def crop_selected_rois_sorted_by_temperature(
    sample_dir,
    selected_roi_names,
    frame_temp_df,
    output_dir,
    genotype,
    sample,
    tif_name=None,
    roi_zip_name="RoiSet.zip",
    box_size=20,
    save_format="tif",
    frame_col="frame",
    temp_col="temp_mean",
    filename_suffix="_sorted_by_temp",
):
    """
    Igual que crop_selected_rois, pero los frames de cada recorte quedan
    reordenados de menor a mayor temperatura (heating y cooling mezclados,
    ordenados solo por temp_col) en vez del orden de adquisicion original.
    Util para ver la intensidad de senal en funcion de la temperatura como
    pelicula en ImageJ.

    Ver sort_stack_by_temperature para el formato esperado de frame_temp_df
    (una fila por frame del stack original, comun a todas las ROI de la
    muestra).

    Parameters
    ----------
    frame_temp_df : pandas.DataFrame
        Con columnas frame_col y temp_col, ej.
        preprocessed_all[preprocessed_all["source_folder"] == folder][["frame", "temp_mean"]].drop_duplicates()
    filename_suffix : str
        Sufijo agregado al nombre de archivo (antes de la extension) para
        distinguir estos recortes de los de crop_selected_rois.

    Returns
    -------
    dict
        {roi_name: output_path or None}, igual que crop_selected_rois().
    """
    sample_dir = Path(sample_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    tif_path = (sample_dir / tif_name) if tif_name is not None else find_sample_tif(sample_dir)
    roi_zip_path = sample_dir / roi_zip_name

    rois = load_rois_ordered(roi_zip_path)
    stack = _load_stack(tif_path)

    results = {}

    for roi_name in selected_roi_names:
        roi_index = _resolve_roi_index(roi_name, len(rois))
        if roi_index is None:
            results[roi_name] = None
            continue

        crop, error = _crop_roi_stack(stack, rois[roi_index], box_size)
        if crop is None:
            print(f"{roi_name}: {error}. Se omite.")
            results[roi_name] = None
            continue

        try:
            sorted_crop, _frame_order, _temp_order = sort_stack_by_temperature(
                crop, frame_temp_df, frame_col=frame_col, temp_col=temp_col
            )
        except ValueError as exc:
            print(f"{roi_name}: {exc}. Se omite.")
            results[roi_name] = None
            continue

        filename = format_crop_filename(
            genotype=genotype,
            sample=sample,
            roi_name=roi_name,
            extension=save_format,
        )
        stem, dot, ext = filename.rpartition(".")
        out_path = output_dir / f"{stem}{filename_suffix}.{ext}"

        if save_format == "tif":
            tifffile.imwrite(str(out_path), sorted_crop, imagej=True)
        elif save_format == "png":
            import imageio.v3 as iio
            iio.imwrite(str(out_path), sorted_crop)
        else:
            raise ValueError(f"save_format no soportado: {save_format}")

        results[roi_name] = out_path

    n_ok = sum(1 for v in results.values() if v is not None)
    print(f"Recortes ordenados por temperatura generados: {n_ok} / {len(selected_roi_names)} en {output_dir}")

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


def list_rois_from_crops_folder(files, crops_subfolder="ROIs", extension="tif"):
    """
    Lee la carpeta de recortes (generada por crop_included_rois_from_notebook)
    y devuelve la lista de ROI que siguen ahí, en formato 'ROI<numero>'.

    Pensado para el flujo: generar recortes de un set amplio de candidatas,
    borrar a mano en el Finder/Explorador los recortes que no se ven bien, y
    releer la carpeta para obtener el roi_valid_manual final sin transcribir
    números a mano.

    Parameters
    ----------
    files : dict
        Dict devuelto por ttl.process_sample()["files"] (o equivalente),
        debe contener la clave "output_dir" (carpeta de la muestra).
    crops_subfolder : str
        Nombre de la subcarpeta donde viven los recortes (ver
        crop_included_rois_from_notebook).
    extension : str
        Extensión de los recortes a leer ('tif' o 'png').

    Returns
    -------
    list of str
        ROI ordenadas numéricamente, ej. ['ROI2', 'ROI8', 'ROI10'].
    """
    crops_dir = Path(files["output_dir"]) / crops_subfolder
    if not crops_dir.exists():
        print(f"No existe la carpeta de recortes: {crops_dir}")
        return []

    roi_numbers = []
    for crop_path in sorted(crops_dir.glob(f"*.{extension}")):
        match = re.search(r"_R(\d+)\.", crop_path.name)
        if match is None:
            print(f"Nombre inesperado, se omite: {crop_path.name}")
            continue
        roi_numbers.append(int(match.group(1)))

    roi_names = [f"ROI{n}" for n in sorted(set(roi_numbers))]
    print(f"ROIs encontradas en {crops_dir}: {len(roi_names)}")
    return roi_names
