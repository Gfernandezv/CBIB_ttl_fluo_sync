import pyabf
import numpy as np
import re
import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
import fnmatch
from pathlib import Path


def load_abf_temp_ttl(file_path, temp_channel=0, ttl_channel=1):
    """
    Carga temperatura y TTL desde un ABF.
    """
    abf = pyabf.ABF(file_path)

    abf.setSweep(0, channel=temp_channel)
    time = abf.sweepX.copy()
    temp = abf.sweepY.copy()

    abf.setSweep(0, channel=ttl_channel)
    ttl = abf.sweepY.copy()

    return {
        "time": time,
        "temperature": temp,
        "ttl": ttl,
        "metadata": {
            "file_path": abf.abfFilePath,
            "protocol": abf.protocol,
            "sample_rate_hz": abf.dataRate,
            "duration_sec": abf.sweepLengthSec,
            "channel_names": abf.adcNames,
            "channel_units": abf.adcUnits,
        }
    }



def mean_temperature_during_ttl(temperature, time, ttl_events):
    """
    Calcula la temperatura promedio dentro de cada pulso TTL.

    Parámetros
    ----------
    temperature : array-like
        Señal de temperatura.
    time : array-like
        Vector de tiempo en segundos.
    ttl_events : pandas.DataFrame
        DataFrame con columnas:
        - onset_time
        - offset_time
        - onset_idx
        - offset_idx

    Retorna
    -------
    result_df : pandas.DataFrame
        Copia de ttl_events con columnas añadidas:
        - temp_mean
        - temp_min
        - temp_max
        - n_points
    """
    temperature = np.asarray(temperature, dtype=float)
    time = np.asarray(time, dtype=float)

    results = []

    for _, row in ttl_events.iterrows():
        t0 = row["onset_time"]
        t1 = row["offset_time"]

        mask = (time >= t0) & (time <= t1)
        temp_segment = temperature[mask]

        if temp_segment.size == 0:
            temp_mean = np.nan
            temp_min = np.nan
            temp_max = np.nan
            n_points = 0
        else:
            temp_mean = np.nanmean(temp_segment)
            temp_min = np.nanmin(temp_segment)
            temp_max = np.nanmax(temp_segment)
            n_points = temp_segment.size

        results.append({
            "temp_mean": temp_mean,
            "temp_min": temp_min,
            "temp_max": temp_max,
            "n_points": n_points
        })

    results_df = pd.DataFrame(results)
    return pd.concat([ttl_events.reset_index(drop=True), results_df], axis=1)

def detect_ttl_events(ttl, time, threshold=None, min_duration_s=0):
    """
    Detecta inicios y finales de pulsos TTL.

    Parámetros
    ----------
    ttl : array-like
        Señal TTL.
    time : array-like
        Tiempo en segundos, mismo largo que ttl.
    threshold : float or None
        Umbral para binarizar la señal. Si es None, usa el punto medio
        entre el mínimo y el máximo de la señal.
    min_duration_s : float
        Duración mínima de pulso para conservarlo.

    Retorna
    -------
    events_df : pandas.DataFrame
        Tabla con:
        - onset_idx
        - offset_idx
        - onset_time
        - offset_time
        - duration_s
    ttl_binary : np.ndarray
        Señal binaria usada para la detección.
    """
    ttl = np.asarray(ttl).astype(float)
    time = np.asarray(time).astype(float)

    if ttl.shape != time.shape:
        raise ValueError("ttl y time deben tener la misma longitud.")

    if threshold is None:
        threshold = (np.nanmin(ttl) + np.nanmax(ttl)) / 2

    ttl_binary = ttl > threshold
    diff = np.diff(ttl_binary.astype(int))

    onset_idx = np.where(diff == 1)[0] + 1
    offset_idx = np.where(diff == -1)[0] + 1

    # si la señal parte en alto, agregar inicio en 0
    if ttl_binary[0]:
        onset_idx = np.r_[0, onset_idx]

    # si la señal termina en alto, agregar final en último índice
    if ttl_binary[-1]:
        offset_idx = np.r_[offset_idx, len(ttl_binary) - 1]

    if len(onset_idx) != len(offset_idx):
        raise ValueError(
            f"Número desigual de inicios ({len(onset_idx)}) y finales ({len(offset_idx)})."
        )

    onset_time = time[onset_idx]
    offset_time = time[offset_idx]
    duration_s = offset_time - onset_time

    events_df = pd.DataFrame({
        "onset_idx": onset_idx,
        "offset_idx": offset_idx,
        "onset_time": onset_time,
        "offset_time": offset_time,
        "duration_s": duration_s,
    })

    if min_duration_s > 0:
        events_df = events_df.loc[events_df["duration_s"] >= min_duration_s].reset_index(drop=True)

    return events_df, ttl_binary

def long_format(temp_ttl_df, fluo_file, start_idx=70, drop_frames_without_ttl=True):
    mean_temp = temp_ttl_df["temp_mean"]

    df = pd.read_csv(fluo_file)

    intden_cols = [c for c in df.columns if re.fullmatch(r"IntDen\d+", c)]
    csv = df[intden_cols].copy()

    csv["frame"] = range(1, len(csv) + 1)
    csv["TTL_index"] = range(start_idx, start_idx + len(csv))

    ttl_table = temp_ttl_df.copy()
    ttl_table["TTL_index"] = ttl_table.index

    merged = pd.merge(
        csv,
        ttl_table,
        on="TTL_index",
        how="left"
    )

    if drop_frames_without_ttl:
        missing_ttl = merged["temp_mean"].isna()
        n_missing = int(missing_ttl.sum())
        if n_missing > 0:
            first_missing = merged.loc[missing_ttl, "TTL_index"].min()
            last_missing = merged.loc[missing_ttl, "TTL_index"].max()
            print(
                f"Se eliminaron {n_missing} frames sin TTL "
                f"(TTL_index {first_missing} a {last_missing})."
            )
            merged = merged.loc[~missing_ttl].reset_index(drop=True)

    intden_cols = [c for c in merged.columns if re.fullmatch(r"IntDen\d+", c)]

    long_df = merged.melt(
        id_vars=[c for c in merged.columns if c not in intden_cols],
        value_vars=intden_cols,
        var_name="ROI",
        value_name="IntDen"
    )

    long_df["ROI"] = long_df["ROI"].str.replace("IntDen", "ROI", regex=False)

    return long_df

def Normalization(long_df, target_temp=20, tol=0.5, norm_type=2):
    df_norm = long_df.copy()

    df_norm["NormSignal"] = np.nan
    df_norm["I_baseline"] = np.nan
    df_norm["Imax"] = np.nan
    df_norm["T_at_Imax"] = np.nan

    for roi in df_norm["ROI"].unique():
        mask_roi = df_norm["ROI"] == roi
        sub = df_norm.loc[mask_roi].copy()

        baseline_mask = (
            (sub["temp_mean"] >= target_temp - tol) &
            (sub["temp_mean"] <= target_temp + tol)
        )

        baseline_vals = sub.loc[baseline_mask, "IntDen"]

        if len(baseline_vals) == 0:
            print(f"{roi}: sin datos cerca de {target_temp}°C")
            continue

        I_baseline = baseline_vals.mean()
        Imax = sub["IntDen"].max()
        T_at_Imax = sub.loc[sub["IntDen"].idxmax(), "temp_mean"]

        if norm_type == 1:
            norm = (sub["IntDen"] - I_baseline) / Imax
            msg = "Normalizacion por Maximo"
        elif norm_type == 2:
            norm = (sub["IntDen"] - I_baseline) / I_baseline
            msg = "Normalizacion por intensidad basal"
        else:
            msg = "Tipo no valido"
            return None

        df_norm.loc[mask_roi, "NormSignal"] = norm.values
        df_norm.loc[mask_roi, "I_baseline"] = I_baseline
        df_norm.loc[mask_roi, "Imax"] = Imax
        df_norm.loc[mask_roi, "T_at_Imax"] = T_at_Imax

    return df_norm


def graph_rois_before_filter(df_phase, phase=None, value_col="NormSignal", alpha=0.25):
    plot_df = df_phase.copy()

    if phase is not None:
        if "phase" not in plot_df.columns:
            raise ValueError("df_phase debe contener una columna 'phase' para filtrar por fase.")
        plot_df = plot_df[plot_df["phase"] == phase].copy()

    if plot_df.empty:
        print("Advertencia: no hay datos para graficar antes del filtro.")
        return plot_df

    required_cols = {"ROI", "temp_mean", value_col}
    missing = required_cols - set(plot_df.columns)
    if missing:
        raise ValueError(f"Faltan columnas en df_phase: {sorted(missing)}")

    plt.figure(figsize=(8, 5))

    for roi in plot_df["ROI"].dropna().unique():
        sub = plot_df[plot_df["ROI"] == roi].sort_values("temp_mean")
        plt.plot(
            sub["temp_mean"],
            sub[value_col],
            linewidth=0.8,
            alpha=alpha
        )

    plt.axhline(0, linestyle="--", alpha=0.4)
    plt.xlabel("Temperature (°C)")
    plt.ylabel(value_col)
    title = "ROI responses before filtering"
    if phase is not None:
        title += f" ({phase})"
    plt.title(title)
    plt.tight_layout()
    plt.show()

    print(f"ROIs graficadas antes del filtro: {plot_df['ROI'].nunique()}")
    print(f"Filas graficadas: {len(plot_df)}")

    return plot_df


def selected_rois_from_long(
    df,
    folder,
    folder_col="source_folder",
    status_col="ROI_status",
    included_status=1,
):
    """
    Devuelve la lista de ROI con status_col == included_status para una
    carpeta/muestra, leyendo directamente el ROI_status ya fijado en el
    *_preprocessed_long.csv (compilado en preprocessed_all).

    Pensado como alternativa a curar la lista a mano o vía recortes de
    imagen: si ROI_status ya refleja la decisión final de
    01_preprocessing.ipynb, esta función la reutiliza tal cual, por ejemplo
    para alimentar selected_rois en process_sample() o para graficar.

    folder acepta patrones estilo shell (fnmatch), ej. 'mut27_*' junta
    'mut27_image10' y 'mut27_image11'. OJO: el nombre de ROI se repite entre
    muestras, así que la lista resultante mezcla ROI de distintas carpetas
    bajo el mismo nombre — sirve para inspección/graficado, pero no como
    selected_rois de un único process_sample() (que opera sobre una sola
    carpeta/archivo crudo).

    Parameters
    ----------
    df : pandas.DataFrame
        Tabla larga con folder_col, ROI y status_col (ej. preprocessed_all).
    folder : str
        Valor (o patrón fnmatch) de folder_col a consultar, ej.
        'mut45_image13' o 'mut27_*'.
    folder_col : str
        Columna que identifica la carpeta/muestra de origen.
    status_col : str
        Columna de estatus de inclusión.
    included_status : int
        Valor de status_col que indica ROI incluida (por defecto 1).

    Returns
    -------
    list of str
        ROI incluidas, ordenadas numéricamente, ej. ['ROI2', 'ROI8', 'ROI10'].
    """
    required_cols = {folder_col, "ROI", status_col}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Faltan columnas en df: {sorted(missing)}")

    folder_mask = df[folder_col].astype(str).apply(lambda v: fnmatch.fnmatch(v, folder))
    matched_folders = sorted(df.loc[folder_mask, folder_col].dropna().unique())
    if len(matched_folders) > 1:
        print(
            f"'{folder}' coincide con {len(matched_folders)} carpetas: {matched_folders}. "
            "La lista de ROI resultante mezcla nombres repetidos entre muestras."
        )

    subset = df[folder_mask & (df[status_col] == included_status)]
    rois = subset["ROI"].dropna().astype(str).unique().tolist()
    rois.sort(key=lambda roi: int("".join(ch for ch in roi if ch.isdigit()) or 0))

    print(f"{folder}: {len(rois)} ROIs con {status_col}={included_status}")
    return rois


def graph_selected_rois_by_folder(
    df,
    folder=None,
    rois=None,
    folder_col="source_folder",
    status_col="ROI_status",
    included_status=1,
    value_col="NormSignal",
    x_col="temp_mean",
    alpha=0.8,
    show_mean_se=False,
    n_bins=20,
):
    """
    Grafica value_col vs x_col para un conjunto de ROI de una o varias muestras.

    Acepta dos fuentes de datos:

    - preprocessed_all (el compilado de todos los *_preprocessed_long.csv en
      02_processing.ipynb): pasa folder para filtrar por folder_col, y las
      ROI se seleccionan por status_col == included_status (o por rois, si
      se entrega explícitamente). folder acepta patrones fnmatch, ej.
      'mut27_*' junta 'mut27_image10' y 'mut27_image11'. Como el nombre de
      ROI se repite entre muestras, si el patrón matchea más de una carpeta
      cada línea se etiqueta "folder:ROI" en vez de solo "ROI".
    - processed["df_phase"] (la salida de ttl.process_sample()/reprocesos
      con selected_rois): no tiene folder_col ni hace falta filtrar por
      status_col, ya viene restringido a las ROI reprocesadas. En ese caso
      no pases folder (o pasa folder solo para el título del gráfico).

    Parameters
    ----------
    df : pandas.DataFrame
        Tabla larga con ROI, x_col, value_col (y folder_col/status_col si
        corresponde filtrar por ellos).
    folder : str, optional
        Valor de folder_col a graficar, ej. 'mut45_image13'. Si df no tiene
        folder_col (ej. viene de process_sample()), se usa solo para el
        título del gráfico.
    rois : list of str, optional
        Si se entrega, grafica exactamente estas ROI (ignora status_col).
        Útil para graficar selected_rois directamente, sin depender de que
        df tenga status_col.
    folder_col : str
        Columna que identifica la carpeta/muestra de origen, si existe en df.
    status_col : str
        Columna de estatus de inclusión, usada solo si rois es None y
        status_col existe en df.
    included_status : int
        Valor de status_col que indica ROI incluida (por defecto 1).
    value_col : str
        Columna a graficar en el eje Y (por defecto 'NormSignal').
    x_col : str
        Columna a graficar en el eje X (por defecto 'temp_mean').
    alpha : float
        Transparencia de cada línea de ROI.
    show_mean_se : bool
        Si True, agrega una curva de media ± SEM sobre las ROI graficadas
        (todas juntas, sin distinguir folder/ROI), calculada agrupando x_col
        en n_bins intervalos. Útil para ver la tendencia general en un batch
        de varias muestras (folder con patrón) sin leer línea por línea.
    n_bins : int
        Número de intervalos de x_col usados para la media ± SEM, si
        show_mean_se=True.

    Returns
    -------
    pandas.DataFrame
        Subconjunto graficado.
    """
    required_cols = {"ROI", x_col, value_col}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Faltan columnas en df: {sorted(missing)}")

    plot_df = df.copy()
    has_folder_col = folder_col in plot_df.columns
    multi_folder = False

    if folder is not None and has_folder_col:
        folder_mask = plot_df[folder_col].astype(str).apply(lambda v: fnmatch.fnmatch(v, folder))
        plot_df = plot_df[folder_mask]
        multi_folder = plot_df[folder_col].nunique() > 1

    if rois is not None:
        selected_set = {str(roi).strip() for roi in rois}
        plot_df = plot_df[plot_df["ROI"].astype(str).str.strip().isin(selected_set)]
    elif status_col in plot_df.columns:
        plot_df = plot_df[plot_df[status_col] == included_status]

    plot_df = plot_df.copy()
    if plot_df.empty:
        print(f"No hay ROIs para graficar (folder={folder!r}, rois={rois!r}).")
        return plot_df

    # Si folder es un patron que junta varias carpetas, el nombre de ROI se
    # repite entre ellas: agrupar solo por ROI mezclaria lineas de muestras
    # distintas bajo la misma etiqueta/color.
    group_cols = [folder_col, "ROI"] if multi_folder else ["ROI"]

    plt.figure(figsize=(8, 5))
    for key, sub in plot_df.groupby(group_cols):
        sub = sub.sort_values(x_col)
        label = ":".join(key) if multi_folder else key
        plt.plot(sub[x_col], sub[value_col], linewidth=1.2, alpha=alpha, label=label)

    if show_mean_se:
        bins = pd.cut(plot_df[x_col], bins=n_bins)
        bin_summary = plot_df.groupby(bins, observed=True)[value_col].agg(
            mean="mean",
            sem=lambda x: x.std(ddof=1) / (len(x) ** 0.5) if len(x) > 1 else 0.0,
            n="size",
        )
        bin_centers = bin_summary.index.map(lambda iv: iv.mid).astype(float)
        plt.errorbar(
            bin_centers,
            bin_summary["mean"],
            yerr=bin_summary["sem"],
            color="black",
            marker="o",
            linewidth=2,
            capsize=3,
            label="mean ± SEM",
        )

    plt.axhline(0, linestyle="--", alpha=0.4)
    plt.xlabel(x_col)
    plt.ylabel(value_col)
    title = f"{folder}: " if folder is not None else ""
    n_folders = plot_df[folder_col].nunique() if has_folder_col else 1
    folders_label = f", {n_folders} muestras" if multi_folder else ""
    plt.title(f"{title}{value_col} de ROI seleccionadas ({plot_df['ROI'].nunique()} ROIs{folders_label})")
    plt.legend(title="folder:ROI" if multi_folder else "ROI", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    plt.tight_layout()
    plt.show()

    print(f"{folder}: {plot_df['ROI'].nunique()} ROIs graficadas, {len(plot_df)} filas")

    return plot_df


def graph_selected_rois_by_phase(
    df_phase,
    selected,
    value_col="NormSignal",
    phases=None,
    phase=None,
):
    if phases is None:
        phases = {"heating": "tomato", "cooling": "royalblue"}

    if phase is not None:
        if phase not in phases:
            raise ValueError(f"phase debe ser una de: {list(phases)}")
        phases = {phase: phases[phase]}

    if selected.empty:
        print("Advertencia: selected está vacío. No hay ROIs para graficar.")
        return pd.DataFrame()

    required_cols = {"ROI", "temp_mean", "phase", value_col}
    missing = required_cols - set(df_phase.columns)
    if missing:
        raise ValueError(f"Faltan columnas en df_phase: {sorted(missing)}")

    plot_df = df_phase[df_phase["ROI"].isin(selected["ROI"])].copy()
    if plot_df.empty:
        print("Advertencia: no hay filas para las ROIs seleccionadas.")
        return plot_df

    for roi in plot_df["ROI"].dropna().unique():
        sub = plot_df[plot_df["ROI"] == roi]

        plt.figure(figsize=(5, 4))
        for phase, color in phases.items():
            part = sub[sub["phase"] == phase]
            if part.empty:
                continue
            plt.scatter(
                part["temp_mean"],
                part[value_col],
                s=20,
                alpha=0.7,
                label=phase,
                color=color,
            )

        plt.xlabel("Temperatura promedio durante TTL/Frame")
        plt.ylabel(value_col)
        plt.title(f"{roi}: {value_col} vs Temperatura")
        plt.legend()
        plt.tight_layout()
        plt.show()

    return plot_df


def detect_roi_signal_outliers(
    df,
    value_col="NormSignal",
    group_col="ROI",
    order_col="TTL_index",
    window=11,
    z_thresh=4.0,
    min_abs_residual=0.05,
):
    """
    Detecta saltos/puntos anómalos respecto a la tendencia local de cada ROI.

    El método usa una mediana móvil por ROI como tendencia robusta. Luego calcula
    el residual punto - tendencia y lo escala usando MAD (median absolute
    deviation). Un punto se marca como outlier si:

    - abs(robust_z) >= z_thresh
    - abs(residual) >= min_abs_residual

    No elimina datos: devuelve una copia con columnas nuevas:
    - local_trend
    - residual
    - robust_z
    - is_outlier
    """
    if window < 3:
        raise ValueError("window debe ser >= 3.")
    if window % 2 == 0:
        window += 1

    required_cols = {group_col, order_col, value_col}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Faltan columnas en df: {sorted(missing)}")

    result = df.copy()
    result["local_trend"] = np.nan
    result["residual"] = np.nan
    result["robust_z"] = np.nan
    result["is_outlier"] = False

    for _, idx in result.groupby(group_col).groups.items():
        sub = result.loc[idx].sort_values(order_col)
        values = sub[value_col].astype(float)

        trend = values.rolling(
            window=window,
            center=True,
            min_periods=max(3, window // 2)
        ).median()
        trend = trend.bfill().ffill()

        residual = values - trend
        mad = np.nanmedian(np.abs(residual - np.nanmedian(residual)))
        scale = 1.4826 * mad

        if not np.isfinite(scale) or scale == 0:
            scale = np.nanstd(residual)
        if not np.isfinite(scale) or scale == 0:
            robust_z = pd.Series(0.0, index=sub.index)
        else:
            robust_z = residual / scale

        is_outlier = (
            robust_z.abs().ge(z_thresh) &
            residual.abs().ge(min_abs_residual)
        )

        result.loc[sub.index, "local_trend"] = trend.values
        result.loc[sub.index, "residual"] = residual.values
        result.loc[sub.index, "robust_z"] = robust_z.values
        result.loc[sub.index, "is_outlier"] = is_outlier.values

    print(f"Outliers detectados: {int(result['is_outlier'].sum())} / {len(result)} puntos")
    if group_col in result.columns:
        n_rois = result.loc[result["is_outlier"], group_col].nunique()
        print(f"ROIs con al menos un outlier: {n_rois}")

    return result


def graph_roi_outlier_detection(
    df_outliers,
    roi,
    value_col="NormSignal",
    x_col="temp_mean",
):
    sub = df_outliers[df_outliers["ROI"] == roi].copy()
    if sub.empty:
        print(f"No hay datos para {roi}.")
        return sub

    normal = sub[~sub["is_outlier"]]
    outliers = sub[sub["is_outlier"]]

    plt.figure(figsize=(6, 4))
    plt.scatter(normal[x_col], normal[value_col], s=20, alpha=0.55, label="normal")
    if not outliers.empty:
        plt.scatter(
            outliers[x_col],
            outliers[value_col],
            s=45,
            alpha=0.9,
            color="crimson",
            label="outlier"
        )

    if "local_trend" in sub.columns and x_col == "TTL_index":
        ordered = sub.sort_values("TTL_index")
        plt.plot(ordered[x_col], ordered["local_trend"], color="black", linewidth=1.5, label="tendencia")

    plt.xlabel(x_col)
    plt.ylabel(value_col)
    plt.title(f"{roi}: detección de saltos")
    plt.legend()
    plt.tight_layout()
    plt.show()

    return sub


def filter_analysis_window(
    df_phase,
    ttl_range=None,
    time_range=None,
):
    """
    Filtra una ventana de análisis sin cambiar la alineación frame-TTL.

    Usa esta función después de process_sample() para excluir segmentos que no
    quieres analizar. Esto conserva start_ttl como parámetro de alineación y
    evita desplazar los frames.

    Parameters
    ----------
    df_phase : pandas.DataFrame
        Tabla procesada con columnas TTL_index, onset_time, offset_time y phase.
    ttl_range : tuple or None
        Rango inclusivo de TTL_index, por ejemplo (11, 350).
    time_range : tuple or None
        Rango de tiempo en segundos, por ejemplo (30, 700). Conserva frames cuyo
        onset_time cae dentro de ese rango.

    Returns
    -------
    pandas.DataFrame
        Subconjunto de df_phase dentro de la ventana elegida.
    """
    if ttl_range is not None and time_range is not None:
        raise ValueError("Usa ttl_range o time_range, no ambos a la vez.")

    df = df_phase.copy()

    if ttl_range is not None:
        ttl_start, ttl_end = ttl_range
        df = df[(df["TTL_index"] >= ttl_start) & (df["TTL_index"] <= ttl_end)].copy()
        label = f"TTL_index {ttl_start} a {ttl_end}"
    elif time_range is not None:
        t_start, t_end = time_range
        if "onset_time" not in df.columns:
            raise ValueError("df_phase debe contener 'onset_time' para filtrar por tiempo.")
        df = df[(df["onset_time"] >= t_start) & (df["onset_time"] <= t_end)].copy()
        label = f"tiempo {t_start} a {t_end} s"
    else:
        label = "sin filtro"

    print(f"Ventana de análisis: {label}")
    print(f"  frames únicos: {df['frame'].nunique() if 'frame' in df.columns else 'NA'}")
    print(f"  TTL únicos: {df['TTL_index'].nunique()}")
    print(f"  ROIs: {df['ROI'].nunique()}")
    if "phase" in df.columns and not df.empty:
        print("  fases:")
        print(df[["TTL_index", "phase"]].drop_duplicates()["phase"].value_counts())

    return df


def plot_temp_ttl(
    data,
    temp_ttl_df=None,
    frame_counts=None,
    only_used_frames=False,
):
    time = data["time"]
    temperature = data["temperature"]
    ttl = data["ttl"]

    if only_used_frames:
        if temp_ttl_df is None or frame_counts is None:
            raise ValueError("temp_ttl_df y frame_counts son necesarios si only_used_frames=True.")

        first_ttl = frame_counts["assigned_ttl_first"]
        last_ttl = min(
            frame_counts["assigned_ttl_last"],
            frame_counts["detected_ttl_last"]
        )

        used_ttl = temp_ttl_df.loc[first_ttl:last_ttl]
        if used_ttl.empty:
            raise ValueError("No hay TTLs usados para graficar.")

        t0 = used_ttl["onset_time"].min()
        t1 = used_ttl["offset_time"].max()
        mask = (time >= t0) & (time <= t1)
        title_suffix = " (frames con TTL asociado)"
    else:
        mask = np.ones_like(time, dtype=bool)
        title_suffix = " (registro completo)"

    fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True)

    axes[0].plot(time[mask], temperature[mask], lw=0.8)
    axes[0].set_ylabel("Temperature")
    axes[0].set_title(f"Temperature channel{title_suffix}")

    axes[1].plot(time[mask], ttl[mask], lw=0.8)
    axes[1].set_ylabel("TTL")
    axes[1].set_xlabel("Time (s)")
    axes[1].set_title(f"TTL channel{title_suffix}")

    plt.tight_layout()
    plt.show()

    return fig, axes


def summarize_rois_by_temp_ranges(
    df_phase,
    phase=None,
    temp_ranges=None,
    value_col="NormSignal",
    min_points=3,
    min_abs_delta=0.01,
    verbose=True,
):
    if temp_ranges is None:
        temp_ranges = {
            "low": (22, 28),
            "mid": (32, 35),
            "high": (36, 42),
        }

    df = df_phase.copy()
    if phase is not None:
        if "phase" not in df.columns:
            raise ValueError("df_phase debe contener una columna 'phase' para filtrar por fase.")
        df = df[df["phase"] == phase].copy()

    required_cols = {"ROI", "temp_mean", value_col}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Faltan columnas en df_phase: {sorted(missing)}")

    rows = []
    for roi in df["ROI"].dropna().unique():
        sub = df[df["ROI"] == roi]
        row = {"ROI": roi}

        enough_points = True
        for name, (t_min, t_max) in temp_ranges.items():
            in_range = sub[(sub["temp_mean"] >= t_min) & (sub["temp_mean"] <= t_max)]
            row[f"{name}_mean"] = in_range[value_col].mean()
            row[f"{name}_sd"] = in_range[value_col].std()
            row[f"{name}_n"] = len(in_range)
            if len(in_range) < min_points:
                enough_points = False

        if {"low_mean", "high_mean"} <= set(row):
            row["delta_high_low"] = row["high_mean"] - row["low_mean"]
        else:
            first_name = next(iter(temp_ranges))
            last_name = next(reversed(temp_ranges))
            row["delta_high_low"] = row[f"{last_name}_mean"] - row[f"{first_name}_mean"]

        if not enough_points or pd.isna(row["delta_high_low"]):
            row["trend"] = "insufficient"
        elif row["delta_high_low"] > min_abs_delta:
            row["trend"] = "increase"
        elif row["delta_high_low"] < -min_abs_delta:
            row["trend"] = "decrease"
        else:
            row["trend"] = "stable"

        rows.append(row)

    summary = pd.DataFrame(rows)
    if summary.empty:
        print("No hay ROIs para resumir.")
        return summary

    if verbose:
        print("Resumen por tendencia:")
        print(summary["trend"].value_counts(dropna=False))

    return summary


def graph_temp_range_summary(
    roi_summary,
    temp_ranges=None,
    value_prefixes=None,
):
    if roi_summary.empty:
        print("Advertencia: roi_summary está vacío. No hay datos para graficar.")
        return pd.DataFrame()

    if temp_ranges is None:
        temp_ranges = {
            "low": (22, 28),
            "mid": (32, 35),
            "high": (36, 42),
        }

    if value_prefixes is None:
        value_prefixes = list(temp_ranges.keys())

    mean_cols = [f"{name}_mean" for name in value_prefixes]
    missing = [col for col in mean_cols if col not in roi_summary.columns]
    if missing:
        raise ValueError(f"Faltan columnas en roi_summary: {missing}")

    x = np.arange(len(value_prefixes))
    labels = [
        f"{name}\n{temp_ranges[name][0]}-{temp_ranges[name][1]}°C"
        if name in temp_ranges else name
        for name in value_prefixes
    ]

    values = roi_summary[mean_cols].to_numpy(dtype=float)
    group_mean = np.nanmean(values, axis=0)
    group_sem = np.nanstd(values, axis=0, ddof=1) / np.sqrt(np.sum(~np.isnan(values), axis=0))

    plt.figure(figsize=(7, 5))
    for _, row in roi_summary.iterrows():
        plt.plot(x, row[mean_cols].to_numpy(dtype=float), color="0.75", linewidth=0.8, alpha=0.5)

    plt.errorbar(
        x,
        group_mean,
        yerr=group_sem,
        color="black",
        linewidth=2.5,
        marker="o",
        capsize=4,
    )
    plt.axhline(0, linestyle="--", alpha=0.4)
    plt.xticks(x, labels)
    plt.ylabel("Mean normalized signal")
    plt.title("Average ROI response by temperature range")
    plt.tight_layout()
    plt.show()

    group_summary = pd.DataFrame({
        "range": value_prefixes,
        "mean": group_mean,
        "sem": group_sem,
        "n_rois": np.sum(~np.isnan(values), axis=0),
    })

    return group_summary


def filter_rois_by_trend(
    roi_summary,
    trend="increase",
    min_abs_delta=0.01,
    min_points=3,
    verbose=True,
):
    if trend not in {"increase", "decrease", "stable", "insufficient"}:
        raise ValueError("trend debe ser 'increase', 'decrease', 'stable' o 'insufficient'.")

    if roi_summary.empty:
        print("Advertencia: roi_summary está vacío. No hay ROIs para filtrar.")
        return roi_summary.copy()

    required_cols = {"ROI", "delta_high_low"}
    missing = required_cols - set(roi_summary.columns)
    if missing:
        raise ValueError(f"Faltan columnas en roi_summary: {sorted(missing)}")

    n_cols = [col for col in roi_summary.columns if col.endswith("_n")]
    if n_cols:
        enough_points = (roi_summary[n_cols] >= min_points).all(axis=1)
    else:
        enough_points = True

    if "trend" in roi_summary.columns:
        trend_mask = roi_summary["trend"] == trend
    elif trend == "increase":
        trend_mask = roi_summary["delta_high_low"] > min_abs_delta
    elif trend == "decrease":
        trend_mask = roi_summary["delta_high_low"] < -min_abs_delta
    elif trend == "stable":
        trend_mask = roi_summary["delta_high_low"].abs() <= min_abs_delta
    else:
        trend_mask = ~enough_points | roi_summary["delta_high_low"].isna()

    if trend == "insufficient":
        selected = roi_summary[trend_mask].copy()
    else:
        selected = roi_summary[enough_points & trend_mask].copy()

    if verbose:
        print(
            f"ROIs seleccionadas con tendencia {trend}: "
            f"{len(selected)} / {len(roi_summary)}"
        )

    return selected


def add_temperature_phase(df):
    temp_phase = (
        df[["TTL_index", "temp_mean"]]
        .drop_duplicates()
        .sort_values("TTL_index")
        .copy()
    )

    temp_phase["dT"] = temp_phase["temp_mean"].diff()
    temp_phase["phase"] = np.where(temp_phase["dT"] >= 0, "heating", "cooling")
    temp_phase.loc[temp_phase["dT"].isna(), "phase"] = "heating"

    return df.drop(columns=["phase", "dT"], errors="ignore").merge(
        temp_phase[["TTL_index", "dT", "phase"]],
        on="TTL_index",
        how="left"
    )


def process_sample(
    input_dir,
    start_ttl=1,
    threshold=2.0,
    target_temp=25,
    norm_type=1,
    temp_channel=0,
    ttl_channel=1,
    drop_frames_without_ttl=True,
    phase_filter=None,
    temp_range=None,
    selected_rois="all",
):
    """
    Procesa una muestra de imagen + ABF manteniendo la alineación frame-TTL.

    Cada TTL marca el inicio de un frame. Si la cámara comienza con el primer
    TTL generado por el sistema, usa start_ttl=0:

        frame 1 -> TTL 0
        frame 2 -> TTL 1
        ...

    start_ttl corrige la alineación entre el CSV de ImageJ y los TTL detectados.
    No debe usarse para eliminar segmentos iniciales si la cámara realmente
    comenzó en TTL 0, porque eso desalinearía imagen, temperatura y TTL.

    El orden del procesamiento es:

        formato largo -> detección heating/cooling -> filtro opcional de fase
        -> normalización

    Si phase_filter es "cooling" o "heating", la normalización se calcula solo
    con los puntos de esa fase. Por lo tanto, I_baseline, Imax, T_at_Imax y
    NormSignal describen el subconjunto filtrado, no el registro completo.

    Si temp_range=(t_min, t_max), la normalización se calcula solo con puntos
    cuya temp_mean cae dentro de ese rango (se aplica después de phase_filter).
    Esto evita que Imax quede definido por temperaturas fuera del rango que
    después se usa en analyze_temp_trends / temp_ranges. target_temp debe caer
    dentro de temp_range, o I_baseline quedará sin datos.

    Si selected_rois es una lista de nombres de ROI (ej. ["ROI2", "ROI8"]),
    la normalización se calcula solo con esas ROIs (se aplica después de
    phase_filter y temp_range). Si es "all" (por defecto), se usan todas las
    ROIs. Como Normalization() calcula I_baseline/Imax por ROI de forma
    independiente, restringir selected_rois no cambia el NormSignal de las
    ROIs incluidas respecto de correr con "all" y filtrar después; sirve para
    reprocesar/verificar un subconjunto ya curado (ej. roi_valid_from_crops en
    01_preprocessing.ipynb) sin volver a normalizar ROIs descartadas.

    Para excluir segmentos no deseados del análisis, filtra después de procesar:

        df_phase_clean = df_phase[df_phase["TTL_index"] >= 11].copy()

    Si además necesitas solo cooling:

        cooling_df_clean = df_phase_clean[df_phase_clean["phase"] == "cooling"].copy()

    En ese caso, cooling_df_clean contiene solo puntos que:
    - tienen TTL_index dentro del rango elegido
    - pertenecen a la fase cooling

    Si hay más frames en el CSV que TTL detectados, los frames sin TTL asociado
    se eliminan cuando drop_frames_without_ttl=True.

    Returns
    -------
    dict
        Diccionario con los datos intermedios y finales:
        - df_long: formato largo antes de clasificar fase.
        - df_long_phase: formato largo con dT y phase.
        - df_for_norm: datos usados para normalizar después de phase_filter y temp_range.
        - df_norm / df_phase: datos normalizados, conservando phase.
        - heating_df / cooling_df: subconjuntos normalizados por fase.
        - frame_counts: resumen de frames/TTL usados.
    """
    if phase_filter not in {None, "cooling", "heating"}:
        raise ValueError("phase_filter debe ser None, 'cooling' o 'heating'.")

    if temp_range is not None:
        t_min, t_max = temp_range
        if t_min >= t_max:
            raise ValueError("temp_range debe ser (t_min, t_max) con t_min < t_max.")

    files = find_files(input_dir)
    roi_csv = pd.read_csv(files["roi_csv"])
    frames_csv_total = len(roi_csv)
    assigned_ttl_first = start_ttl
    assigned_ttl_last = start_ttl + frames_csv_total - 1

    data = load_abf_temp_ttl(
        files["abf_file"],
        temp_channel=temp_channel,
        ttl_channel=ttl_channel
    )

    events_df, ttl_binary = detect_ttl_events(
        data["ttl"],
        data["time"],
        threshold=threshold
    )

    temp_ttl_df = mean_temperature_during_ttl(
        data["temperature"],
        data["time"],
        events_df
    )

    df_long = long_format(
        temp_ttl_df,
        files["roi_csv"],
        start_idx=start_ttl,
        drop_frames_without_ttl=drop_frames_without_ttl
    )

    frames_with_ttl = df_long["frame"].nunique()
    frames_without_ttl = frames_csv_total - frames_with_ttl

    df_long_phase = add_temperature_phase(df_long)

    if phase_filter is None:
        df_for_norm = df_long_phase.copy()
    else:
        df_for_norm = df_long_phase[df_long_phase["phase"] == phase_filter].copy()

    frames_after_phase_filter = df_for_norm["frame"].nunique()

    if temp_range is not None:
        t_min, t_max = temp_range
        df_for_norm = df_for_norm[
            (df_for_norm["temp_mean"] >= t_min) & (df_for_norm["temp_mean"] <= t_max)
        ].copy()
        if df_for_norm["frame"].nunique() == 0:
            raise ValueError(
                f"temp_range={temp_range} no dejó ningún frame "
                f"(había {frames_after_phase_filter} antes del filtro)."
            )

    if selected_rois != "all":
        selected_set = {str(roi).strip() for roi in selected_rois}
        df_for_norm = df_for_norm[
            df_for_norm["ROI"].astype(str).str.strip().isin(selected_set)
        ].copy()
        if df_for_norm.empty:
            raise ValueError(
                "Ninguna de las selected_rois quedó en df_for_norm; revisa los "
                "nombres de ROI o los filtros de phase_filter/temp_range."
            )

    df_norm = Normalization(
        df_for_norm,
        norm_type=norm_type,
        target_temp=target_temp
    )

    df_phase = df_norm.copy()

    outputs = {
        "files": files,
        "data": data,
        "events_df": events_df,
        "ttl_binary": ttl_binary,
        "temp_ttl_df": temp_ttl_df,
        "df_long": df_long,
        "df_long_phase": df_long_phase,
        "df_for_norm": df_for_norm,
        "df_norm": df_norm,
        "df_phase": df_phase,
        "heating_df": df_phase[df_phase["phase"] == "heating"].copy(),
        "cooling_df": df_phase[df_phase["phase"] == "cooling"].copy(),
        "phase_filter": phase_filter,
        "temp_range": temp_range,
        "selected_rois": selected_rois,
        "frame_counts": {
            "frames_csv_total": frames_csv_total,
            "frames_with_ttl": frames_with_ttl,
            "frames_without_ttl": frames_without_ttl,
            "frames_after_phase_filter": frames_after_phase_filter,
            "frames_after_temp_range": df_for_norm["frame"].nunique(),
            "assigned_ttl_first": assigned_ttl_first,
            "assigned_ttl_last": assigned_ttl_last,
            "detected_ttl_first": int(temp_ttl_df.index.min()),
            "detected_ttl_last": int(temp_ttl_df.index.max()),
        },
    }

    print("Procesamiento completo:")
    print(f"  eventos TTL: {len(events_df)}")
    print(f"  frames CSV total: {frames_csv_total}")
    print(f"  frames con TTL asociado: {frames_with_ttl}")
    print(f"  frames sin TTL asociado: {frames_without_ttl}")
    print(f"  filtro de fase previo a normalizar: {phase_filter or 'ninguno'}")
    print(f"  filtro de rango de temperatura previo a normalizar: {temp_range or 'ninguno'}")
    print(f"  frames usados para normalizar: {df_for_norm['frame'].nunique()}")
    print(f"  TTL_index asignados a frames: {assigned_ttl_first} a {assigned_ttl_last}")
    print(f"  TTL_index detectados: {temp_ttl_df.index.min()} a {temp_ttl_df.index.max()}")
    print(f"  filas df_norm: {len(df_norm)}")
    print(f"  ROIs: {df_norm['ROI'].nunique()}")
    print("  fases:")
    print(df_phase[["TTL_index", "phase"]].drop_duplicates()["phase"].value_counts())

    # Añadir información temporal: tiempo total y tiempo entre frames
    total_time = None
    try:
        total_time = data.get("metadata", {}).get("duration_sec", None)
    except Exception:
        total_time = None

    if total_time is None or (isinstance(total_time, float) and np.isnan(total_time)):
        try:
            times = np.asarray(data.get("time", []), dtype=float)
            if times.size:
                total_time = float(np.nanmax(times) - np.nanmin(times))
            else:
                total_time = float('nan')
        except Exception:
            total_time = float('nan')

    # tiempo entre frames: diferencias de onset_time en temp_ttl_df
    frame_intervals = np.array([])
    try:
        if "onset_time" in temp_ttl_df.columns and temp_ttl_df["onset_time"].notna().sum() > 1:
            frame_intervals = np.diff(temp_ttl_df["onset_time"].values.astype(float))
    except Exception:
        frame_intervals = np.array([])

    # imprimir resultados
    try:
        if np.isfinite(total_time):
            print(f"  tiempo total (s): {total_time:.2f}")
        else:
            print("  tiempo total (s): N/A")

        if frame_intervals.size > 0:
            print(f"  tiempo entre frames (mediana, s): {float(np.nanmedian(frame_intervals)):.4f}")
            print(f"  tiempo entre frames (media, s): {float(np.nanmean(frame_intervals)):.4f}")
        else:
            print("  tiempo entre frames: N/A")
    except Exception:
        pass

    return outputs


def process_selected_folders(base_dir, folders, preprocessed_all, **process_sample_kwargs):
    """
    Reprocesa varias muestras desde los archivos crudos y junta los resultados.

    process_sample() solo puede leer un ABF/CSV crudo a la vez, así que un
    folder con patrón (ej. 'mut27_*') no le sirve directamente: esta función
    recorre folders con un for, y para cada carpeta deriva su propia lista de
    ROI activas con selected_rois_from_long(preprocessed_all, folder=<esa
    carpeta>) antes de llamar a process_sample(). Como cada carpeta se
    resuelve por separado, no hay riesgo de mezclar ROI del mismo nombre
    entre muestras.

    Parameters
    ----------
    base_dir : str or Path
        Carpeta que contiene una subcarpeta por muestra (ej. Proc_data).
    folders : list of str
        Carpetas a reprocesar, ej. ['mut27_image10', 'mut27_image11']. Se
        puede construir con fnmatch, ej.:
            fnmatch.filter(preprocessed_all['source_folder'].unique(), 'mut27_*')
    preprocessed_all : pandas.DataFrame
        Compilado de *_preprocessed_long.csv (ver load_all_preprocessed_long),
        usado para derivar selected_rois de cada carpeta vía ROI_status.
    **process_sample_kwargs
        Argumentos comunes para process_sample() (start_ttl, threshold,
        target_temp, norm_type, phase_filter, temp_range, etc.), aplicados
        igual a todas las carpetas.

    Returns
    -------
    dict
        - "df_phase_all": concat de df_phase de cada carpeta reprocesada, con
          columna source_folder agregada (compatible con folder="patron" en
          graph_selected_rois_by_folder).
        - "processed_by_folder": {folder: dict de process_sample()} por si se
          necesita algo más además de df_phase.
        - "skipped": lista de (folder, motivo) para carpetas sin ROI activas
          o donde process_sample() falló.
    """
    df_phase_list = []
    processed_by_folder = {}
    skipped = []

    for folder in folders:
        selected_rois = selected_rois_from_long(preprocessed_all, folder=folder)
        if not selected_rois:
            skipped.append((folder, "sin ROIs con ROI_status=1"))
            continue

        input_dir = Path(base_dir) / folder
        try:
            processed = process_sample(input_dir, selected_rois=selected_rois, **process_sample_kwargs)
        except Exception as exc:
            skipped.append((folder, str(exc)))
            continue

        df_phase = processed["df_phase"].copy()
        df_phase["source_folder"] = folder
        df_phase_list.append(df_phase)
        processed_by_folder[folder] = processed

    df_phase_all = pd.concat(df_phase_list, ignore_index=True) if df_phase_list else pd.DataFrame()

    print(f"Reprocesadas {len(processed_by_folder)}/{len(folders)} carpetas.")
    if skipped:
        print("Carpetas omitidas:")
        for folder, reason in skipped:
            print(f"  {folder}: {reason}")

    return {
        "df_phase_all": df_phase_all,
        "processed_by_folder": processed_by_folder,
        "skipped": skipped,
    }


def analyze_temp_trends(
    df_phase,
    phase="cooling",
    temp_ranges=None,
    value_col="NormSignal",
    min_abs_delta=0.01,
    min_points=3,
    verbose=True,
):
    """
    Resume y clasifica el comportamiento de cada ROI por rangos de temperatura.

    Para cada ROI calcula el promedio de value_col dentro de cada rango definido
    en temp_ranges. La tendencia se define con delta_high_low:

        delta_high_low = high_mean - low_mean

    Por lo tanto:
    - increase: señal mayor en el rango high que en low
    - decrease: señal menor en el rango high que en low
    - stable: cambio absoluto menor o igual a min_abs_delta
    - insufficient: faltan puntos en algún rango o delta_high_low es NaN

    Si phase no es None, primero filtra df_phase por esa fase. Con phase="cooling",
    la clasificación describe cómo cambia la señal entre rangos de temperatura
    dentro de los puntos cooling, no necesariamente el orden temporal de la curva.
    """
    roi_summary = summarize_rois_by_temp_ranges(
        df_phase,
        phase=phase,
        temp_ranges=temp_ranges,
        value_col=value_col,
        min_points=min_points,
        min_abs_delta=min_abs_delta,
        verbose=verbose,
    )

    selected_increase = filter_rois_by_trend(
        roi_summary,
        trend="increase",
        min_abs_delta=min_abs_delta,
        min_points=min_points,
        verbose=verbose,
    )

    selected_decrease = filter_rois_by_trend(
        roi_summary,
        trend="decrease",
        min_abs_delta=min_abs_delta,
        min_points=min_points,
        verbose=verbose,
    )

    selected_stable = filter_rois_by_trend(
        roi_summary,
        trend="stable",
        min_abs_delta=min_abs_delta,
        min_points=min_points,
        verbose=verbose,
    )

    selected_insufficient = filter_rois_by_trend(
        roi_summary,
        trend="insufficient",
        min_abs_delta=min_abs_delta,
        min_points=min_points,
        verbose=verbose,
    )

    return {
        "roi_summary": roi_summary,
        "selected_increase": selected_increase,
        "selected_decrease": selected_decrease,
        "selected_stable": selected_stable,
        "selected_insufficient": selected_insufficient,
    }

import os
import pandas as pd

def export_results(
    df_norm,
    selected,
    sample="sample_07",
    genotype="WT",
    output_dir="outputs",
    export_filtered=True,
    roi_summary=None,
):
    """
    Exporta los datos normalizados completos y, opcionalmente,
    una versión filtrada con las ROI seleccionadas.

    processed contiene todas las filas presentes en df_norm. Si df_norm viene
    de process_sample(phase_filter="cooling"), processed contiene todas las
    ROI durante cooling, ya normalizadas sobre cooling.

    filtered contiene un subconjunto de processed: solo las ROI presentes en
    selected["ROI"]. Por ejemplo, si selected es selected_increase, filtered
    contiene solo las ROI clasificadas como increase.

    Si se entrega roi_summary, agrega las columnas de comportamiento por ROI
    (por ejemplo trend y delta_high_low) a processed y filtered.

    Parameters
    ----------
    df_norm : pandas.DataFrame
        DataFrame normalizado, idealmente salida de Normalization().
    selected : pandas.DataFrame
        DataFrame con ROI seleccionadas, idealmente salida de filter_rois().
        Debe contener una columna 'ROI'.
    sample : str
        Identificador de la muestra. Se agrega como columna si no existe en
        df_norm.
    genotype : str
        Genotipo o condición experimental. Se agrega como columna si no existe
        en df_norm.
    output_dir : str
        Carpeta de salida.
    export_filtered : bool
        Si True, exporta además el subconjunto filtrado.
    roi_summary : pandas.DataFrame, optional
        Resumen por ROI, idealmente salida de analyze_temp_trends()["roi_summary"].
        Debe contener la columna 'ROI'. Sus columnas se agregan a cada fila
        del export según la ROI correspondiente.
    
    Returns
    -------
    dict
        Rutas de los archivos exportados.
    """
    os.makedirs(output_dir, exist_ok=True)

    # columnas mínimas esperadas
    required_cols = {"ROI", "TTL_index", "temp_mean", "IntDen", "NormSignal"}
    missing = required_cols - set(df_norm.columns)
    if missing:
        raise ValueError(f"Faltan columnas en df_norm: {sorted(missing)}")

    # copia completa
    export_df = df_norm.copy()
    if "sample" not in export_df.columns:
        export_df["sample"] = sample
    if "genotype" not in export_df.columns:
        export_df["genotype"] = genotype

    if roi_summary is not None:
        if "ROI" not in roi_summary.columns:
            raise ValueError("roi_summary debe contener una columna 'ROI'.")

        summary_cols = [
            c for c in roi_summary.columns
            if c == "ROI" or c not in export_df.columns
        ]
        export_df = export_df.merge(
            roi_summary[summary_cols].drop_duplicates(subset=["ROI"]),
            on="ROI",
            how="left",
        )

    # columnas a exportar
    preferred_cols = [
        "sample",
        "genotype",
        "frame",
        "ROI",
        "TTL_index",
        "temp_mean",
        "IntDen",
        "NormSignal",
        "I_baseline",
        "Imax",
        "T_at_Imax",
        "onset_time",
        "offset_time",
        "duration_s",
        "temp_min",
        "temp_max",
        "n_points",
        "low_mean",
        "low_sd",
        "low_n",
        "mid_mean",
        "mid_sd",
        "mid_n",
        "high_mean",
        "high_sd",
        "high_n",
        "delta_high_low",
        "trend",
    ]

    cols_to_export = [c for c in preferred_cols if c in export_df.columns]
    export_df = export_df[cols_to_export]

    file_all = os.path.join(output_dir, f"{sample}_{genotype}_processed.csv")
    export_df.to_csv(file_all, index=False)

    output_files = {"processed": file_all}

    # subconjunto filtrado
    if export_filtered:
        if "ROI" not in selected.columns:
            raise ValueError("selected debe contener una columna 'ROI'.")

        filtered_rois = selected["ROI"].dropna().unique().tolist()

        filtered_export = export_df[export_df["ROI"].isin(filtered_rois)].copy()
        file_filtered = os.path.join(output_dir, f"{sample}_{genotype}_filtered.csv")
        filtered_export.to_csv(file_filtered, index=False)

        output_files["filtered"] = file_filtered

    print("Exportado:")
    for k, v in output_files.items():
        print(f"  {k}: {v}")

    return output_files


def load_all_preprocessed_long(base_dir, pattern="*/*_preprocessed_long.csv"):
    """
    Concatena todos los *_preprocessed_long.csv encontrados bajo base_dir,
    sin depender de un Excel de experimentos activos.

    Cada CSV ya trae 'sample' y 'genotype' nativamente (fijados en
    01_preprocessing.ipynb), así que no hace falta metadata externa para
    identificar una fila. Solo se agregan 'source_folder' y 'source_file'
    a partir de la ruta del archivo, para poder rastrear el origen y para
    detectar duplicados (más de un archivo por carpeta).

    Parameters
    ----------
    base_dir : str or Path
        Carpeta que contiene una subcarpeta por muestra (ej. Proc_data).
    pattern : str
        Patrón glob relativo a base_dir para encontrar los CSV.

    Returns
    -------
    dict
        - "preprocessed_all": DataFrame concatenado (vacío si no hay archivos).
        - "load_status": DataFrame con una fila por carpeta encontrada,
          listando sus *_preprocessed_long.csv y si hay más de uno.
    """
    base_dir = Path(base_dir)

    dfs = []
    files_by_folder = {}
    for file_path in sorted(base_dir.glob(pattern)):
        folder = file_path.parent.name
        files_by_folder.setdefault(folder, []).append(file_path.name)

        df = pd.read_csv(file_path)
        df["source_folder"] = folder
        df["source_file"] = file_path.name
        dfs.append(df)

    load_status = pd.DataFrame([
        {
            "folder": folder,
            "preprocessed_files": names,
            "n_files": len(names),
            "status": "ok" if len(names) == 1 else ("multiple" if len(names) > 1 else "faltan csv"),
        }
        for folder, names in sorted(files_by_folder.items())
    ])

    return {
        "preprocessed_all": pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame(),
        "load_status": load_status,
    }


def mark_roi_status(
    df,
    roi_exclude=None,
    roi_col="ROI",
    scoped_exclude=None,
    status_col="ROI_status",
    default_status=1,
    excluded_status=0,
    preserve_existing=True,
):
    """
    Marca ROIs para incluir/excluir del análisis sin eliminar filas.

    Si df ya contiene ROI_status y preserve_existing=True, conserva esas marcas.
    Si no, todas las filas parten con ROI_status=1. Las ROI indicadas en
    roi_exclude o scoped_exclude se marcan con ROI_status=0 solo dentro de ese
    DataFrame. Si ya filtraste por genotype antes de llamar esta función, la
    marca aplica solo a ese genotype/subconjunto.

    Parameters
    ----------
    df : pandas.DataFrame
        Tabla con una columna de ROI.
    roi_exclude : list or str, optional
        ROI a marcar como excluidas globalmente en todos los experimentos. Ejemplo:
        ["ROI7", "ROI10"].
    roi_col : str
        Nombre de la columna que identifica las ROI.
    scoped_exclude : list of dict, optional
        Marcas específicas por metadatos. Cada dict debe contener roi_col
        y cualquier otra columna que delimite el contexto. Ejemplos:

        [{"source_folder": "mut36_image8", "ROI": "ROI7"}]
        [{"sample": "sample_08", "genotype": "36", "ROI": "ROI7"}]
    status_col : str
        Nombre de la columna de estado.
    default_status : int
        Valor para ROI incluidas.
    excluded_status : int
        Valor para ROI excluidas.
    preserve_existing : bool
        Si True, conserva ROI_status existente y solo agrega nuevas marcas.

    Returns
    -------
    pandas.DataFrame
        Copia de df con la columna ROI_status.
    """
    if df.empty:
        result = df.copy()
        if status_col not in result.columns:
            result[status_col] = default_status
        return result

    if roi_col not in df.columns:
        raise ValueError(f"Falta columna de ROI: {roi_col}")

    result = df.copy()
    if status_col not in result.columns or not preserve_existing:
        result[status_col] = default_status
    else:
        result[status_col] = result[status_col].fillna(default_status).astype(int)
    exclude_mask = pd.Series(False, index=result.index)

    if roi_exclude is not None:
        if isinstance(roi_exclude, str):
            roi_exclude = [roi_exclude]
        excluded = {str(roi).strip() for roi in roi_exclude}
        exclude_mask |= result[roi_col].astype(str).str.strip().isin(excluded)

    if scoped_exclude is not None:
        if isinstance(scoped_exclude, dict):
            scoped_exclude = [scoped_exclude]

        for rule in scoped_exclude:
            if roi_col not in rule:
                raise ValueError(f"Cada regla de scoped_exclude debe incluir '{roi_col}'.")

            missing = set(rule) - set(result.columns)
            if missing:
                raise ValueError(f"Columnas no encontradas en scoped_exclude: {sorted(missing)}")

            rule_mask = pd.Series(True, index=result.index)
            for col, value in rule.items():
                rule_mask &= result[col].astype(str).str.strip().eq(str(value).strip())
            exclude_mask |= rule_mask

    result.loc[exclude_mask, status_col] = excluded_status

    marked_rows = int(exclude_mask.sum())
    if marked_rows:
        marked_rois = result.loc[exclude_mask, roi_col].dropna().astype(str).unique()
        print(f"ROIs marcadas ROI_status={excluded_status}: {sorted(marked_rois)}")
        print(f"Filas marcadas: {marked_rows}")
    else:
        print(f"ROIs marcadas ROI_status={excluded_status}: ninguna")

    return result


def attach_roi_summary_to_long_df(
    long_df,
    roi_summary,
    id_cols=None,
    roi_col="ROI",
    status_col="ROI_status",
    default_status=1,
):
    """
    Agrega métricas por ROI a un DataFrame largo.

    long_df tiene una fila por ROI/frame/TTL. roi_summary tiene una fila por
    ROI/muestra con columnas como trend, delta_high_low, low_mean, etc.
    """
    if id_cols is None:
        id_cols = [roi_col]

    missing_long = set(id_cols) - set(long_df.columns)
    missing_summary = set(id_cols) - set(roi_summary.columns)
    if missing_long:
        raise ValueError(f"Faltan columnas en long_df: {sorted(missing_long)}")
    if missing_summary:
        raise ValueError(f"Faltan columnas en roi_summary: {sorted(missing_summary)}")

    summary_cols = [
        col for col in roi_summary.columns
        if col in id_cols or col not in long_df.columns
    ]
    summary_unique = roi_summary[summary_cols].drop_duplicates(subset=id_cols)

    out = long_df.merge(
        summary_unique,
        on=id_cols,
        how="left",
    )
    if status_col not in out.columns:
        out[status_col] = default_status
    out[status_col] = out[status_col].fillna(default_status).astype(int)

    return out


def preprocess_long_for_batch(
    df_phase,
    sample=None,
    genotype=None,
    phase="cooling",
    temp_ranges=None,
    value_col="NormSignal",
    min_abs_delta=0.01,
    min_points=3,
    roi_summary=None,
    norm_temp_range=None,
):
    """
    Construye el CSV largo enriquecido usado como preprocesamiento batch.

    El resultado conserva una fila por ROI/frame/TTL e incluye:
    - phase y dT
    - trend y métricas por rango de temperatura
    - ROI_status: conserva el valor que ya traiga df_phase (por ejemplo de
      apply_manual_roi_exclusion) y solo usa 1 por defecto donde no haya nada.
    - temp_range_<nombre>_min / temp_range_<nombre>_max: metadata (mismo valor
      repetido en toda la muestra) que registra qué temp_ranges se usó
      efectivamente para calcular trend/delta_high_low, para poder detectar en
      02_processing.ipynb si dos muestras se preprocesaron con rangos distintos.
    - norm_temp_range_min / norm_temp_range_max: metadata que registra el
      temp_range (ventana de corte antes de normalizar) con el que se llamó
      process_sample() para esta muestra. Es distinto de temp_ranges
      (low/mid/high, usado para clasificar trend): temp_range acota qué
      puntos entran a Normalization() y por lo tanto qué I_baseline/Imax se
      calculan. Si norm_temp_range es None (no se usó ventana de corte), estas
      columnas quedan en NaN, para poder distinguir "no se usó ventana" de
      "CSV generado antes de este registro" (columnas ausentes).
    """
    long_df = df_phase.copy()
    if sample is not None and "sample" not in long_df.columns:
        long_df["sample"] = sample
    if genotype is not None and "genotype" not in long_df.columns:
        long_df["genotype"] = genotype

    if "phase" not in long_df.columns or "dT" not in long_df.columns:
        long_df = add_temperature_phase(long_df)

    analysis_df = long_df
    if phase is not None:
        if "phase" not in analysis_df.columns:
            raise ValueError("df_phase debe contener 'phase' para filtrar por fase.")
        analysis_df = analysis_df[analysis_df["phase"] == phase].copy()

    temp_ranges_used = temp_ranges if temp_ranges is not None else {
        "low": (22, 28),
        "mid": (32, 35),
        "high": (36, 42),
    }

    if roi_summary is None:
        trend_results = analyze_temp_trends(
            analysis_df,
            phase=None,
            temp_ranges=temp_ranges_used,
            value_col=value_col,
            min_abs_delta=min_abs_delta,
            min_points=min_points,
        )
        roi_summary = trend_results["roi_summary"]

    result = attach_roi_summary_to_long_df(
        analysis_df,
        roi_summary,
        id_cols=["ROI"],
    )

    for name, (t_min, t_max) in temp_ranges_used.items():
        result[f"temp_range_{name}_min"] = t_min
        result[f"temp_range_{name}_max"] = t_max

    result["norm_temp_range_min"] = norm_temp_range[0] if norm_temp_range is not None else np.nan
    result["norm_temp_range_max"] = norm_temp_range[1] if norm_temp_range is not None else np.nan

    return result


def export_preprocessed_long(
    df_phase,
    sample="sample",
    genotype="genotype",
    output_dir="outputs",
    phase="cooling",
    temp_ranges=None,
    value_col="NormSignal",
    min_abs_delta=0.01,
    min_points=3,
    roi_summary=None,
    filename=None,
    norm_temp_range=None,
):
    """
    Exporta un CSV largo enriquecido para análisis batch.

    norm_temp_range debe ser el mismo temp_range (t_min, t_max) con el que se
    llamó process_sample() para esta muestra, para que quede registrado en el
    CSV exportado (ver docstring de preprocess_long_for_batch).
    """
    export_df = preprocess_long_for_batch(
        df_phase,
        sample=sample,
        genotype=genotype,
        phase=phase,
        temp_ranges=temp_ranges,
        value_col=value_col,
        min_abs_delta=min_abs_delta,
        min_points=min_points,
        roi_summary=roi_summary,
        norm_temp_range=norm_temp_range,
    )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if filename is None:
        phase_label = phase if phase is not None else "all"
        filename = f"{sample}_{genotype}_{phase_label}_preprocessed_long.csv"

    output_path = output_dir / filename
    export_df.to_csv(output_path, index=False)
    print(f"Preprocesado exportado: {output_path}")

    return {
        "preprocessed": str(output_path),
        "df": export_df,
    }



def find_files(input_dir):
    """
    Busca automáticamente:
    - ROI*.csv
    - *.abf

    dentro de input_dir.

    Returns
    -------
    dict
    """

    if not os.path.isdir(input_dir):
        raise ValueError(f"No existe carpeta: {input_dir}")

    roi_files = sorted(
        glob.glob(os.path.join(input_dir, "ROI*.csv"))
    )

    abf_files = sorted(
        glob.glob(os.path.join(input_dir, "*.abf"))
    )

    print(f"Archivos encontrados en: {input_dir}")
    print("  ROI CSV:")
    if roi_files:
        for file in roi_files:
            print(f"    - {os.path.basename(file)}")
    else:
        print("    - ninguno")

    print("  ABF:")
    if abf_files:
        for file in abf_files:
            print(f"    - {os.path.basename(file)}")
    else:
        print("    - ninguno")

    if len(roi_files) == 0:
        raise FileNotFoundError("No se encontró ROI*.csv")

    if len(abf_files) == 0:
        raise FileNotFoundError("No se encontró *.abf")

    if len(roi_files) > 1:
        print(
            "Advertencia: múltiples ROI csv encontrados. "
            f"Se usará: {os.path.basename(roi_files[0])}"
        )

    if len(abf_files) > 1:
        print(
            "Advertencia: múltiples ABF encontrados. "
            f"Se usará: {os.path.basename(abf_files[0])}"
        )

    return {
        "roi_csv": roi_files[0],
        "abf_file": abf_files[0],
        "output_dir": input_dir
    }
