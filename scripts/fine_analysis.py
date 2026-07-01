"""
Fine-grained ROI analysis for temperature-synchronized fluorescence data.

This module adds a reproducible layer on top of the existing preprocessing:
- temperature-bin response curves per ROI
- curve-level metrics per ROI
- sample/genotype summaries
- lightweight QC flags

It uses only numpy/pandas so it works with the current conda environment.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_ID_COLS = [
    "source_folder",
    "source_file",
    "sample",
    "genotype",
    "genotype_meta",
    "nickname_meta",
    "ROI",
]


def load_processing_table(
    batch_dir="data/Proc_data/batch_analysis",
    filename="processing_active.csv",
):
    """Load the long active processing table produced by notebook 02."""
    path = Path(batch_dir) / filename
    if not path.exists():
        raise FileNotFoundError(f"No existe la tabla de procesamiento: {path}")
    return pd.read_csv(path, low_memory=False)


def _present_columns(df, cols):
    return [col for col in cols if col in df.columns]


def filter_roi_rows(
    df,
    phase="cooling",
    roi_status=1,
    min_temp=None,
    max_temp=None,
    value_col="NormSignal",
):
    """Apply the common filters used before fine ROI profiling."""
    out = df.copy()
    if phase is not None and "phase" in out.columns:
        out = out[out["phase"].astype(str).eq(str(phase))].copy()
    if roi_status is not None and "ROI_status" in out.columns:
        out = out[out["ROI_status"].fillna(1).astype(int).eq(int(roi_status))].copy()
    if min_temp is not None:
        out = out[out["temp_mean"] >= min_temp].copy()
    if max_temp is not None:
        out = out[out["temp_mean"] <= max_temp].copy()

    required = {"temp_mean", value_col, "ROI"}
    missing = required - set(out.columns)
    if missing:
        raise ValueError(f"Faltan columnas requeridas: {sorted(missing)}")

    out["temp_mean"] = pd.to_numeric(out["temp_mean"], errors="coerce")
    out[value_col] = pd.to_numeric(out[value_col], errors="coerce")
    return out.dropna(subset=["temp_mean", value_col, "ROI"]).copy()


def build_temperature_bins(df, bin_width=0.5, min_temp=None, max_temp=None):
    """Return deterministic temperature bin edges covering the data."""
    if df.empty:
        raise ValueError("No se pueden crear bins con una tabla vacia.")
    t_min = float(df["temp_mean"].min() if min_temp is None else min_temp)
    t_max = float(df["temp_mean"].max() if max_temp is None else max_temp)
    start = np.floor(t_min / bin_width) * bin_width
    stop = np.ceil(t_max / bin_width) * bin_width + bin_width
    return np.arange(start, stop + bin_width / 10, bin_width)


def summarize_roi_temperature_bins(
    df,
    bin_width=0.5,
    value_col="NormSignal",
    min_temp=None,
    max_temp=None,
    id_cols=None,
):
    """
    Summarize each ROI in small temperature bins.

    Returns one row per ROI x temperature bin with mean, sd, sem, n, median and
    selected quantiles. This table is useful for smooth curves and heatmaps.
    """
    if id_cols is None:
        id_cols = DEFAULT_ID_COLS
    id_cols = _present_columns(df, id_cols)
    if not id_cols:
        id_cols = ["ROI"]

    bins = build_temperature_bins(df, bin_width=bin_width, min_temp=min_temp, max_temp=max_temp)
    out = df.copy()
    out["temp_bin"] = pd.cut(out["temp_mean"], bins=bins, include_lowest=True, right=False)

    summary = (
        out.dropna(subset=["temp_bin"])
        .groupby(id_cols + ["temp_bin"], observed=True, dropna=False)
        .agg(
            temp_center=("temp_mean", "mean"),
            temp_min=("temp_mean", "min"),
            temp_max=("temp_mean", "max"),
            signal_mean=(value_col, "mean"),
            signal_sd=(value_col, "std"),
            signal_median=(value_col, "median"),
            signal_q25=(value_col, lambda x: x.quantile(0.25)),
            signal_q75=(value_col, lambda x: x.quantile(0.75)),
            n_points=(value_col, "size"),
        )
        .reset_index()
    )
    summary["signal_sem"] = summary["signal_sd"] / np.sqrt(summary["n_points"])
    summary["temp_bin_left"] = summary["temp_bin"].map(lambda interval: interval.left).astype(float)
    summary["temp_bin_right"] = summary["temp_bin"].map(lambda interval: interval.right).astype(float)
    summary["temp_bin"] = summary["temp_bin"].astype(str)
    return summary


def _first_crossing_temperature(sub, threshold, direction):
    if threshold is None or sub.empty:
        return np.nan
    ordered = sub.sort_values("temp_center")
    signal = ordered["signal_mean"]
    if direction == "decrease":
        crossed = ordered[signal <= threshold]
    else:
        crossed = ordered[signal >= threshold]
    if crossed.empty:
        return np.nan
    return float(crossed["temp_center"].iloc[0])


def _slope_by_temperature(sub):
    clean = sub.dropna(subset=["temp_center", "signal_mean"])
    if len(clean) < 2 or clean["temp_center"].nunique() < 2:
        return np.nan
    return float(np.polyfit(clean["temp_center"], clean["signal_mean"], deg=1)[0])


def _auc_by_temperature(sub):
    clean = sub.dropna(subset=["temp_center", "signal_mean"]).sort_values("temp_center")
    if len(clean) < 2:
        return np.nan
    return float(np.trapezoid(clean["signal_mean"], clean["temp_center"]))


def summarize_roi_curve_metrics(
    binned_df,
    response_threshold=0.05,
    edge_fraction=0.25,
    value_col="signal_mean",
    id_cols=None,
):
    """
    Convert binned ROI curves into one metrics row per ROI.

    Metrics include dynamic range, temperature at max/min response, AUC, slope,
    first threshold crossings, high-vs-low temperature delta, and a simple
    response class based on that high-vs-low delta.
    """
    if id_cols is None:
        id_cols = DEFAULT_ID_COLS
    id_cols = _present_columns(binned_df, id_cols)
    if not id_cols:
        id_cols = ["ROI"]

    rows = []
    for keys, sub in binned_df.groupby(id_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        sub = sub.rename(columns={value_col: "signal_mean"}).sort_values("temp_center")
        row = dict(zip(id_cols, keys))
        clean = sub.dropna(subset=["temp_center", "signal_mean"])

        if clean.empty:
            row.update(
                n_bins=0,
                temp_min=np.nan,
                temp_max=np.nan,
                signal_min=np.nan,
                signal_max=np.nan,
                signal_range=np.nan,
                temp_at_signal_min=np.nan,
                temp_at_signal_max=np.nan,
                auc_signal_temp=np.nan,
                slope_signal_per_c=np.nan,
                low_temp_signal_mean=np.nan,
                high_temp_signal_mean=np.nan,
                delta_high_low_curve=np.nan,
                threshold_temp_increase=np.nan,
                threshold_temp_decrease=np.nan,
                response_class="insufficient",
            )
            rows.append(row)
            continue

        max_idx = clean["signal_mean"].idxmax()
        min_idx = clean["signal_mean"].idxmin()
        signal_max = float(clean.loc[max_idx, "signal_mean"])
        signal_min = float(clean.loc[min_idx, "signal_mean"])
        slope = _slope_by_temperature(clean)
        n_edge_bins = max(1, int(np.ceil(len(clean) * edge_fraction)))
        low_temp_signal_mean = float(clean.head(n_edge_bins)["signal_mean"].mean())
        high_temp_signal_mean = float(clean.tail(n_edge_bins)["signal_mean"].mean())
        delta_high_low_curve = high_temp_signal_mean - low_temp_signal_mean

        if len(clean) < 3:
            response_class = "insufficient"
        elif delta_high_low_curve >= response_threshold:
            response_class = "increase"
        elif delta_high_low_curve <= -response_threshold:
            response_class = "decrease"
        else:
            response_class = "stable"

        row.update(
            n_bins=int(len(clean)),
            n_points_total=int(clean["n_points"].sum()) if "n_points" in clean.columns else np.nan,
            temp_min=float(clean["temp_center"].min()),
            temp_max=float(clean["temp_center"].max()),
            signal_min=signal_min,
            signal_max=signal_max,
            signal_range=signal_max - signal_min,
            signal_mean_over_curve=float(clean["signal_mean"].mean()),
            temp_at_signal_min=float(clean.loc[min_idx, "temp_center"]),
            temp_at_signal_max=float(clean.loc[max_idx, "temp_center"]),
            auc_signal_temp=_auc_by_temperature(clean),
            slope_signal_per_c=slope,
            low_temp_signal_mean=low_temp_signal_mean,
            high_temp_signal_mean=high_temp_signal_mean,
            delta_high_low_curve=delta_high_low_curve,
            threshold_temp_increase=_first_crossing_temperature(clean, response_threshold, "increase"),
            threshold_temp_decrease=_first_crossing_temperature(clean, -response_threshold, "decrease"),
            response_class=response_class,
        )
        rows.append(row)

    return pd.DataFrame(rows)


def add_roi_qc_flags(
    metrics_df,
    min_bins=6,
    min_points_total=15,
    max_abs_signal=5.0,
):
    """Add conservative QC flags without removing any ROI."""
    out = metrics_df.copy()
    out["qc_low_bin_count"] = out["n_bins"].fillna(0) < min_bins
    if "n_points_total" in out.columns:
        out["qc_low_point_count"] = out["n_points_total"].fillna(0) < min_points_total
    else:
        out["qc_low_point_count"] = False
    out["qc_extreme_signal"] = (
        out[["signal_min", "signal_max"]].abs().max(axis=1).fillna(0) > max_abs_signal
    )
    qc_cols = ["qc_low_bin_count", "qc_low_point_count", "qc_extreme_signal"]
    out["qc_flagged"] = out[qc_cols].any(axis=1)
    out["qc_reason"] = out[qc_cols].apply(
        lambda row: ";".join(col for col, flagged in row.items() if flagged),
        axis=1,
    )
    return out


def summarize_metrics_by_group(
    metrics_df,
    group_cols=("genotype_meta", "sample"),
):
    """Aggregate ROI-level fine metrics by sample/genotype."""
    group_cols = _present_columns(metrics_df, group_cols)
    if metrics_df.empty or not group_cols:
        return pd.DataFrame()

    valid = metrics_df[~metrics_df.get("qc_flagged", False)].copy()
    if valid.empty:
        valid = metrics_df.copy()

    summary = (
        valid.groupby(group_cols, dropna=False)
        .agg(
            n_roi=("ROI", "nunique"),
            signal_range_mean=("signal_range", "mean"),
            signal_range_sem=("signal_range", lambda x: x.std() / np.sqrt(len(x)) if len(x) else np.nan),
            slope_mean=("slope_signal_per_c", "mean"),
            slope_sem=("slope_signal_per_c", lambda x: x.std() / np.sqrt(len(x)) if len(x) else np.nan),
            auc_mean=("auc_signal_temp", "mean"),
            threshold_increase_median=("threshold_temp_increase", "median"),
            threshold_decrease_median=("threshold_temp_decrease", "median"),
            qc_flagged_n=("qc_flagged", "sum"),
        )
        .reset_index()
    )

    classes = (
        valid.groupby(group_cols + ["response_class"], dropna=False)
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    return summary.merge(classes, on=group_cols, how="left")


def run_fine_analysis(
    input_csv="data/Proc_data/batch_analysis/processing_active.csv",
    output_dir="data/Proc_data/batch_analysis/fine_analysis",
    phase="cooling",
    roi_status=1,
    bin_width=0.5,
    min_temp=None,
    max_temp=None,
    response_threshold=0.05,
):
    """Run the full fine-analysis export workflow."""
    input_csv = Path(input_csv)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_csv, low_memory=False)
    filtered = filter_roi_rows(
        df,
        phase=phase,
        roi_status=roi_status,
        min_temp=min_temp,
        max_temp=max_temp,
    )
    binned = summarize_roi_temperature_bins(
        filtered,
        bin_width=bin_width,
        min_temp=min_temp,
        max_temp=max_temp,
    )
    metrics = summarize_roi_curve_metrics(
        binned,
        response_threshold=response_threshold,
    )
    metrics = add_roi_qc_flags(metrics)
    group_summary = summarize_metrics_by_group(metrics)

    binned_path = output_dir / "roi_temperature_bins.csv"
    metrics_path = output_dir / "roi_curve_metrics.csv"
    group_path = output_dir / "roi_curve_metrics_by_group.csv"
    filtered_path = output_dir / "filtered_points_used.csv"

    filtered.to_csv(filtered_path, index=False)
    binned.to_csv(binned_path, index=False)
    metrics.to_csv(metrics_path, index=False)
    group_summary.to_csv(group_path, index=False)

    print("Analisis fino exportado:")
    print(f"  puntos filtrados: {len(filtered)} -> {filtered_path}")
    print(f"  bins ROI-temperatura: {len(binned)} -> {binned_path}")
    print(f"  metricas por ROI: {len(metrics)} -> {metrics_path}")
    print(f"  resumen por grupo: {len(group_summary)} -> {group_path}")
    return {
        "filtered": filtered,
        "binned": binned,
        "metrics": metrics,
        "group_summary": group_summary,
        "paths": {
            "filtered": filtered_path,
            "binned": binned_path,
            "metrics": metrics_path,
            "group_summary": group_path,
        },
    }


if __name__ == "__main__":
    run_fine_analysis()
