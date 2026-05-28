import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def filter_values(df, column, values=None):
    """
    Filtra por valores exactos. values=None conserva todas las filas.
    """
    if values is None or column not in df.columns:
        return df.copy()
    if isinstance(values, str):
        values = [values]
    selected = {str(value).strip() for value in values}
    return df[df[column].astype(str).str.strip().isin(selected)].copy()


def add_high_low_ratio(
    roi_summary,
    numerator_col="high_mean",
    denominator_col="low_mean",
    ratio_col="high_low_ratio",
    min_abs_denominator=1e-9,
    ratio_mode="relative_total",
):
    """
    Agrega high/low ratio y delta high-low a una tabla resumen por ROI.

    ratio_mode="relative_total" asume que low/high son DeltaF/F0 y calcula:

        (1 + high_mean) / (1 + low_mean)

    Esto compara F/F0 estimado entre rangos de temperatura.
    """
    required_cols = {numerator_col, denominator_col, "ROI"}
    missing_cols = required_cols - set(roi_summary.columns)
    if missing_cols:
        raise ValueError(f"Faltan columnas para calcular el ratio: {sorted(missing_cols)}")

    out = roi_summary.copy()
    out[numerator_col] = pd.to_numeric(out[numerator_col], errors="coerce")
    out[denominator_col] = pd.to_numeric(out[denominator_col], errors="coerce")

    if ratio_mode == "relative_total":
        out["high_relative_total"] = 1 + out[numerator_col]
        out["low_relative_total"] = 1 + out[denominator_col]
        denominator = out["low_relative_total"]
        numerator = out["high_relative_total"]
    elif ratio_mode == "raw_normsignal":
        denominator = out[denominator_col]
        numerator = out[numerator_col]
    else:
        raise ValueError("ratio_mode debe ser 'relative_total' o 'raw_normsignal'.")

    valid_denominator = denominator.abs() > min_abs_denominator
    out[ratio_col] = np.where(
        valid_denominator,
        numerator / denominator,
        np.nan,
    )
    out["raw_normsignal_high_low_ratio"] = np.where(
        out[denominator_col].abs() > min_abs_denominator,
        out[numerator_col] / out[denominator_col],
        np.nan,
    )
    out["high_minus_low"] = out[numerator_col] - out[denominator_col]
    out = out.replace([np.inf, -np.inf], np.nan)
    out["ratio_valid"] = out[ratio_col].notna()
    return out


def prepare_ratio_dataframe(
    roi_summary,
    genotype_filter=None,
    trend_filter=None,
    roi_status_filter=1,
    numerator_col="high_mean",
    denominator_col="low_mean",
    ratio_col="high_low_ratio",
    min_abs_denominator=1e-9,
    ratio_mode="relative_total",
):
    """
    Filtra el resumen por ROI y calcula el high/low ratio.
    """
    ratio_df = roi_summary.copy()
    ratio_df = filter_values(ratio_df, "genotype_meta", genotype_filter)
    ratio_df = filter_values(ratio_df, "trend", trend_filter)

    if roi_status_filter is not None and "ROI_status" in ratio_df.columns:
        ratio_df = ratio_df[
            ratio_df["ROI_status"].fillna(1).astype(int) == int(roi_status_filter)
        ].copy()

    ratio_df = add_high_low_ratio(
        ratio_df,
        numerator_col=numerator_col,
        denominator_col=denominator_col,
        ratio_col=ratio_col,
        min_abs_denominator=min_abs_denominator,
        ratio_mode=ratio_mode,
    )

    sort_cols = [col for col in ["genotype_meta", "trend", "sample", "ROI"] if col in ratio_df.columns]
    return ratio_df.sort_values(sort_cols).reset_index(drop=True)


def summarize_ratio_by_group(
    ratio_df,
    group_cols=("genotype_meta", "sample", "trend"),
    ratio_col="high_low_ratio",
    numerator_col="high_mean",
    denominator_col="low_mean",
):
    """
    Resume el ratio por grupo.
    """
    group_cols = [col for col in group_cols if col in ratio_df.columns]
    if ratio_df.empty or not group_cols:
        return pd.DataFrame()

    return (
        ratio_df[ratio_df["ratio_valid"]]
        .groupby(group_cols, dropna=False)
        .agg(
            n_roi=(ratio_col, "size"),
            ratio_mean=(ratio_col, "mean"),
            ratio_sd=(ratio_col, "std"),
            ratio_sem=(ratio_col, lambda x: x.std() / np.sqrt(len(x)) if len(x) else np.nan),
            high_mean=(numerator_col, "mean"),
            low_mean=(denominator_col, "mean"),
            delta_mean=("high_minus_low", "mean"),
        )
        .reset_index()
    )


def plot_ratio_by_trend(
    ratio_df,
    ratio_col="high_low_ratio",
    trend_col="trend",
    color_col="sample",
    trend_order=("increase", "stable", "decrease", "insufficient"),
    jitter_width=0.18,
    point_size=42,
    alpha=0.72,
    show_mean_sem=True,
    ax=None,
):
    """
    Category plot de high/low ratio por trend.
    """
    plot_df = ratio_df[ratio_df["ratio_valid"]].dropna(subset=[trend_col, ratio_col]).copy()
    if plot_df.empty:
        print("No hay ratios válidos para graficar con los filtros actuales.")
        return ax

    categories = [trend for trend in trend_order if trend in set(plot_df[trend_col].astype(str))]
    extras = sorted(set(plot_df[trend_col].astype(str)) - set(categories))
    categories.extend(extras)
    x_pos = {trend: i for i, trend in enumerate(categories)}

    if ax is None:
        fig_width = max(6, 1.35 * len(categories) + 3)
        _, ax = plt.subplots(figsize=(fig_width, 5))

    color_col = color_col if color_col in plot_df.columns else None
    if color_col is None:
        color_values = ["ROI"]
        colors = {"ROI": "tab:blue"}
    else:
        color_values = sorted(plot_df[color_col].dropna().astype(str).unique())
        cmap = plt.colormaps.get_cmap("tab20").resampled(max(len(color_values), 1))
        colors = {value: cmap(i) for i, value in enumerate(color_values)}

    rng = np.random.default_rng(7)
    for label, sub in plot_df.groupby(color_col, dropna=False) if color_col else [("ROI", plot_df)]:
        label_key = str(label)
        xs = sub[trend_col].astype(str).map(x_pos).to_numpy(dtype=float)
        jitter = rng.uniform(-jitter_width, jitter_width, size=len(sub))
        ax.scatter(
            xs + jitter,
            sub[ratio_col],
            s=point_size,
            alpha=alpha,
            color=colors.get(label_key, "0.4"),
            edgecolor="white",
            linewidth=0.4,
            label=label_key,
        )

    if show_mean_sem:
        summary = (
            plot_df.groupby(plot_df[trend_col].astype(str), dropna=False)[ratio_col]
            .agg(["mean", "std", "count"])
            .reindex(categories)
        )
        sem = summary["std"] / np.sqrt(summary["count"])
        ax.errorbar(
            range(len(categories)),
            summary["mean"],
            yerr=sem,
            color="black",
            marker="D",
            markersize=5,
            linewidth=1.8,
            capsize=4,
            label="mean +/- SEM",
        )

    ax.axhline(1, color="black", linestyle="--", linewidth=1, alpha=0.55)
    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels(categories)
    ax.set_xlabel("Trend")
    ax.set_ylabel("High / low mean NormSignal")
    ax.set_title("High/low ratio by trend")
    if color_col is not None:
        ax.legend(title=color_col, bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.figure.tight_layout()
    return ax


def plot_ratio_boxplots_by_genotype(
    ratio_df,
    ratio_col="high_low_ratio",
    trend_col="trend",
    genotype_col="genotype_meta",
    color_col="sample",
    trend_order=("increase", "stable", "decrease", "insufficient"),
    jitter_width=0.12,
    point_size=28,
    alpha=0.65,
    label_roi_points=False,
    label_col="ROI",
    legend_loc="upper left",
    legend_fontsize=5,
    legend_title_fontsize=6,
    ylim=None,
    figsize=None,
):
    """
    Genera un boxplot de high/low ratio por trend para cada genotipo presente.

    Retorna un diccionario {genotype: ax}.
    """
    plot_df = ratio_df[ratio_df["ratio_valid"]].dropna(subset=[trend_col, ratio_col]).copy()
    if plot_df.empty:
        print("No hay ratios válidos para graficar con los filtros actuales.")
        return {}

    if genotype_col not in plot_df.columns:
        plot_df[genotype_col] = "all"

    genotypes = sorted(plot_df[genotype_col].dropna().astype(str).unique())
    axes = {}

    for genotype in genotypes:
        sub_genotype = plot_df[plot_df[genotype_col].astype(str) == genotype].copy()
        categories = [trend for trend in trend_order if trend in set(sub_genotype[trend_col].astype(str))]
        extras = sorted(set(sub_genotype[trend_col].astype(str)) - set(categories))
        categories.extend(extras)
        if not categories:
            continue

        values = [
            sub_genotype.loc[sub_genotype[trend_col].astype(str) == trend, ratio_col].dropna().to_numpy()
            for trend in categories
        ]

        if figsize is None:
            fig_width = max(6, 1.35 * len(categories) + 3)
            figsize_current = (fig_width, 5)
        else:
            figsize_current = figsize
        fig, ax = plt.subplots(figsize=figsize_current)
        box = ax.boxplot(
            values,
            positions=range(len(categories)),
            widths=0.55,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 1.4},
            boxprops={"facecolor": "0.88", "edgecolor": "0.25", "linewidth": 1.1},
            whiskerprops={"color": "0.25", "linewidth": 1.1},
            capprops={"color": "0.25", "linewidth": 1.1},
        )
        for patch in box["boxes"]:
            patch.set_alpha(0.75)

        color_col_active = color_col if color_col in sub_genotype.columns else None
        if color_col_active is None:
            color_values = ["ROI"]
            colors = {"ROI": "tab:blue"}
        else:
            color_values = sorted(sub_genotype[color_col_active].dropna().astype(str).unique())
            cmap = plt.colormaps.get_cmap("tab20").resampled(max(len(color_values), 1))
            colors = {value: cmap(i) for i, value in enumerate(color_values)}

        x_pos = {trend: i for i, trend in enumerate(categories)}
        rng = np.random.default_rng(7)
        grouped = sub_genotype.groupby(color_col_active, dropna=False) if color_col_active else [("ROI", sub_genotype)]
        for label, sub in grouped:
            label_key = str(label)
            xs = sub[trend_col].astype(str).map(x_pos).to_numpy(dtype=float)
            jitter = rng.uniform(-jitter_width, jitter_width, size=len(sub))
            plot_x = xs + jitter
            ax.scatter(
                plot_x,
                sub[ratio_col],
                s=point_size,
                alpha=alpha,
                color=colors.get(label_key, "0.4"),
                edgecolor="white",
                linewidth=0.35,
                label=label_key,
                zorder=3,
            )
            if label_roi_points and label_col in sub.columns:
                for x, y, roi in zip(plot_x, sub[ratio_col], sub[label_col]):
                    ax.annotate(
                        str(roi),
                        (x, y),
                        xytext=(2, 2),
                        textcoords="offset points",
                        fontsize=5,
                        alpha=0.65,
                        zorder=4,
                    )

        counts = sub_genotype.groupby(sub_genotype[trend_col].astype(str))[ratio_col].size()
        labels = [f"{trend}\nn={int(counts.get(trend, 0))}" for trend in categories]

        ax.axhline(1, color="black", linestyle="--", linewidth=1, alpha=0.55)
        ax.set_xticks(range(len(categories)))
        ax.set_xticklabels(labels)
        ax.set_xlabel("Trend")
        ax.set_ylabel("(1 + high mean) / (1 + low mean)")
        ax.set_title(str(genotype))
        if ylim is not None:
            ax.set_ylim(ylim)
        if color_col_active is not None:
            ax.legend(
                title=color_col_active,
                loc=legend_loc,
                fontsize=legend_fontsize,
                title_fontsize=legend_title_fontsize,
                frameon=True,
                framealpha=0.82,
                borderpad=0.35,
                labelspacing=0.25,
                handlelength=1,
                handletextpad=0.35,
            )
        fig.tight_layout()
        axes[genotype] = ax

    return axes


def plot_low_vs_high_by_sample_trend(
    ratio_df,
    low_col="low_mean",
    high_col="high_mean",
    sample_col="sample",
    trend_col="trend",
    genotype_col="genotype_meta",
    point_size=52,
    alpha=0.8,
    figsize=(6.5, 5),
):
    """
    Scatter low_mean vs high_mean con color por sample y forma por trend.

    Genera un gráfico separado por genotipo y retorna {genotype: ax}.
    """
    required = {low_col, high_col}
    missing = required - set(ratio_df.columns)
    if missing:
        raise ValueError(f"Faltan columnas para graficar: {sorted(missing)}")

    plot_df = ratio_df.dropna(subset=[low_col, high_col]).copy()
    if plot_df.empty:
        print("No hay datos válidos para graficar low vs high.")
        return {}

    if genotype_col not in plot_df.columns:
        plot_df[genotype_col] = "all"

    sample_col_active = sample_col if sample_col in plot_df.columns else None
    trend_col_active = trend_col if trend_col in plot_df.columns else None

    if sample_col_active is None:
        samples = ["sample"]
        colors = {"sample": "tab:blue"}
        plot_df["_sample_plot"] = "sample"
        sample_col_active = "_sample_plot"
    else:
        samples = sorted(plot_df[sample_col_active].dropna().astype(str).unique())
        cmap = plt.colormaps.get_cmap("tab20").resampled(max(len(samples), 1))
        colors = {sample: cmap(i) for i, sample in enumerate(samples)}

    marker_cycle = ["o", "s", "^", "D", "v", "P", "X", "*"]
    if trend_col_active is None:
        trends = ["trend"]
        markers = {"trend": "o"}
        plot_df["_trend_plot"] = "trend"
        trend_col_active = "_trend_plot"
    else:
        preferred_order = ["increase", "stable", "decrease", "insufficient"]
        present = set(plot_df[trend_col_active].dropna().astype(str))
        trends = [trend for trend in preferred_order if trend in present]
        trends.extend(sorted(present - set(trends)))
        markers = {
            trend: marker_cycle[i % len(marker_cycle)]
            for i, trend in enumerate(trends)
        }

    genotypes = sorted(plot_df[genotype_col].dropna().astype(str).unique())
    axes = {}

    for genotype in genotypes:
        sub_genotype = plot_df[plot_df[genotype_col].astype(str).eq(str(genotype))].copy()
        _, ax = plt.subplots(figsize=figsize)

        for sample in samples:
            for trend in trends:
                sub = sub_genotype[
                    sub_genotype[sample_col_active].astype(str).eq(str(sample))
                    & sub_genotype[trend_col_active].astype(str).eq(str(trend))
                ]
                if sub.empty:
                    continue
                ax.scatter(
                    sub[low_col],
                    sub[high_col],
                    s=point_size,
                    alpha=alpha,
                    color=colors.get(str(sample), "0.4"),
                    marker=markers.get(str(trend), "o"),
                    edgecolor="white",
                    linewidth=0.45,
                )

        lim_min = np.nanmin([sub_genotype[low_col].min(), sub_genotype[high_col].min()])
        lim_max = np.nanmax([sub_genotype[low_col].max(), sub_genotype[high_col].max()])
        pad = (lim_max - lim_min) * 0.08 if lim_max > lim_min else 0.05
        lims = [lim_min - pad, lim_max + pad]
        ax.plot(lims, lims, color="black", linestyle="--", linewidth=1, alpha=0.5)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_xlabel("Low mean NormSignal")
        ax.set_ylabel("High mean NormSignal")
        ax.set_title(str(genotype))

        sample_handles = [
            plt.Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                markerfacecolor=colors[sample],
                markeredgecolor="white",
                markersize=6,
                label=sample,
            )
            for sample in samples
            if sample in set(sub_genotype[sample_col_active].astype(str))
        ]
        trend_handles = [
            plt.Line2D(
                [0],
                [0],
                marker=markers[trend],
                linestyle="none",
                color="0.25",
                markerfacecolor="0.65",
                markeredgecolor="white",
                markersize=6,
                label=trend,
            )
            for trend in trends
            if trend in set(sub_genotype[trend_col_active].astype(str))
        ]

        legend_samples = ax.legend(
            handles=sample_handles,
            title=sample_col,
            loc="upper left",
            fontsize=6,
            title_fontsize=7,
            framealpha=0.85,
        )
        ax.add_artist(legend_samples)
        ax.legend(
            handles=trend_handles,
            title=trend_col,
            loc="lower right",
            fontsize=6,
            title_fontsize=7,
            framealpha=0.85,
        )

        ax.figure.tight_layout()
        axes[genotype] = ax

    return axes
