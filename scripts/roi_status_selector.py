from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.path import Path as MplPath
from matplotlib.widgets import LassoSelector


def apply_roi_status(
    df,
    status_df,
    id_cols,
    status_col="ROI_status",
):
    """
    Propaga ROI_status desde una tabla resumen hacia otra tabla por id_cols.
    """
    missing_df = set(id_cols) - set(df.columns)
    missing_status = set(id_cols + [status_col]) - set(status_df.columns)
    if missing_df:
        raise ValueError(f"Faltan columnas en df: {sorted(missing_df)}")
    if missing_status:
        raise ValueError(f"Faltan columnas en status_df: {sorted(missing_status)}")

    status_map = status_df[id_cols + [status_col]].drop_duplicates(subset=id_cols)
    out = df.drop(columns=[status_col], errors="ignore").merge(
        status_map,
        on=id_cols,
        how="left",
    )
    out[status_col] = out[status_col].fillna(1).astype(int)
    return out


def update_roi_status_registry(
    registry,
    updates,
    id_cols,
    status_col="ROI_status",
):
    """
    Actualiza/acumula un registro de ROI_status sin perder marcas previas.

    Las filas de updates reemplazan solo las mismas claves id_cols en registry.
    Las marcas de otros genotipos, muestras, tendencias o filtros se conservan.
    """
    required = set(id_cols + [status_col])
    missing_updates = required - set(updates.columns)
    if missing_updates:
        raise ValueError(f"Faltan columnas en updates: {sorted(missing_updates)}")

    updates_clean = updates[id_cols + [status_col]].drop_duplicates(subset=id_cols).copy()
    updates_clean[status_col] = updates_clean[status_col].fillna(1).astype(int)

    if registry is None or registry.empty:
        return updates_clean.reset_index(drop=True)

    missing_registry = required - set(registry.columns)
    if missing_registry:
        raise ValueError(f"Faltan columnas en registry: {sorted(missing_registry)}")

    registry_clean = registry[id_cols + [status_col]].drop_duplicates(subset=id_cols).copy()
    key_index = pd.MultiIndex.from_frame(updates_clean[id_cols].astype(str))
    registry_index = pd.MultiIndex.from_frame(registry_clean[id_cols].astype(str))
    registry_keep = registry_clean.loc[~registry_index.isin(key_index)].copy()

    return pd.concat([registry_keep, updates_clean], ignore_index=True)


def apply_roi_status_rules(
    df,
    exclude_rules=None,
    include_rules=None,
    status_col="ROI_status",
    default_status=1,
):
    """
    Aplica reglas reproducibles para incluir/excluir ROIs.

    Cada regla es un dict con columnas y valores exactos. Ejemplos:

        {"sample": "sample_19", "ROI": "ROI3"}
        {"genotype": "m65", "trend": "stable"}

    exclude_rules marca ROI_status=0; include_rules marca ROI_status=1.
    Las reglas se aplican en ese orden, por lo que include_rules puede rescatar
    casos excluidos por una regla más general.
    """
    out = df.copy()
    if status_col not in out.columns:
        out[status_col] = default_status
    out[status_col] = out[status_col].fillna(default_status).astype(int)

    def apply_rules(rules, value):
        nonlocal out
        if rules is None:
            return
        if isinstance(rules, dict):
            rules = [rules]

        for rule in rules:
            missing = set(rule) - set(out.columns)
            if missing:
                raise ValueError(f"Columnas no encontradas en regla ROI_status: {sorted(missing)}")

            mask = pd.Series(True, index=out.index)
            for col, expected in rule.items():
                if isinstance(expected, (list, tuple, set)):
                    selected = {str(item).strip() for item in expected}
                    mask &= out[col].astype(str).str.strip().isin(selected)
                else:
                    mask &= out[col].astype(str).str.strip().eq(str(expected).strip())
            out.loc[mask, status_col] = value

    apply_rules(exclude_rules, 0)
    apply_rules(include_rules, 1)
    return out


def roi_status_registry_from_summary(
    roi_summary,
    id_cols,
    status_col="ROI_status",
):
    """
    Construye un registry acumulable desde un resumen por ROI.
    """
    required = set(id_cols + [status_col])
    missing = required - set(roi_summary.columns)
    if missing:
        raise ValueError(f"Faltan columnas para construir registry: {sorted(missing)}")

    registry = roi_summary[id_cols + [status_col]].drop_duplicates(subset=id_cols).copy()
    registry[status_col] = registry[status_col].fillna(1).astype(int)
    return registry.reset_index(drop=True)


def load_roi_status_registry(
    batch_dir,
    id_cols,
    status_col="ROI_status",
    registry_filename="roi_status_registry.csv",
    seed_filename="roi_temp_summary_with_roi_status.csv",
):
    """
    Carga el registro acumulado de ROI_status.

    Si no existe, intenta inicializarlo desde una selección previa guardada por
    el selector. Esto permite migrar corridas antiguas sin perder marcas.
    """
    batch_dir = Path(batch_dir)
    batch_dir.mkdir(parents=True, exist_ok=True)
    registry_path = batch_dir / registry_filename

    required_cols = list(id_cols) + [status_col]
    if registry_path.exists():
        registry = pd.read_csv(registry_path)
        print(f"ROI_status registry cargado: {registry_path}")
    else:
        seed_path = batch_dir / seed_filename
        if seed_path.exists():
            seed = pd.read_csv(seed_path)
            if set(required_cols).issubset(seed.columns):
                registry = (
                    seed[required_cols]
                    .drop_duplicates(subset=id_cols)
                    .copy()
                )
                registry[status_col] = registry[status_col].fillna(1).astype(int)
                registry.to_csv(registry_path, index=False)
                print(f"ROI_status registry creado desde selección previa: {seed_path}")
            else:
                registry = pd.DataFrame(columns=required_cols)
                print("ROI_status registry iniciado vacío; la selección previa no tiene las columnas esperadas.")
        else:
            registry = pd.DataFrame(columns=required_cols)
            print("ROI_status registry iniciado vacío.")

    if not registry.empty:
        missing = set(required_cols) - set(registry.columns)
        if missing:
            raise ValueError(f"Faltan columnas en registry: {sorted(missing)}")
        registry = registry[required_cols].drop_duplicates(subset=id_cols).copy()
        registry[status_col] = registry[status_col].fillna(1).astype(int)

    return registry, registry_path


def summarize_roi_temperature_table(
    df,
    id_cols,
    group_cols=None,
    status_col="ROI_status",
):
    """
    Reduce una tabla larga a una fila por ROI/muestra con las métricas de temperatura.
    """
    if group_cols is None:
        group_cols = [
            "source_folder",
            "source_file",
            "sample",
            "genotype",
            "genotype_meta",
            "nickname_meta",
        ]

    temp_range_cols = [
        col for col in df.columns
        if col.startswith("temp_range_") or col.startswith("norm_temp_range_")
    ]

    summary_cols = list(group_cols) + [
        "ROI",
        "low_mean", "low_sd", "low_n",
        "mid_mean", "mid_sd", "mid_n",
        "high_mean", "high_sd", "high_n",
        "delta_high_low", "trend", status_col,
    ] + temp_range_cols
    summary_cols = [col for col in summary_cols if col in df.columns]

    missing_ids = set(id_cols) - set(df.columns)
    if missing_ids:
        raise ValueError(f"Faltan columnas id en df: {sorted(missing_ids)}")

    summary = (
        df[summary_cols]
        .drop_duplicates(subset=id_cols)
        .reset_index(drop=True)
    )
    if status_col not in summary.columns:
        summary[status_col] = 1
    summary[status_col] = summary[status_col].fillna(1).astype(int)
    return summary


def prepare_roi_processing_selection(
    preprocessed_all,
    registry,
    id_cols,
    filter_func,
    temp_summary_to_long_func,
    genotype_filter=None,
    phase_filter="cooling",
    trend_filter=None,
    temp_ranges=None,
    status_col="ROI_status",
):
    """
    Aplica filtros de procesamiento y propaga ROI_status acumulado al subset.
    """
    if temp_ranges is None:
        temp_ranges = {
            "low": (22, 28),
            "mid": (32, 35),
            "high": (36, 42),
        }

    selected = preprocessed_all.copy()
    selected = filter_func(selected, "genotype_meta", genotype_filter)
    selected = filter_func(selected, "phase", phase_filter)
    selected = filter_func(selected, "trend", trend_filter)

    if registry is not None and not registry.empty:
        selected = apply_roi_status(
            selected,
            registry,
            id_cols=id_cols,
            status_col=status_col,
        )
    elif status_col not in selected.columns:
        selected[status_col] = 1

    selected[status_col] = selected[status_col].fillna(1).astype(int)
    active = selected[selected[status_col] == 1].copy()

    roi_summary = summarize_roi_temperature_table(
        selected,
        id_cols=id_cols,
        status_col=status_col,
    )
    roi_summary_active = roi_summary[roi_summary[status_col] == 1].copy()

    range_long = temp_summary_to_long_func(
        roi_summary_active,
        id_cols=list(id_cols) + ["delta_high_low", "trend"],
        value_prefixes=list(temp_ranges.keys()),
    )

    return {
        "processing_selected": selected,
        "processing_active": active,
        "roi_temp_summary": roi_summary,
        "roi_temp_summary_active": roi_summary_active,
        "range_long": range_long,
        "temp_ranges": temp_ranges,
    }


def select_and_update_roi_status(
    roi_summary,
    registry,
    id_cols,
    registry_path,
    x_col="low_mean",
    y_col="high_mean",
    label_col="ROI",
    title="ROI selector: high vs low response",
    annotate=True,
    color_col="sample",
    save_path=None,
    status_col="ROI_status",
):
    """
    Ejecuta el selector interactivo y persiste el registry acumulado.
    """
    selected_summary = select_roi_status(
        roi_summary,
        x_col=x_col,
        y_col=y_col,
        label_col=label_col,
        id_cols=id_cols,
        title=title,
        annotate=annotate,
        color_col=color_col,
        save_path=save_path,
        status_col=status_col,
    )

    updated_registry = update_roi_status_registry(
        registry,
        selected_summary,
        id_cols=id_cols,
        status_col=status_col,
    )
    Path(registry_path).parent.mkdir(parents=True, exist_ok=True)
    updated_registry.to_csv(registry_path, index=False)
    print(f"ROI_status registry actualizado y guardado: {registry_path}")

    return selected_summary, updated_registry


def apply_registry_to_processing_tables(
    tables,
    registry,
    id_cols,
    temp_summary_to_long_func,
    temp_ranges,
    status_col="ROI_status",
):
    """
    Reaplica el registry a las tablas de procesamiento después de usar el selector.
    """
    selected = tables["processing_selected"]
    if registry is not None and not registry.empty:
        selected = apply_roi_status(
            selected,
            registry,
            id_cols=id_cols,
            status_col=status_col,
        )
    selected[status_col] = selected[status_col].fillna(1).astype(int)

    active = selected[selected[status_col] == 1].copy()
    roi_summary = summarize_roi_temperature_table(
        selected,
        id_cols=id_cols,
        status_col=status_col,
    )
    roi_summary_active = roi_summary[roi_summary[status_col] == 1].copy()
    range_long = temp_summary_to_long_func(
        roi_summary_active,
        id_cols=list(id_cols) + ["delta_high_low", "trend"],
        value_prefixes=list(temp_ranges.keys()),
    )

    out = dict(tables)
    out.update({
        "processing_selected": selected,
        "processing_active": active,
        "roi_temp_summary": roi_summary,
        "roi_temp_summary_active": roi_summary_active,
        "range_long": range_long,
    })
    return out


def build_master_roi_outputs(
    preprocessed_all,
    registry,
    id_cols,
    temp_summary_to_long_func,
    temp_ranges,
    status_col="ROI_status",
):
    """
    Construye las salidas globales acumuladas usadas como archivos maestros.
    """
    all_with_status = preprocessed_all.copy()
    if registry is not None and not registry.empty:
        all_with_status = apply_roi_status(
            all_with_status,
            registry,
            id_cols=id_cols,
            status_col=status_col,
        )
    elif status_col not in all_with_status.columns:
        all_with_status[status_col] = 1

    all_with_status[status_col] = all_with_status[status_col].fillna(1).astype(int)
    all_active = all_with_status[all_with_status[status_col] == 1].copy()
    roi_summary_all = summarize_roi_temperature_table(
        all_with_status,
        id_cols=id_cols,
        status_col=status_col,
    )
    roi_summary_all_active = roi_summary_all[roi_summary_all[status_col] == 1].copy()
    range_long_all_active = temp_summary_to_long_func(
        roi_summary_all_active,
        id_cols=list(id_cols) + ["delta_high_low", "trend"],
        value_prefixes=list(temp_ranges.keys()),
    )

    return {
        "all_with_roi_status": all_with_status,
        "all_active": all_active,
        "roi_temp_summary_all": roi_summary_all,
        "roi_temp_summary_all_active": roi_summary_all_active,
        "range_long_all_active": range_long_all_active,
    }


def save_roi_processing_outputs(
    current_tables,
    master_tables,
    save_func,
    base_dir,
):
    """
    Guarda la salida principal y las tablas auxiliares del procesamiento ROI.
    """
    saved_paths = {}

    output_map = {
        "roi_temp_summary_active.csv": master_tables["roi_temp_summary_all_active"],
        "processing_current_filter_with_roi_status.csv": current_tables["processing_selected"],
        "processing_current_filter_active.csv": current_tables["processing_active"],
        "roi_temp_summary_current_filter_with_roi_status.csv": current_tables["roi_temp_summary"],
        "roi_temp_summary_current_filter_active.csv": current_tables["roi_temp_summary_active"],
        "roi_temp_summary_current_filter_long.csv": current_tables["range_long"],
        "processing_all_with_roi_status.csv": master_tables["all_with_roi_status"],
        "processing_all_active.csv": master_tables["all_active"],
        "roi_temp_summary_all_with_roi_status.csv": master_tables["roi_temp_summary_all"],
        "roi_temp_summary_all_active.csv": master_tables["roi_temp_summary_all_active"],
        "roi_temp_summary_long.csv": master_tables["range_long_all_active"],
    }

    for filename, df in output_map.items():
        saved_paths[filename] = save_func(df, filename, base_dir=base_dir)

    print("Salida principal: roi_temp_summary_active.csv")
    print("ROI_status acumulado en todos los datos:")
    print(master_tables["roi_temp_summary_all"]["ROI_status"].value_counts().sort_index())
    print("ROIs activas en salida principal:", master_tables["roi_temp_summary_all_active"].shape[0])

    return saved_paths


def select_roi_status(
    df,
    x_col="low_mean",
    y_col="high_mean",
    label_col="ROI",
    status_col="ROI_status",
    id_cols=None,
    title="ROI selector",
    annotate=False,
    color_col=None,
    save_path=None,
):
    """
    Selector interactivo para marcar ROI_status en un scatter plot.

    Uso:
    - Dibuja un lazo/polígono para seleccionar puntos.
    - Tecla 0: los siguientes lazos marcan ROI_status=0.
    - Tecla 1: los siguientes lazos marcan ROI_status=1.
    - Tecla u: deshacer el último cambio.
    - Tecla s: guardar CSV si save_path fue entregado.
    - Tecla q: cerrar la ventana y devolver el DataFrame actualizado.

    Si color_col se entrega, los puntos activos se colorean por esa columna.
    Los puntos con ROI_status=0 se muestran en gris.

    La función bloquea la celda hasta que se cierra la ventana.
    """
    required = {x_col, y_col, label_col}
    if color_col is not None:
        required.add(color_col)
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Faltan columnas para graficar: {sorted(missing)}")

    out = df.copy()
    if status_col not in out.columns:
        out[status_col] = 1
    out[status_col] = out[status_col].fillna(1).astype(int)

    plot_df = out.dropna(subset=[x_col, y_col]).copy()
    if plot_df.empty:
        print("No hay puntos para seleccionar.")
        return out

    current_status = {"value": 0}
    history = []

    fig, ax = plt.subplots(figsize=(7, 6))

    if color_col is None:
        color_map = {}
    else:
        color_values = sorted(plot_df[color_col].dropna().astype(str).unique())
        cmap = plt.get_cmap("tab20", max(len(color_values), 1))
        color_map = {
            value: cmap(i)
            for i, value in enumerate(color_values)
        }

    def colors():
        if color_col is None:
            return [
                "tab:blue" if status == 1 else "0.75"
                for status in plot_df[status_col].to_numpy()
            ]

        return [
            color_map.get(str(value), "tab:blue") if status == 1 else "0.75"
            for value, status in zip(
                plot_df[color_col].astype(str),
                plot_df[status_col].to_numpy(),
            )
        ]

    scatter = ax.scatter(
        plot_df[x_col],
        plot_df[y_col],
        c=colors(),
        s=35,
        alpha=0.8,
        edgecolors="black",
        linewidths=0.3,
    )

    if annotate:
        for _, row in plot_df.iterrows():
            ax.annotate(
                str(row[label_col]),
                (row[x_col], row[y_col]),
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=7,
                alpha=0.75,
            )

    ax.axline((0, 0), slope=1, color="black", linestyle="--", alpha=0.35)
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title(title)

    if color_col is not None and color_map:
        legend_handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor=color,
                markeredgecolor="black",
                markersize=7,
                label=label,
            )
            for label, color in color_map.items()
        ]
        legend_handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor="0.75",
                markeredgecolor="black",
                markersize=7,
                label=f"{status_col}=0",
            )
        )
        ax.legend(
            handles=legend_handles,
            title=color_col,
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
        )

    status_text = ax.text(
        0.01,
        0.99,
        "",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "0.8"},
    )

    def sync_main_df():
        out.loc[plot_df.index, status_col] = plot_df[status_col].to_numpy()

    def save():
        if save_path is None:
            print("save_path no definido; no se guardó archivo.")
            return
        sync_main_df()
        path = Path(save_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(path, index=False)
        print(f"Guardado: {path}")

    def refresh():
        scatter.set_color(colors())
        n_active = int((plot_df[status_col] == 1).sum())
        n_excluded = int((plot_df[status_col] == 0).sum())
        status_text.set_text(
            "Lazo: marca ROI_status={status}\n"
            "0 excluir | 1 incluir | u deshacer | s guardar | q cerrar\n"
            "incluidas: {active} | excluidas: {excluded}".format(
                status=current_status["value"],
                active=n_active,
                excluded=n_excluded,
            )
        )
        fig.canvas.draw_idle()

    def on_select(vertices):
        path = MplPath(vertices)
        points = plot_df[[x_col, y_col]].to_numpy(dtype=float)
        selected = path.contains_points(points)
        if not selected.any():
            return

        history.append(plot_df[status_col].copy())
        plot_df.loc[selected, status_col] = current_status["value"]
        sync_main_df()

        selected_labels = sorted(plot_df.loc[selected, label_col].astype(str).unique())
        print(f"ROI_status={current_status['value']}:", selected_labels)
        refresh()

    def on_key(event):
        if event.key == "0":
            current_status["value"] = 0
            refresh()
        elif event.key == "1":
            current_status["value"] = 1
            refresh()
        elif event.key == "u" and history:
            plot_df[status_col] = history.pop()
            sync_main_df()
            refresh()
        elif event.key == "s":
            save()
        elif event.key == "q":
            plt.close(fig)

    lasso = LassoSelector(ax, on_select)
    fig.canvas.mpl_connect("key_press_event", on_key)
    refresh()
    plt.show(block=True)

    # Keep a reference alive until the window closes.
    _ = lasso
    sync_main_df()
    return out
