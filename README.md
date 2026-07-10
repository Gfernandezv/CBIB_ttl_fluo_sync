# TRPM3 Fluorescence TTL Temperature Analysis

Analysis workflow for TRPM3 fluorescence imaging synchronized with TTL and temperature recordings. The project takes raw ABF (temperature/TTL) + ImageJ ROI fluorescence recordings, aligns them frame by frame, normalizes and curates the signal per ROI, and produces genotype/sample-level plots of fluorescence response across temperature — plus an exploratory tool to visualize how ROI intensity looks across the temperature range directly on the images.

## Project Structure

```text
.
├── environment.yml
├── notebooks/
│   ├── 01_preprocessing.ipynb   # ABF+ROI alignment, normalization, curation, export
│   ├── 02_processing.ipynb      # compile all samples, batch reprocessing utilities
│   ├── 03_analysis.ipynb        # genotype-level NormSignal-vs-temperature curves
│   └── 04_temperature_tiff.ipynb  # exploratory: per-ROI images ordered by temperature
├── scripts/
│   ├── analisis_ttl.py          # core pipeline: ABF/TTL, normalization, curation, export
│   └── roi_image_crops.py       # ROI image crops, temperature-sorted stacks, montages
├── data/              # local data, not tracked by git
├── graphs/            # optional exported figures, not tracked by git
└── processed/         # local intermediate outputs, not tracked by git
```


## Environment

Create the conda environment:

```bash
conda env create -f environment.yml
conda activate ttl_imagej
```

Then open Jupyter from the project folder:

```bash
jupyter lab
```

## Expected Data Layout

Raw and processed data are expected locally under:

```text
data/Proc_data/
```

Each sample folder should contain the ABF file, the ImageJ ROI results CSV, `RoiSet.zip`, and the `.tif`/`.tiff` image stack:

```text
data/Proc_data/
├── mut27_image10/
│   ├── <recording>.abf
│   ├── ROI_Results.csv
│   ├── RoiSet.zip
│   └── <image_stack>.tif
├── mut27_image11/
├── mut36_image8/
└── ...
```

Large data files are intentionally ignored by git.

## Workflow

The pipeline is a straight chain: each notebook's output feeds the next one directly through `*_preprocessed_long.csv`, no external Excel/registry file is required.

### 1. Preprocessing — `notebooks/01_preprocessing.ipynb`

Per sample folder, this notebook:

- loads the ABF temperature/TTL recording (`ttl.load_abf_temp_ttl`) and detects TTL pulses (`ttl.detect_ttl_events`);
- aligns ImageJ ROI fluorescence frames with TTL/temperature, frame by frame (`ttl.process_sample`, wrapping `ttl.long_format` + `ttl.Normalization` + `ttl.add_temperature_phase`);
- normalizes each ROI's signal to a baseline temperature (`NormSignal = (IntDen - I_baseline) / I_baseline`);
- assigns `phase` (heating/cooling) and `trend` (increase/decrease/stable) labels per ROI;
- curates ROIs: marks `ROI_status` (1 = included) via `ttl.mark_roi_status`, optionally helped by visual crops of each candidate ROI (`crops.crop_included_rois_from_notebook` + `crops.list_rois_from_crops_folder`);
- exports one `*_preprocessed_long.csv` per sample (`ttl.export_preprocessed_long`) — one row per ROI/frame/TTL, with `phase`, `trend`, low/mid/high range metrics, and the final `ROI_status` already baked in.

This is the only notebook that touches raw ABF/TIFF/ROI files; everything downstream reads the exported `_preprocessed_long.csv`.

### 2. Compilation and batch reprocessing — `notebooks/02_processing.ipynb`

- compiles every `*_preprocessed_long.csv` under `data/Proc_data/` into a single DataFrame (`ttl.load_all_preprocessed_long`), flagging folders with more than one exported CSV (leftover from a previous run of notebook 01);
- reads back the `ROI_status == 1` ROIs of a given sample straight from the compiled table (`ttl.selected_rois_from_long`), e.g. to feed `selected_rois` into a reprocessing call;
- reprocesses one or several sample folders at once from the raw ABF/TIFF, by folder-name pattern (`ttl.process_selected_folders`), for cases where parameters need to change after the fact.

### 3. Analysis — `notebooks/03_analysis.ipynb`

- starts from `preprocessed_all` (the compiled table from step 2) filtered to `ROI_status == 1` — no recomputation of normalization/phase/trend;
- filters by `phase`, `trend`, genotype, and temperature range (parameters at the top of the notebook);
- bins `NormSignal` by temperature per ROI, then summarizes mean ± SEM per genotype;
- plots per-ROI fine curves plus the genotype-level mean curve across temperature.

There is intentionally no Arrhenius-style (`ln` vs `1/T`) phase here: `NormSignal` is a steady-state fluorescence amplitude, not a reaction rate constant, so that analysis wouldn't have a valid physical interpretation.

### 4. Temperature-sorted ROI viewer (exploratory) — `notebooks/04_temperature_tiff.ipynb`

For one chosen sample folder:

- plots `NormSignal` vs temperature for its `ROI_status == 1` ROIs (`ttl.graph_selected_rois_by_folder`);
- builds, per ROI, a multi-frame `.tif` with the **same crop as the original acquisition order, but frames reordered from lowest to highest temperature** (heating/cooling mixed, sorted only by `temp_mean`) — `crops.crop_selected_rois_sorted_by_temperature` — so it plays back like a "temperature movie" in ImageJ/Fiji;
- since Python has no built-in interactive TIFF viewer, also builds a **montage collage**: one representative frame every `temp_bin_width` °C, laid out as a single wide strip (`crops.select_indices_by_temp_bin` + `crops.save_montage_png`), so the whole temperature range can be inspected in one PNG image without opening ImageJ.

Outputs are written to `data/Proc_data/<folder>/ROIs_by_temp/`.

## Main Outputs

Per sample, under `data/Proc_data/<folder>/`:

```text
<sample>_<genotype>_<phase>_preprocessed_long.csv   # canonical output of notebook 01
ROIs/                                                # optional: original-order ROI crops (manual curation aid)
ROIs_by_temp/                                        # optional: temperature-sorted crops + montages (notebook 04)
```

## Notes

- `data/`, `raw/`, `processed/`, and `graphs/` are ignored by git to avoid uploading large local files.
- The interactive ROI curation cells in `01_preprocessing.ipynb` use a local Matplotlib Qt window (`%matplotlib qt`).
- Analysis logic lives in `scripts/`, so notebooks stay focused on workflow and parameter choices. See the function reference below for what each module provides.

## Function Reference

Full parameter documentation lives in each function's docstring; this is a map of what's available and where, grouped by purpose.

### `scripts/analisis_ttl.py` — core pipeline

**ABF / TTL / temperature**

| Function | Purpose |
|---|---|
| `load_abf_temp_ttl(file_path, temp_channel, ttl_channel)` | Load temperature and TTL channels from an ABF file. |
| `detect_ttl_events(ttl, time, threshold, min_duration_s)` | Detect TTL pulse onsets/offsets. |
| `mean_temperature_during_ttl(temperature, time, ttl_events)` | Average temperature within each detected TTL pulse. |
| `plot_temp_ttl(data, temp_ttl_df, frame_counts, only_used_frames)` | Quick-look plot of temperature/TTL alongside frame usage. |

**Alignment, normalization, phase/trend**

| Function | Purpose |
|---|---|
| `long_format(temp_ttl_df, fluo_file, start_idx, drop_frames_without_ttl)` | Merge the ImageJ Multi Measure CSV with TTL/temperature into one long DataFrame (one row per ROI/frame). |
| `Normalization(long_df, target_temp, tol, norm_type)` | Compute `NormSignal`, `I_baseline`, `Imax`, `T_at_Imax` per ROI from a baseline temperature window. |
| `add_temperature_phase(df)` | Derive `phase` (heating/cooling) from the sign of the temperature slope. |
| `filter_analysis_window(df_phase, ttl_range, time_range)` | Restrict rows to a TTL/time window without breaking frame-TTL alignment. |

**ROI curation and inclusion status**

| Function | Purpose |
|---|---|
| `mark_roi_status(df, roi_exclude, roi_col, scoped_exclude, status_col, default_status, excluded_status, preserve_existing)` | Flag ROIs as included/excluded via a `status_col` (default `ROI_status`) without dropping rows. |
| `filter_rois_by_trend(roi_summary, trend, min_abs_delta, min_points, verbose)` | Keep only ROIs whose classified trend matches a target (increase/decrease/stable). |
| `detect_roi_signal_outliers(df, value_col, group_col, order_col, window, z_thresh, min_abs_residual)` | Flag anomalous points per ROI using a robust rolling-median residual (MAD-scaled z-score). |
| `selected_rois_from_long(df, folder, folder_col, status_col, included_status)` | Read back the `ROI_status == 1` ROI names for a folder directly from a compiled `_long` table. |

**Summaries and trend classification**

| Function | Purpose |
|---|---|
| `summarize_rois_by_temp_ranges(df_phase, phase, temp_ranges, value_col, min_points, min_abs_delta, verbose)` | Per-ROI low/mid/high range means/SDs and deltas. |
| `analyze_temp_trends(df_phase, phase, temp_ranges, value_col, min_abs_delta, min_points, verbose)` | Wraps the above and classifies each ROI's `trend`. |
| `attach_roi_summary_to_long_df(long_df, roi_summary, id_cols, roi_col, status_col, default_status)` | Join per-ROI summary metrics back onto the long (frame-level) DataFrame. |

**Sample-level processing and export**

| Function | Purpose |
|---|---|
| `process_sample(input_dir, start_ttl, threshold, target_temp, norm_type, temp_channel, ttl_channel, drop_frames_without_ttl, phase_filter, temp_range, selected_rois)` | End-to-end: ABF+ROI CSV in a folder → aligned, normalized, phase/trend-labeled long DataFrame. The main entry point used interactively in `01_preprocessing.ipynb`. |
| `process_selected_folders(base_dir, folders, preprocessed_all, **process_sample_kwargs)` | Run `process_sample` over several folders (by name or fnmatch pattern) and concatenate results. |
| `preprocess_long_for_batch(...)` / `export_preprocessed_long(...)` | Build/export the enriched long CSV (`*_preprocessed_long.csv`) that downstream notebooks read. |
| `export_results(df_norm, selected, sample, genotype, output_dir, export_filtered, roi_summary)` | Export full normalized data and, optionally, a filtered subset. |
| `load_all_preprocessed_long(base_dir, pattern)` | Concatenate every `*_preprocessed_long.csv` under `base_dir` into one DataFrame (`preprocessed_all`), with per-folder duplicate-file diagnostics (`load_status`). |
| `find_files(input_dir)` | Locate the ABF, ROI CSV, and `RoiSet.zip` expected inside a sample folder. |

**Plotting**

| Function | Purpose |
|---|---|
| `graph_rois_before_filter(df_phase, phase, value_col, alpha)` | Plot every ROI (pre-curation) to eyeball candidates before filtering. |
| `graph_selected_rois_by_folder(df, folder, rois, folder_col, status_col, included_status, value_col, x_col, alpha, show_mean_se, n_bins)` | Plot `value_col` vs `x_col` for one or more samples/ROIs, with an optional binned mean ± SEM overlay. Used in notebooks 03 and 04. |
| `graph_selected_rois_by_phase(df_phase, selected, value_col, phases, phase)` | Per-ROI scatter of `value_col` vs temperature, colored by phase. |
| `graph_roi_outlier_detection(df_outliers, roi, value_col, x_col)` | Visualize the local-trend/outlier flags from `detect_roi_signal_outliers`. |
| `graph_temp_range_summary(roi_summary, temp_ranges, value_prefixes)` | Bar/summary plot of the low/mid/high range statistics. |

### `scripts/roi_image_crops.py` — ROI image crops and temperature-sorted views

Independent of the analysis pipeline above (doesn't read/write CSVs); works directly on the `.tif` stack + `RoiSet.zip` of a sample folder.

| Function | Purpose |
|---|---|
| `find_sample_tif(sample_dir)` | Locate the single `.tif`/`.tiff` in a sample folder. |
| `load_rois_ordered(roi_zip_path)` | Load `RoiSet.zip` in ROI-Manager order (must match the `IntDen1..N` column order of the exported CSV). |
| `verify_roi_count(roi_zip_path, expected_n_rois)` | Sanity check: ROI count in the `.zip` vs. distinct ROIs reported by the pipeline. |
| `format_crop_filename(genotype, sample, roi_name, extension)` | Build the lab's crop filename convention, e.g. `m65_i19_R008.tif`. |
| `crop_selected_rois(sample_dir, selected_roi_names, output_dir, genotype, sample, tif_name, roi_zip_name, box_size, save_format)` | Save a `box_size × box_size` multi-frame crop per ROI, all frames in original acquisition order. |
| `crop_included_rois_from_notebook(df_phase, files, sample, genotype, status_col, included_status, roi_zip_name, box_size, save_format, crops_subfolder)` | Shortcut for `01_preprocessing.ipynb`: crop every `ROI_status == 1` ROI right after curation, into the sample's `ROIs/` folder. |
| `list_rois_from_crops_folder(files, crops_subfolder, extension)` | Re-read a crops folder (after manually deleting bad-looking ones) to recover the final ROI list. |
| `sort_stack_by_temperature(stack, frame_temp_df, frame_col, temp_col)` | Reorder a stack's frames from lowest to highest temperature (heating/cooling mixed), dropping frames with no known temperature. |
| `crop_selected_rois_sorted_by_temperature(sample_dir, selected_roi_names, frame_temp_df, output_dir, genotype, sample, tif_name, roi_zip_name, box_size, save_format, frame_col, temp_col, filename_suffix)` | Same as `crop_selected_rois`, but frames are saved in temperature order instead of acquisition order. Used in `04_temperature_tiff.ipynb`. |
| `select_indices_by_temp_bin(temps, bin_width)` | Pick one representative index per `bin_width`-degree temperature bin (closest to the bin center), to sub-sample an over-dense temperature-sorted stack. |
| `montage_from_stack(stack, ncols, padding, pad_value)` | Arrange a stack's frames into a single 2D grid image (via `skimage.util.montage`). |
| `save_montage_png(stack, out_path, ncols, padding)` | Build the montage and save it as an 8-bit PNG, min-max scaled across the whole stack. |
