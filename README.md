# TRPM3 Fluorescence TTL Temperature Analysis

Analysis workflow for TRPM3 fluorescence imaging synchronized with TTL and temperature recordings. The project includes notebooks and reusable Python scripts to preprocess ROI fluorescence traces, align frames with TTL/temperature data, curate active ROIs, quantify responses across temperature ranges, and generate genotype/sample-level plots.

## Project Structure

```text
.
├── environment.yml
├── notebooks/
│   ├── 01_preprocessing.ipynb
│   ├── 02_processing.ipynb
│   └── 03_graphs.ipynb
├── scripts/
│   ├── analisis_ttl.py
│   ├── roi_status_selector.py
│   └── graphs.py
├── data/              # local data, not tracked by git
├── graphs/            # optional exported figures
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

Each sample folder should contain the ABF file, ROI results, and image-related files needed by the preprocessing notebook.

Example:

```text
data/Proc_data/
├── TRPM3_imagenes.xlsx
├── mut27_image10/
├── mut36_image8/
├── mut57_image2/
└── mut65_image19/
```

Large data files are intentionally ignored by git.

## Workflow

### 1. Preprocessing

Notebook:

```text
notebooks/01_preprocessing.ipynb
```

This step:

- loads ABF temperature/TTL recordings;
- detects TTL events;
- aligns ImageJ ROI fluorescence frames with TTL/temperature;
- normalizes fluorescence traces;
- assigns temperature phase and trend labels;
- exports `*_preprocessed_long.csv` files for each sample.

### 2. ROI Processing and Curation

Notebook:

```text
notebooks/02_processing.ipynb
```

This step:

- loads active experiments from `TRPM3_imagenes.xlsx`;
- filters data by genotype, phase, and trend;
- applies or updates ROI inclusion/exclusion status;
- uses an interactive Matplotlib selector to curate ROIs;
- saves a cumulative `roi_status_registry.csv`;
- exports the main curated table:

```text
data/Proc_data/batch_analysis/roi_temp_summary_active.csv
```

This is the main output for downstream plotting and analysis.

### 3. Graphs

Notebook:

```text
notebooks/03_graphs.ipynb
```

This step:

- loads `roi_temp_summary_active.csv`;
- computes the high/low response ratio as:

```text
(1 + high_mean) / (1 + low_mean)
```

Because `low_mean` and `high_mean` are normalized fluorescence signals (`DeltaF/F0`), this ratio compares estimated `F/F0` between high and low temperature ranges.

The notebook generates:

- boxplots of high/low ratio by trend, one figure per genotype;
- low vs high response scatter plots, colored by sample and shaped by trend;
- ratio summary CSV files and exported figures.

### 4. Fine-Grained ROI Metrics

Script:

```text
scripts/fine_analysis.py
```

This optional step builds a more detailed ROI profile from
`data/Proc_data/batch_analysis/processing_active.csv`:

- response curves in small temperature bins;
- ROI-level curve metrics such as signal range, slope, AUC, max/min response
  temperature, high-vs-low temperature delta, and threshold-crossing temperature;
- conservative QC flags;
- sample/genotype summaries.

Run it from the project root:

```bash
python scripts/fine_analysis.py
```

It exports:

```text
data/Proc_data/batch_analysis/fine_analysis/
├── filtered_points_used.csv
├── roi_temperature_bins.csv
├── roi_curve_metrics.csv
└── roi_curve_metrics_by_group.csv
```

## Main Outputs

Important files produced under:

```text
data/Proc_data/batch_analysis/
```

include:

```text
roi_status_registry.csv
roi_temp_summary_active.csv
roi_temp_summary_long.csv
roi_high_low_ratio_*.csv
roi_high_low_ratio_summary_*.csv
figures/
```

## Notes

- `data/`, `raw/`, and `processed/` are ignored by git to avoid uploading large local files.
- The interactive ROI selector uses a local Matplotlib Qt window.
- Most analysis logic lives in `scripts/` so notebooks stay focused on workflow and parameter choices.
