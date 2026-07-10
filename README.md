# TRPM3 Fluorescence TTL Temperature Analysis

Analysis workflow for TRPM3 fluorescence imaging synchronized with TTL and temperature recordings. The project includes notebooks and reusable Python scripts to preprocess ROI fluorescence traces, align frames with TTL/temperature data, curate active ROIs, quantify responses across temperature ranges, and generate genotype/sample-level plots.

## Project Structure

```text
.
├── environment.yml
├── notebooks/
│   ├── 01_preprocessing.ipynb
│   ├── 02_processing.ipynb
│   └── 03_analysis.ipynb
├── scripts/
│   ├── analisis_ttl.py
│   └── roi_status_selector.py
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

### 3. Analysis

Notebook:

```text
notebooks/03_analysis.ipynb
```

This step:

- compiles every `*_preprocessed_long.csv` under `data/Proc_data/` via `ttl.load_all_preprocessed_long` (no external Excel needed);
- keeps only ROIs already marked as included (`ROI_status == 1`);
- filters by `phase`, `trend`, genotype, and temperature range (parameters at the top of the notebook);
- bins `NormSignal` by temperature per ROI, then summarizes mean ± SEM per genotype;
- plots per-ROI fine curves plus the genotype-level mean curve across temperature.

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
