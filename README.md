# Transient Grating Project

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
![Jupyter Notebook](https://img.shields.io/badge/Jupyter-Notebook-F37626?logo=jupyter&logoColor=white)
![Status](https://img.shields.io/badge/status-active-success)

Repository for transient grating (TG) signal analysis, including:

- scan metadata management from JSON files,
- loading and preprocessing of experimental `.mat` scans,
- filtering and normalization of TG traces,
- exploratory plotting for phase-space coverage and stacked signals.

This project is organized to keep analysis code and external experimental inputs clearly separated.

## Repository layout

```text
transient_grating_project/
├── README.md
├── .gitignore
├── src/
│   ├── README.md
│   ├── tg_analysis.py
│   └── results_tg.ipynb
└── external_files/
    ├── README.md
    ├── parameters_CoOx12.json
    ├── Co2plus_absorbance.txt
    ├── Co3plus_absorbance.txt
    └── Co3O4_absorbance.txt
```

## Project goals

- Provide a reusable Python analysis class for TG scans.
- Keep metadata and external reference files in a dedicated data folder.
- Make quick visual diagnostics easy (phase-space maps and stacked traces).
- Support later extension to fit models and publication-grade plots.

## Current analysis module

The main analysis implementation is in `src/tg_analysis.py` with class `TGAnalysis`.

Main capabilities currently available:

- Load scan metadata from a JSON parameter file.
- Build a scan table (`pandas.DataFrame`) from all `ScanXXX` entries.
- Filter scans by:
  - constant energy (`get_data_E_const`),
  - constant intensity (`get_data_I_const`).
- Load scan signals from MATLAB files (`scipy.io.loadmat`) using `main_path`.
- Apply optional time filtering (for scans with discontinuities).
- Normalize TG signals and align time to signal maximum.
- Plot:
  - phase-space points (`Energy_eV`, `XUV_intensity_uJ`),
  - stacked TG traces with scan labels and experimental conditions.

## Data and metadata contract

The JSON file in `external_files/` is expected to contain:

- top-level `main_path` pointing to the folder with scan `.mat` files,
- one object per scan (for example `Scan038`, `Scan039`, ...),
- scan fields used by the code:
  - `Scan`,
  - `Energy_eV`,
  - `XUV_intensity_uJ`,
  - `flag_secure`,
  - `repeated`,
  - `filter` (used when filtering time traces).

Example workflow assumptions:

- `main_path + f"{Scan}.mat"` exists for each selected scan.
- secure and non-repeated scans are selected by default in analysis filters.

## Getting started

### 1. Clone and enter project

```bash
git clone <your-repo-url>
cd transient_grating_project
```

### 2. Create a Python environment

```bash
python -m venv .venv
source .venv/bin/activate
pip install numpy pandas scipy matplotlib jupyter
```

### 3. Prepare external files

- Put metadata JSON and reference text files in `external_files/`.
- Ensure `main_path` in JSON points to your local `.mat` scan directory.

### 4. Run analysis

Use either:

- `src/results_tg.ipynb` for interactive exploration, or
- a Python script importing `TGAnalysis` from `src/tg_analysis.py`.

## Recommended development workflow

- Keep reusable logic in `src/tg_analysis.py`.
- Keep exploratory steps and temporary plots in notebooks.
- When a notebook block becomes stable, migrate it to `.py` code.
- Use clear scan-selection criteria and document them in notebook cells.

## External data policy

- `external_files/` is intended for lightweight metadata and reference files.
- Large raw data files should stay outside git when possible.
- If large data must be tracked, prefer Git LFS and document it explicitly.

## Notes and limitations

- Plot rendering currently uses LaTeX text mode (`text.usetex = True`), which
requires a working LaTeX installation on your system.
- Some methods assume specific MATLAB internal structure under key `R`.
- The current implementation is focused on TG preprocessing/visualization;
dedicated model-fitting interfaces can be added incrementally.

## Next improvements (suggested)

- Add a requirements file (`requirements.txt` or `pyproject.toml`).
- Add tests for JSON parsing and scan selection behavior.
- Add explicit error handling for missing scan files and malformed metadata.
- Add fit models as class methods with uncertainty estimation.

## License

This project is distributed under the MIT License.
See `LICENSE` for the full text.

