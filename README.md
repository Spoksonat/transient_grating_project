# Transient Grating Project

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
![Jupyter Notebook](https://img.shields.io/badge/Jupyter-Notebook-F37626?logo=jupyter&logoColor=white)

Python workflow for transient grating (TG) analysis from scan metadata and MATLAB outputs. The project is centered on `src/tg_analysis.py` and mirrored in `src/results_tg.ipynb`.

## Repository layout

```text
transient_grating_project/
├── README.md
├── pyproject.toml
├── requirements.txt
├── .github/workflows/
│   └── docs.yml          # Sphinx → GitHub Pages on push to main
├── docs/                  # Sphinx sources (see docs/README.md)
├── report/                # LaTeX report (see report/README.md)
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

`report/` holds the LaTeX companion article (`main.tex`, figures, bibliography). HTML API documentation is generated under `docs/` and built output is kept in `docs/_build/` (ignored by git).

## What changed in the current analysis workflow

The current code/notebook workflow includes:

- TG signal preprocessing with baseline subtraction and normalization by total area (`signal / sum(signal)`).
- Mode-based scan selection at fixed intensity or fixed energy from JSON metadata.
- Three fit models with parameter extraction and propagated uncertainties.
- Generic parameter-vs-energy plotting via `plot_params_vs_energy(param_name=...)`.
- Parameter-vs-intensity plotting via `plot_params_vs_intensity(param_name=...)`.
- Multi-model comparison via `plot_params_all_models(..., mode=...)`.

## Data contract

The metadata file `external_files/parameters_CoOx12.json` must include:

- `main_path`: folder that contains `ScanXXX.mat` files.
- one `ScanXXX` object per scan.
- fields used by `TGAnalysis`: `Scan`, `Energy_eV`, `XUV_intensity_uJ`, `flag_secure`, `repeated`, `filter`.

The `.mat` files are expected to contain the MATLAB structure under key `R` used in `get_data_scan`.

## Quick start

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Example usage (same logic as notebook):

```python
from tg_analysis import TGAnalysis

analysis = TGAnalysis("external_files/parameters_CoOx12.json")
mode = "constant_I"
param_name = "tau"

analysis.plot_phase_space(errors_bool=True)

if mode == "constant_E":
    analysis.get_data_scan({"E": 63.6, "I": "all"})
else:
    analysis.get_data_scan({"E": "all", "I": 2.0})

analysis.get_fit_parameters(model_idxs=2, initial_guess_bool=True, bounds=True)

analysis.plot_fits()
if mode == "constant_E":
    analysis.plot_params_vs_intensity(param_name=param_name, errors_bool=False)
else:
    analysis.plot_params_vs_energy(param_name=param_name, errors_bool=False)
```

## Documentation

Build the Sphinx manual locally:

```bash
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e ".[docs]"
sphinx-build -b html docs docs/_build/html
```

Open `docs/_build/html/index.html`. The build directory `docs/_build/` is **gitignored** by design; regenerate it whenever sources change. Pushes to `main` rebuild and deploy HTML via GitHub Actions (`.github/workflows/docs.yml`).

More detail: `docs/README.md`.

## External files and Git tracking

`external_files/parameters_CoOx12.json` is intentionally tracked in git and must be committed for reproducible scan selection and metadata provenance.

## Notes

- Plotting is configured with `matplotlib` LaTeX rendering (`text.usetex=True`), so a local LaTeX installation is required.
- Some fit/plot methods assume `get_data_scan(...)` and `get_fit_parameters(...)` have already been executed.

## License

Distributed under the MIT License. See `LICENSE`.