# Source Code

This folder contains the analysis implementation and interactive notebooks.

## Files

- `tg_analysis.py`: Main Python module with class `TGAnalysis`.
- `results_tg.ipynb`: Notebook that executes the end-to-end workflow used in the current analysis.

## Responsibilities of this folder

- Store reusable analysis logic in Python modules.
- Keep notebooks as analysis drivers and experiment logs.
- Avoid placing raw external data here (use `../external_files/`).

## Suggested coding pattern

- Add reusable functions/methods in `tg_analysis.py`.
- Keep notebook cells focused on:
  - loading parameters,
  - selecting scans,
  - plotting/fit diagnostics,
  - saving summarized outputs.

## Current notebook flow (`results_tg.ipynb`)

1. Load metadata JSON and initialize `TGAnalysis`.
2. Inspect phase space (`plot_phase_space`).
3. Select scans (`get_data_scan`).
4. Fit one model or per-scan models (`get_fit_parameters`).
5. Visualize fits (`plot_fits`).
6. Plot any fitted parameter vs energy (`plot_params_vs_energy`).
7. Compare models (`plot_params_all_models`).

## Import note

When running notebooks from this directory, use relative paths such as:

- `../external_files/parameters_CoOx12.json`

to access metadata and reference files.

