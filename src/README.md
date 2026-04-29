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
2. Define analysis mode (`constant_E` or `constant_I`) and target fit parameter.
3. Inspect phase space (`plot_phase_space`, optionally with error bars).
4. Select scans (`get_data_scan`) according to the chosen mode.
5. Fit one model or per-scan models (`get_fit_parameters`).
6. Visualize fits (`plot_fits`).
7. Plot any fitted parameter vs energy (`plot_params_vs_energy`) or vs intensity (`plot_params_vs_intensity`), depending on mode.
8. Compare models (`plot_params_all_models`) with explicit mode selection.

## Import note

When running notebooks from this directory, use relative paths such as:

- `../external_files/parameters_CoOx12.json`

to access metadata and reference files.

