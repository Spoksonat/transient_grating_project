# Source Code

This folder contains the analysis implementation and interactive notebooks.

## Files

- `tg_analysis.py`: Main Python module with class `TGAnalysis`.
- `results_tg.ipynb`: Notebook for exploratory analysis and plotting.

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

## Import note

When running notebooks from this directory, use relative paths such as:

- `../external_files/parameters_CoOx12.json`

to access metadata and reference files.

