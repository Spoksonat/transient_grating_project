# External Files

This folder stores metadata and external reference inputs used by the analysis code.

## Current contents

- `parameters_CoOx12.json`: Scan metadata and `main_path` to MATLAB scan files.
- `Co2plus_absorbance.txt`: Reference absorbance spectrum.
- `Co3plus_absorbance.txt`: Reference absorbance spectrum.
- `Co3O4_absorbance.txt`: Reference absorbance spectrum.

## Usage in code

`src/tg_analysis.py` reads `parameters_CoOx12.json` to:

- locate scan files through `main_path`,
- select scans by energy/intensity,
- filter by quality/repetition flags.

## Conventions for new files

- Keep file names descriptive and stable.
- Use plain text formats when possible (`.json`, `.txt`, `.csv`).
- Document units and provenance in file headers or companion notes.

## Version-control policy

- Keep lightweight metadata and small reference files in git.
- Do not commit large raw datasets by default.
- If sharing large data is required, prefer Git LFS or external storage.
