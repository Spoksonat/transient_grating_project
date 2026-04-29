# Transient Grating Results Report (LaTeX)

This folder contains the LaTeX report associated with the results of the transient grating analysis project.

## Contents

- `main.tex`: Main report source.
- `references.bib`: Bibliography database.
- `figures/`: Figures used in the report.
- `logos/`: Optional logos for title pages.
- `.gitignore`: LaTeX auxiliary files ignored in this folder.

## Build

From this folder:

```bash
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

Recommended (automatic dependency tracking):

```bash
latexmk -pdf -interaction=nonstopmode main.tex
```

## Notes

- Write one sentence per line in running text to keep diffs clean.
- Keep figure sources in `figures/` and reference them with relative paths.
- Add citations in `main.tex` with `~\cite{...}` and keep entries in `references.bib`.
