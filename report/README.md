# Transient Grating Results Report (LaTeX)

This folder contains the LaTeX report associated with the results of the transient grating analysis project.

## Contents

- `main.tex`: Main report source.
- `references.bib`: Bibliography database (enable `biblatex`/`bibliography` in `main.tex` when you cite entries).
- `figures/`: Figures produced by `src/results_tg.ipynb` (or equivalent) and included from LaTeX.
- `logos/`: Logos for the title page.
- `.gitignore`: LaTeX auxiliary files ignored in this folder.

Compiled PDFs such as `main.pdf` are usually kept out of commits unless your group policy requires them; regenerate locally after updating figures or text.

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
