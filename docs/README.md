# Sphinx documentation source

This directory contains reStructuredText sources for the project manual. **Built HTML is not committed**: output goes to `docs/_build/html/` (listed in the root `.gitignore`).

## Build locally

From the repository root, with the optional docs extras installed:

```bash
pip install -e ".[docs]"
sphinx-build -b html docs docs/_build/html
```

Open `docs/_build/html/index.html` in a browser.

## Published site

Pushes to `main` trigger `.github/workflows/docs.yml`, which builds the same tree and deploys to GitHub Pages (see workflow file for the canonical base URL).
