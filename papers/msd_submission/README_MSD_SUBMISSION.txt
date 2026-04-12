This directory contains a flat LaTeX submission package for Multibody System Dynamics.

Contents
- `main.tex`: main manuscript source
- `00_abstract_en.tex` to `07_declarations_en.tex`: section files
- `references.bib`, `references_extra.bib`: bibliography databases
- all figure PDF files used by the manuscript

Build
- `latexmk -pdf -interaction=nonstopmode -file-line-error main.tex`

Notes
- The package is flattened to avoid subdirectories during LaTeX submission.
- The manuscript PDF compiled from this directory is intended to match the working paper under `papers/paper en`.
