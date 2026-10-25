# Implementation Plan: Refine PyWFS Pupil Separation, Coronagraph Contrast, and Documentation

**Branch**: `002-optics-refinement` | **Date**: 2026-10-25 | **Spec**: [spec.md](file:///Users/seamus/Research/tablasul/specs/002-optics-refinement/spec.md)

## Summary

1. Refine `PyramidWFS` separation physics: increase default `separation_ratio` to `2.0` and optimize Fourier grid allocation so that the 4 circular pupils are well-separated with wide dark margins matching standard camera displays.
2. Refine `LyotCoronagraph` transmission profile: support both standard hard-edge and smooth apodized FPMs (`edge_softness`), and update visualization tools to display true astronomical normalized contrast units ($I / I_{\text{direct\_peak}}$) on a logarithmic scale from $10^{-6}$ to $1$.
3. Git & Documentation: Initialize a clean Git repository, establish `.gitignore`, write an extensive documentation manual in `docs/` (optical derivations, coordinate systems, synchronization guides), and update `README.md`.

## Changes Planned

- `src/tablasul/optics/pywfs.py`: Set default `separation_ratio = 2.0`, update quadrant window bounding box calculation.
- `src/tablasul/optics/coronagraph.py`: Add `edge_taper_lod: float = 0.5` (cosine/error-function smooth taper option).
- `examples/generate_example_dataset.py`: Update visualization to plot normalized astronomical contrast.
- `.gitignore`: Ignore `.venv/`, `__pycache__/`, `.pytest_cache/`, `*.egg-info/`, `.specify/feature.json`.
- `docs/optical_theory.md` and `docs/api.md`: Detailed documentation.
