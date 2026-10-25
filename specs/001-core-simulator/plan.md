# Implementation Plan: First-Principles PyWFS, Coronagraph, and Science PSF Data Simulator

**Branch**: `001-core-simulator` | **Date**: 2026-10-25 | **Spec**: [spec.md](file:///Users/seamus/Research/tablasul/specs/001-core-simulator/spec.md)

**Input**: Feature specification from `/Users/seamus/Research/tablasul/specs/001-core-simulator/spec.md`

## Summary

Build `tablasul`, a standalone, first-principles wave-optics simulation library in pure Python and vectorized NumPy. `tablasul` models the end-to-end optical train between an entrance pupil, an evolving atmosphere, a Pyramid Wavefront Sensor (PyWFS), and a Science Focal Plane Camera equipped with a Lyot coronagraph. The primary data output is a synchronized suite of multi-dimensional cubes (`wfs_cube`, `psf_cube`, `ground_truth_phase`) designed specifically for training and benchmarking WFS-based PSF reconstruction algorithms. All physical components feature independent toggles enabling seamless scaling from an unaberrated diffraction-limited "toy model" to a full SCExAO-fidelity simulation.

## Technical Context

**Language/Version**: Python 3.11+  
**Primary Dependencies**:
- `numpy` (>=1.24): Core vectorized array math, complex electric fields, 2D FFTs
- `scipy` (>=1.10): Fast interpolation, special functions (Bessel functions for Airy verification)
- `astropy` (>=5.0): Astronomical FITS file I/O with standard header keywords
- `pytest` (>=7.0): Test suite and analytical regression tests

**Target Platform**: macOS (Apple Silicon / Intel), Linux, scalable to HPC clusters.  
**Project Type**: Python library + CLI dataset generation tool.  
**Performance Goals**:
- Toy mode ($128\times 128$ grid, unmodulated): $<50$ ms per frame pair on MacBook CPU.
- Modulated mode ($128\times 128$ grid, $N_{\text{mod}} = 16$): $<500$ ms per frame pair on MacBook CPU.
- 1000-frame synchronized cube generation: $<60$ seconds in toy mode.

**Constraints**:
- First-principles Fourier optics: zero dependencies on black-box monolithic simulators.
- Strict physical units throughout: meters, radians, wavelengths $\lambda$, flux in photons.
- Independent toggles for all 6 subsystems (Atmosphere, Obscuration, NCPA, Coronagraph, Modulation, Noise).

## Constitution Check

- [x] **First-Principles Wave Optics**: Uses paraxial Fraunhofer 2D FFT with centered coordinates ($n // 2$ DC index).
- [x] **Analytical Benchmarks Included**: Unit tests check Airy disk intensity within $<10^{-3}$ and Parseval energy conservation $>99.9\%$.
- [x] **Modular & Toggleable**: Every optical class accepts `enabled: bool` and defaults to bypass mode when disabled.
- [x] **Synchronized Cubes**: Data outputs are multi-dimensional arrays aligned in time and wavelength.
- [x] **No Heavy Black-Box Dependencies**: Pure NumPy/SciPy core.

## Project Structure

```text
/Users/seamus/Research/tablasul/
├── pyproject.toml
├── README.md
├── src/
│   └── tablasul/
│       ├── __init__.py
│       ├── core/
│       │   ├── __init__.py
│       │   ├── field.py          # ComplexElectricField: Fraunhofer propagation, FFT shifts, energy conservation
│       │   ├── aperture.py       # CircularAperture, SubaruAperture (spiders, central obscuration)
│       │   ├── source.py         # PointSource, CompanionSource, CalibrationSource, FilterBand
│       │   └── atmosphere.py     # KolmogorovScreen, FrozenFlowAtmosphere (dynamic phase shifting)
│       ├── optics/
│       │   ├── __init__.py
│       │   ├── pywfs.py          # PyramidWFS: 4-facet phase mask, circular modulation, raw pupils & slopes
│       │   ├── coronagraph.py    # LyotCoronagraph: focal plane mask (FPM), Lyot stop, bypass logic
│       │   └── ncpa.py           # ZernikeNCPA, StaticPhaseNCPA
│       ├── detector/
│       │   ├── __init__.py
│       │   └── camera.py         # DetectorModel: integration, Poisson shot noise, Gaussian read noise
│       ├── engine/
│       │   ├── __init__.py
│       │   ├── simulator.py      # OpticalSystem: pipeline graph orchestrator with component toggles
│       │   └── sync.py           # SynchronizedCubeEngine: temporal oversampling & cube assembly
│       └── io/
│           ├── __init__.py
│           └── export.py         # FITS, HDF5, and NPZ serializers with complete header metadata
└── tests/
    ├── test_field_propagation.py # Analytical Airy disk, Parseval energy conservation
    ├── test_aperture.py          # Clear pupil vs Subaru pupil geometry
    ├── test_pywfs.py             # 4-quadrant balance, modulation linearity (R^2 > 0.99)
    ├── test_coronagraph.py       # On-axis starlight rejection, companion transmission, bypass mode
    ├── test_atmosphere.py        # Kolmogorov variance, frozen-flow phase drift
    ├── test_synchronization.py   # Spatiotemporal cube shapes, temporal oversampling
    └── test_toggles.py           # 64-combination toggle test asserting zero crashes and stable shapes
```
