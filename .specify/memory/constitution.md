# Tablasul Constitution
<!-- Project Constitution for Tablasul: First-Principles PyWFS, Coronagraph, and PSF Data Simulator -->

## Core Principles

### I. First-Principles Physical & Mathematical Rigor
All optical propagations and transformations must be mathematically and physically grounded in wave optics:
- **Fourier Optics Foundation**: Optical field propagation must be implemented from first principles using 2D Fourier transforms (Fraunhofer / Angular Spectrum) operating on complex electric fields $E(\vec{r}) = A(\vec{r}) \exp(i\phi(\vec{r}))$.
- **Explicit Physical Units**: All quantities must have clear, consistent physical units (apertures in meters, wavelengths $\lambda$ in meters, wavefront phase in radians or nanometers of optical path difference (OPD), flux in photons/frame or photons/s).
- **Sampling & Anti-Aliasing**: Spatial grid dimensions and sampling parameters must satisfy the Nyquist criterion (at least 2 pixels per $\lambda/D$ resolution element) to prevent numerical aliasing artifacts in both focal and pupil planes. Energy conservation (Parseval's theorem) must hold across all lossless optical planes.

### II. Composable, Toggleable Architecture ("Toy to Full Scale")
The simulator must accommodate both rapid prototyping ("toy model") and high-fidelity astronomical simulation:
- **Individual Component Toggles**: Every optical and physical element—atmosphere, pupil spiders/obscurations, NCPA, coronagraph mask, Lyot stop, PyWFS modulation, and detector noise—must be independently toggleable (`enabled: bool`) and individually parameterizable.
- **Graceful Baseline**: At minimal configuration (`all toggles disabled`), the engine collapses to a fast, clean, diffraction-limited baseline: clear circular aperture, unmodulated PyWFS, and direct Airy disk PSF without noise.
- **Zero API Churn**: Toggling features on or off must not break the data pipeline interface or change output tensor ranks.

### III. Synchronized Spatiotemporal Data Cubes
The fundamental output contract of Tablasul is producing synchronized multi-dimensional data arrays:
- **Paired Data Streams**: The simulator must generate synchronized streams containing:
  1. `wfs_cube`: PyWFS telemetry (raw 4-pupil intensity images $[T, 4, H_{\text{pup}}, W_{\text{pup}}]$ and/or normalized slopes $[T, 2, H_{\text{pup}}, W_{\text{pup}}]$).
  2. `psf_cube`: Science focal plane images $[T, C, H_{\text{psf}}, W_{\text{psf}}]$.
  3. `ground_truth_phase`: Instantaneous and exposure-integrated wavefront phase screens $[T, H_{\text{pup}}, W_{\text{pup}}]$ at the entrance pupil and focal-plane electric field.
- **Temporal Fidelity**: The temporal engine must accurately capture the relationship between atmospheric turbulence evolution (e.g. wind velocity, frozen flow time steps $dt_{\text{atm}}$), WFS frame rate $f_{\text{wfs}}$, and science camera integration time $\tau_{\text{psf}}$ (including multi-substep integration for long-exposure PSFs).

### IV. Generalized Source & Spectral Mechanics
Optical fields must handle diverse astronomical and bench calibration scenarios:
- **Target Types**: Support stellar sources (blackbody spectra, magnitude-based flux), off-axis point sources (exoplanet companions with specified separation in $\lambda/D$, contrast ratio, and position angle), and internal calibration sources (monochromatic calibration lasers, flat lamps).
- **Wavelength & Filter Separation**: Enable distinct spectral bands between the WFS path (e.g. visible $R/I$ band) and Science path (e.g. NIR $J/H$ band). Broadband sources must be integrated across discrete monochromatic channels weighted by standardized or user-defined transmission curves.

### V. Clean First-Principles Python & Scale Readiness
- **Zero Heavy Dependencies**: Core optical propagation must be self-contained in pure, vectorized NumPy and SciPy without requiring heavy external astronomical frameworks.
- **Vectorized Array Operations**: Code must be written with vectorized NumPy operations and clean functional/class boundaries to enable straightforward future acceleration (CuPy, PyTorch, JAX) and multi-core scaling.
- **Export Standards**: Data products must be exportable to standard scientific formats (NumPy `.npy`/`.npz`, HDF5, and astronomical FITS with comprehensive FITS header metadata).

## Optical & Technical Standards

1. **PyWFS Physics**:
   - Circular tip-tilt modulation simulated by integrating discrete pupil tilts on a circle of radius $r_{\text{mod}}$ in units of $\lambda_{\text{wfs}}/D$. Setting $r_{\text{mod}} = 0$ yields static, unmodulated pyramid sensing.
   - 4-sided pyramid modeled as a 4-quadrant phase mask in the focal plane.
   - Dual output capability: both raw 4-pupil detector images ($I_1, I_2, I_3, I_4$) and normalized quad-cell slopes ($S_x, S_y$).
2. **Coronagraph Physics**:
   - Classical / Apodized Lyot Coronagraph: Entrance pupil $\to$ Focal Plane Mask (opaque central disk with user-defined IWA radius in $\lambda/D$) $\to$ Lyot Stop plane (undersized circular aperture with optional spider masking) $\to$ Science focal plane.
   - Bypass mode: When disabled, the beam bypasses the focal mask and Lyot stop directly to the science detector.
3. **Atmosphere & Aberration Budget**:
   - Kolmogorov / von Kármán turbulence phase screens with Fried parameter $r_0$, outer scale $L_0$, and wind vector $(v_x, v_y)$.
   - Distinct Non-Common Path Aberration (NCPA) phase map injected specifically into the science arm or WFS arm to simulate bench calibration offsets.
4. **Detector Noise Model**:
   - Photon shot noise (Poisson distribution).
   - Read noise (zero-mean Gaussian distribution with standard deviation $\sigma_{\text{RON}}$).
   - Dark current and background photon counts.

## Quality & Verification Gates

- **Analytical Verification**:
  - Unaberrated direct PSF must match the analytical Airy disk pattern to within $<10^{-4}$ relative intensity error.
  - Pyramid slope response under small tip-tilt aberrations must match analytical linear slope models.
  - Energy conservation: total flux in pupil plane must equal total flux in unobstructed focal plane (Parseval theorem).
- **Test Coverage**: All core optical transformations, modulation routines, filter weightings, and noise models must have automated unit tests (`pytest`).

## Governance

- The Constitution defines the non-negotiable architectural and physical guidelines for the `tablasul` package.
- Any pull request or new feature must adhere to the modular toggle pattern, physical unit consistency, and test verification standards.

**Version**: 1.0.0 | **Ratified**: 2026-10-25 | **Author**: SCExAO / Tablasul Team
