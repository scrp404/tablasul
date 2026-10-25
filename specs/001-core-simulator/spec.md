# Feature Specification: First-Principles Synchronized PyWFS, Coronagraph, and Science PSF Data Simulator

**Feature Branch**: `001-core-simulator`  
**Created**: 2026-10-25  
**Status**: Draft  
**Project**: Tablasul (`/Users/seamus/Research/tablasul`)  
**Input**: First-principles physical/mathematical data simulator producing synchronized cubes of Pyramid Wavefront Sensor (PyWFS) and Science Focal Plane Camera (PSF) data with toggleable coronagraph, atmosphere, modulation, sources/filters, and detector noise.

---

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Minimal Diffraction-Limited Benchmark ("Toy Model") (Priority: P1)

As an adaptive optics researcher simulating synchronized WFS and PSF data,  
I want to simulate a baseline, unaberrated, diffraction-limited optical system (clear aperture, unmodulated PyWFS, direct science focal plane, no noise),  
So that I can verify that the first-principles Fourier optics engine produces an analytical Airy disk and correct 4-pupil intensity distributions with zero numerical artifacts before adding complexity.

**Why this priority**: Establishes the foundational wave-optics truth. If the basic unaberrated Fraunhofer propagation and pyramid prism phase splitting fail analytical benchmarks, downstream algorithms cannot evaluate data reliably.

**Independent Test**:
Can be fully tested by generating a single monochromatic frame ($N_{\text{frames}} = 1$) on a clear circular pupil, asserting that the science image matches the theoretical Airy pattern to within $<10^{-3}$ relative intensity error, and checking that the unmodulated pyramid produces 4 identical pupil quadrants of equal integrated flux.

**Acceptance Scenarios**:
1. **Given** a clear circular pupil of diameter $D$ and monochromatic light $\lambda$, **When** simulated in direct imaging mode with all optical toggles disabled, **Then** the science focal plane produces an exact Airy pattern centered on the optical axis with total flux conserved ($>99.9\%$ energy conservation).
2. **Given** an unmodulated PyWFS ($r_{\text{mod}} = 0$) and zero aberration, **When** the wavefront passes through the 4-faceted pyramid, **Then** 4 distinct pupil quadrant images are generated on the WFS detector with identical spatial profile and equal integrated counts ($I_1 = I_2 = I_3 = I_4$).

---

### User Story 2 - Composable Optical Train & Lyot Coronagraph (Priority: P2)

As an AO researcher,  
I want to toggle and configure realistic telescope pupil geometry (Subaru 8.2m aperture with central obscuration and 4 spider vanes) and a Lyot coronagraph (focal plane mask + Lyot stop),  
So that I can study how pupil obscurations and coronagraphic diffraction suppression alter the relationship between WFS signals and science focal plane speckles.

**Why this priority**: Telescope systems operate behind obscured apertures with spiders and use coronagraphs to suppress on-axis starlight. Simulating focal plane data behind a coronagraph requires explicit modeling of the focal plane mask and Lyot stop.

**Independent Test**:
Can be tested by comparing direct science PSF vs. coronagraphic science PSF for an on-axis star: the Lyot coronagraph must suppress the on-axis peak intensity by the expected theoretical rejection factor while preserving an off-axis companion.

**Acceptance Scenarios**:
1. **Given** a Subaru pupil geometry (central obscuration ratio $\approx 0.29$ and 4 spider struts), **When** enabled, **Then** the pupil transmission mask correctly modulates the complex electric field, and diffraction spikes appear in the science focal plane.
2. **Given** an enabled Lyot coronagraph with an opaque focal plane mask of radius $r_{\text{FPM}}$ (e.g. $2.5\,\lambda/D$) and an undersized Lyot stop (e.g. 90% pupil diameter), **When** light propagates through the science arm, **Then** the on-axis starlight is rejected at the Lyot plane, and residual light in the science focal plane exhibits a dark hole / coronagraphic attenuation profile.
3. **Given** the coronagraph toggle is set to `enabled=False`, **When** the simulation executes, **Then** the focal plane mask and Lyot stop are cleanly bypassed, returning a direct unocculted PSF without reconfiguring the pipeline.

---

### User Story 3 - PyWFS Tip-Tilt Modulation & Dual Output Telemetry (Priority: P3)

As an AO modeler,  
I want to configure the circular tip-tilt modulation radius of the PyWFS ($r_{\text{mod}}$ in $\lambda_{\text{wfs}}/D$) and retrieve both raw 4-pupil detector images and normalized slope maps ($S_x, S_y$),  
So that I can analyze WFS telemetry using either raw image data or classical wavefront slope estimates.

**Why this priority**: Astronomical PyWFS instruments typically modulate the PyWFS (e.g., $1.5\text{--}3\,\lambda/D$) to increase dynamic range and linearity. Testing models against both modulated and unmodulated WFS data is crucial for assessing sensitivity to modulation errors.

**Independent Test**:
Can be tested by applying a known small tip/tilt phase tilt across the entrance pupil, varying modulation radius from $r_{\text{mod}} = 0$ to $r_{\text{mod}} = 3\,\lambda/D$, and verifying that the calculated slope signals ($S_x, S_y$) respond linearly to input tilt with slope gain scaling inversely with modulation radius.

**Acceptance Scenarios**:
1. **Given** a user setting for modulation radius $r_{\text{mod}} > 0$ and $N_{\text{mod}}$ modulation points, **When** simulated, **Then** the system integrates the pyramid focal plane phase mask over a circular path of radius $r_{\text{mod}}$, producing modulated 4-pupil detector images.
2. **Given** raw 4-pupil images $(I_1, I_2, I_3, I_4)$, **When** slope calculation is requested, **Then** the simulator outputs normalized quad-cell slopes:
   $$S_x = \frac{(I_1 + I_4) - (I_2 + I_3)}{I_1 + I_2 + I_3 + I_4}, \quad S_y = \frac{(I_1 + I_2) - (I_3 + I_4)}{I_1 + I_2 + I_3 + I_4}$$
   masked to the valid illuminated pupil subapertures.
3. **Given** $r_{\text{mod}} = 0$, **When** simulated, **Then** the engine runs in static unmodulated mode with maximum sensitivity.

---

### User Story 4 - Multi-Source & Spectral Filter Modeling (Priority: P4)

As an optical experimenter,  
I want to configure diverse source types (astronomical stars with blackbody/filter spectra, off-axis planetary companions, or internal monochromatic calibration lasers and flat lamps), with distinct spectral bands for the WFS arm and Science arm,  
So that I can simulate realistic chromatic differences between the WFS sensing wavelength and science imaging wavelength, as well as bench calibration runs.

**Why this priority**: Typically, the PyWFS senses in the visible ($R/I$ band, e.g. 600–900 nm) while science cameras operate in the near-infrared ($J/H/K$ band, 1.1–1.8 $\mu$m). Wavefront phase scales inversely with wavelength ($\phi \propto 1/\lambda$), creating chromatic speckle patterns that simulations must account for.

**Independent Test**:
Can be tested by simulating a calibration laser (single monochromatic line $\lambda = 635$ nm in both arms) vs an astronomical target (visible WFS $\lambda_{\text{wfs}} = 800$ nm, NIR science $\lambda_{\text{sci}} = 1600$ nm). Verifying that an optical path difference $\Delta \text{OPD}$ creates phase $\phi_{\text{wfs}} = 2\pi \Delta \text{OPD}/\lambda_{\text{wfs}}$ and $\phi_{\text{sci}} = 2\pi \Delta \text{OPD}/\lambda_{\text{sci}}$ (a factor of 2 difference).

**Acceptance Scenarios**:
1. **Given** an astronomical star with a specified filter band (e.g. $R$-band for WFS, $H$-band for Science), **When** simulated, **Then** the electric fields in both arms are evaluated at their respective central wavelengths and integrated across discrete filter spectral channels.
2. **Given** an off-axis companion with contrast ratio $C = 10^{-4}$, radial separation $\rho = 5\,\lambda/D$, and position angle $\theta = 45^\circ$, **When** simulated, **Then** an off-axis point source appears in the science PSF at the correct position and contrast.
3. **Given** a calibration source setting (e.g. monochromatic internal laser or flat lamp), **When** selected, **Then** the simulator configures a clean monochromatic calibration beam with zero atmospheric turbulence.

---

### User Story 5 - Synchronized Spatiotemporal Data Cubes (Priority: P5)

As an ML/AI model developer,  
I want to generate synchronized multi-dimensional data cubes of WFS telemetry, Science PSFs, and Ground-Truth wavefront phase screens, with user-defined frame counts and temporal cadences,  
So that I have perfectly synchronized paired training and evaluation sets where every science frame corresponds precisely to the concurrent WFS state(s).

**Why this priority**: Simulating synchronized WFS and PSF telemetry requires time-aligned data. Because science cameras often integrate over longer exposure times than fast WFS loops (e.g. 1 kHz WFS vs. 10 Hz science camera), the simulator must handle oversampling where multiple WFS frames integrate into one science PSF.

**Independent Test**:
Can be tested by generating a cube of $T_{\text{wfs}} = 100$ frames with an exposure ratio of 10:1 (1 science frame per 10 WFS frames), verifying that the science cube has $T_{\text{psf}} = 10$ frames, each equal to the average of the 10 corresponding instantaneous PSF realizations.

**Acceptance Scenarios**:
1. **Given** parameters $N_{\text{frames}}$, WFS frame rate $f_{\text{wfs}}$, and science exposure time $\tau_{\text{psf}}$, **When** the simulation executes over an evolving Kolmogorov phase screen, **Then** the engine outputs:
   - `wfs_cube`: Shape $[T_{\text{wfs}}, 4, H_{\text{wfs}}, W_{\text{wfs}}]$ (or $[T_{\text{wfs}}, 2, H_{\text{pup}}, W_{\text{pup}}]$ for slopes).
   - `psf_cube`: Shape $[T_{\text{psf}}, C_{\text{sci}}, H_{\text{psf}}, W_{\text{psf}}]$.
   - `ground_truth_phase`: Shape $[T_{\text{wfs}}, H_{\text{pup}}, W_{\text{pup}}]$ in radians and nanometers OPD.
2. **Given** an evolving frozen-flow atmosphere with wind velocity $\vec{v} = (v_x, v_y)$, **When** stepping forward in time, **Then** the phase screen shifts smoothly across the pupil aperture by $\Delta \vec{r} = \vec{v} \cdot dt$.

---

### User Story 6 - Non-Common Path Aberrations (NCPA) & Realistic Detector Noise (Priority: P6)

As an AO researcher,  
I want to inject static/quasi-static Non-Common Path Aberrations (NCPA) into the science arm and apply realistic detector noise (Poisson shot noise, Gaussian read noise $\sigma_{\text{RON}}$, dark current),  
So that the synthetic data accurately reflects real-world optical bench and on-sky data limits.

**Why this priority**: NCPA is the primary limiting factor in coronagraphic high-contrast imaging. Simulations without NCPA and detector noise will misrepresent real instruments.

**Independent Test**:
Can be tested by adding a known Zernike NCPA mode (e.g. 50 nm RMS astigmatism) to the science arm with zero atmospheric turbulence. Verifying that the WFS sees zero aberration (common path is flat) while the science focal plane shows degraded PSF with characteristic astigmatic diffraction rings.

**Acceptance Scenarios**:
1. **Given** an NCPA phase map $\phi_{\text{ncpa}}(\vec{r})$ assigned to the science arm, **When** simulated, **Then** the science PSF degrades according to $\phi_{\text{ncpa}}$ while the PyWFS reports a flat wavefront.
2. **Given** detector noise parameters ($\sigma_{\text{RON}} > 0$, dark current, photon noise enabled), **When** simulated, **Then** the output detector counts follow Poisson statistics modulated by additive Gaussian read noise, clipping negative counts or preserving bias offsets as configured.
3. **Given** noise toggles set to `enabled=False`, **When** simulated, **Then** the output data is pure floating-point noise-free intensities.

---

## Edge Cases

- **Zero-Photon Pixels**: In quad-cell slope calculation, illuminated pupil edges or unilluminated subapertures have near-zero total flux $I_{\text{tot}} \approx 0$. The system MUST apply a flux threshold mask to avoid division-by-zero (`NaN` or `Inf`).
- **Nyquist Sampling Violations**: If the user configures a focal plane pixel scale that undersamples the resolution element ($\text{pixels per } \lambda/D < 2$), the system MUST log a clear warning explaining the risk of aliasing.
- **Extreme Modulation Radius**: If $r_{\text{mod}}$ is set very large ($r_{\text{mod}} > 10\,\lambda/D$), modulation points could wrap outside the Fourier grid. The system MUST validate that $r_{\text{mod}} + 1 < N_{\text{grid}} / 4$ and gracefully clamp or warn.
- **Coronagraph Mask Larger than Field of View**: If $r_{\text{FPM}}$ exceeds the focal plane grid radius, the system MUST raise a descriptive configuration error.
- **Broadband Sampling**: When simulating broadband light with $N_{\lambda}$ spectral channels, total execution time scales linearly. The system must provide a default small $N_{\lambda}$ (e.g. 3–5 bins) with progress reporting.

---

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST compute complex 2D Fraunhofer electric field propagation between pupil planes and focal planes using centered 2D Fast Fourier Transforms (FFT) with exact physical grid spacing.
- **FR-002**: System MUST provide a Subaru 8.2m telescope pupil generator supporting central obscuration ($D_{\text{in}}/D_{\text{out}}$), spider struts (width, angles), and clean circular boundary masking.
- **FR-003**: System MUST provide an unocculted clear circular aperture mode as the default toy model.
- **FR-004**: System MUST model the 4-faceted pyramid prism as a 4-quadrant tip-tilt phase mask in the focal plane:
  $$\phi_{\text{pyr}}(u, v) = \alpha \cdot (|u| + |v|)$$
  producing 4 non-overlapping pupil images upon re-imaging.
- **FR-005**: System MUST support circular tip-tilt PyWFS modulation by integrating $N_{\text{mod}}$ sample points along a circle of radius $r_{\text{mod}}$ in units of $\lambda_{\text{wfs}}/D$, with $r_{\text{mod}} = 0$ executing a fast single-pass unmodulated model.
- **FR-006**: System MUST output PyWFS telemetry as either raw 4-pupil intensity images $[4, H_{\text{pup}}, W_{\text{pup}}]$ or normalized $(S_x, S_y)$ slope maps with customizable subaperture valid-pixel thresholding.
- **FR-007**: System MUST provide a classical Lyot coronagraph subsystem consisting of an opaque focal plane mask (FPM radius $r_{\text{FPM}}$ in $\lambda/D$) and an undersized Lyot stop.
- **FR-008**: System MUST support an independent toggle (`enabled: bool`) for every optical component: Atmosphere, Pupil Spiders, NCPA, Coronagraph, PyWFS Modulation, and Detector Noise.
- **FR-009**: System MUST support multi-source scenes: primary star (magnitude/flux), off-axis point companions (contrast, separation, position angle), and internal calibration sources (monochromatic laser, flat field).
- **FR-010**: System MUST support distinct wavelength channels between the WFS path ($\lambda_{\text{wfs}}$) and Science path ($\lambda_{\text{sci}}$), with optional filter bandpass integration over $N_\lambda$ channels.
- **FR-011**: System MUST provide a dynamic Kolmogorov / von Kármán phase screen generator parameterized by Fried parameter $r_0$, outer scale $L_0$, and wind vector $(v_x, v_y)$.
- **FR-012**: System MUST inject independent Non-Common Path Aberrations (NCPA) defined by Zernike coefficients or 2D phase maps into either the science arm or WFS arm.
- **FR-013**: System MUST simulate detector noise including Poisson shot noise, Gaussian read noise ($\sigma_{\text{RON}}$), and dark current on both WFS and Science cameras.
- **FR-014**: System MUST export synchronized data cubes:
  - `wfs_cube`: $[T_{\text{wfs}}, 4, H_{\text{wfs}}, W_{\text{wfs}}]$ or $[T_{\text{wfs}}, 2, H_{\text{pup}}, W_{\text{pup}}]$
  - `psf_cube`: $[T_{\text{psf}}, C_{\text{sci}}, H_{\text{psf}}, W_{\text{psf}}]$
  - `ground_truth_phase`: $[T_{\text{wfs}}, H_{\text{pup}}, W_{\text{pup}}]$
  along with simulation metadata dictionaries.
- **FR-015**: System MUST allow configuring the temporal relationship between WFS frame rate $f_{\text{wfs}}$, science camera exposure time $\tau_{\text{psf}}$, and atmospheric wind motion $dt$.
- **FR-016**: System MUST provide export utilities for NumPy `.npz`, HDF5, and standard astronomical FITS files.

### Key Entities

- **`OpticalSystem`**: The master orchestrator maintaining the configuration, component graph, and optical train execution.
- **`Source`**: Represents an astronomical or calibration light source (spectrum, flux, off-axis position, filter profile).
- **`Aperture`**: Entrance pupil mask (Subaru geometry or circular disk, grid dimensions, pixel scale in meters/pixel).
- **`Atmosphere`**: Dynamic turbulence generator producing evolving 2D phase screens.
- **`PyWFS`**: The pyramid sensor optical model (prism facet angles, modulation parameters, 4-pupil detector layout, slope calculation).
- **`Coronagraph`**: The coronagraph module (FPM mask type/radius, Lyot stop geometry, bypass logic).
- **`Detector`**: Science or WFS camera model (pixel pitch, integration time, quantum efficiency, Poisson & Gaussian noise models).
- **`SimulationDataset`**: The container for synchronized data cubes, telemetry, ground truth phases, and simulation metadata.

---

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001 (Physical Accuracy)**: Direct unocculted PSF for a circular pupil matches the analytical Airy pattern to within $<10^{-3}$ relative intensity error across the first 3 Airy rings.
- **SC-002 (Energy Conservation)**: Total flux integrated over the focal plane equals pupil plane flux to within $>99.9\%$ across all lossless propagation planes.
- **SC-003 (Pyramid Linearity)**: For small tip-tilt aberrations ($\le 0.1\,\text{rad}$), the calculated $(S_x, S_y)$ slopes respond with linear correlation coefficient $R^2 > 0.99$.
- **SC-004 (Execution Performance)**: Generating a single synchronized frame pair (WFS + PSF at $128\times 128$ grid) on a standard MacBook CPU completes in under 50 ms for unmodulated toy mode and under 500 ms for modulated mode ($N_{\text{mod}} = 16$).
- **SC-005 (Synchronization Integrity)**: Generated data cubes preserve exact 1-to-1 temporal alignment between WFS telemetry and science exposure intervals across multi-frame sequences ($N_{\text{frames}} \ge 1000$).
- **SC-006 (Zero-Crash Modular Toggles)**: All 6 optical toggles (Atmosphere, Obscuration, NCPA, Coronagraph, Modulation, Noise) can be toggled in all $2^6 = 64$ combinations without runtime exceptions or output shape mutations.

---

## Assumptions

- Core simulation runs on CPU using NumPy and SciPy; GPU acceleration (CuPy / PyTorch) can be layered on without changing the optical class interfaces.
- The optical propagation treats near-field / pupil-to-focal transformations under paraxial Fraunhofer diffraction, which is standard and physically accurate for astronomical extreme AO systems.
- The default pupil diameter is set to Subaru Telescope ($D = 8.2$ m) with central obscuration ratio $0.29$, but the diameter and obscuration ratio can be scaled arbitrarily.
- The primary data format is in-memory NumPy arrays, with file serializers for FITS, HDF5, and `.npz`.
