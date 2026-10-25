# Tasks: First-Principles PyWFS, Coronagraph, and Science PSF Data Simulator

**Branch**: `001-core-simulator` | **Spec**: [spec.md](file:///Users/seamus/Research/tablasul/specs/001-core-simulator/spec.md) | **Plan**: [plan.md](file:///Users/seamus/Research/tablasul/specs/001-core-simulator/plan.md)

---

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Project scaffolding, package configuration, and dependency setup

- [x] T001 Create Python package directory structure in `src/tablasul/` and `tests/`
- [x] T002 Configure `pyproject.toml` with dependencies (`numpy`, `scipy`, `astropy`, `pytest`)
- [x] T003 Configure test runner in `tests/conftest.py` with common fixtures (grid sizes, wavelength standards)

---

## Phase 2: Foundational (Wave Optics Core)

**Purpose**: First-principles complex electric field and Fraunhofer propagation engine

- [x] T004 Implement `ComplexElectricField` in `src/tablasul/core/field.py` (amplitude, phase, wavelength, grid spacing, energy conservation checks)
- [x] T005 Implement paraxial Fraunhofer 2D Fourier propagation between pupil and focal planes in `src/tablasul/core/field.py`
- [x] T006 Write unit tests in `tests/test_field_propagation.py` verifying Parseval's theorem ($\sum |E_f|^2 = \sum |E_0|^2$) and exact grid scaling

---

## Phase 3: User Story 1 - Minimal Diffraction-Limited Benchmark ("Toy Model") (Priority: P1) 🎯 MVP

**Goal**: Establish an end-to-end unmodulated PyWFS and direct Airy PSF pipeline on a clear circular aperture with zero noise

- [x] T007 [US1] Implement `CircularAperture` in `src/tablasul/core/aperture.py` (diameter, grid resolution, transmission mask with $n // 2$ centering)
- [x] T008 [US1] Implement static unmodulated `PyramidWFS` phase mask in `src/tablasul/optics/pywfs.py` (4-quadrant tip-tilt prism splitting focal field into 4 pupil quadrants)
- [x] T009 [US1] Implement direct focal plane PSF calculation in `src/tablasul/engine/simulator.py`
- [x] T010 [US1] Write test in `tests/test_field_propagation.py` verifying that unaberrated direct PSF matches the analytical J1 Bessel function Airy disk to within $<10^{-3}$ relative intensity error
- [x] T011 [US1] Write test in `tests/test_pywfs.py` verifying that unmodulated PyWFS under zero aberration yields 4 identical pupil quadrants ($I_1 = I_2 = I_3 = I_4$)

---

## Phase 4: User Story 2 - Composable Optical Train & Lyot Coronagraph (Priority: P2)

**Goal**: Implement Subaru 8.2m pupil geometry (spiders, central obscuration) and a classical Lyot coronagraph with bypass toggles

- [x] T012 [US2] Implement `SubaruAperture` in `src/tablasul/core/aperture.py` (central obscuration ratio $0.29$, 4 spider struts, toggleable geometry)
- [x] T013 [US2] Implement `LyotCoronagraph` in `src/tablasul/optics/coronagraph.py` (focal plane mask radius $r_{\text{FPM}}$, undersized Lyot stop, bypass logic when `enabled=False`)
- [x] T014 [US2] Integrate coronagraph and pupil into `src/tablasul/engine/simulator.py`
- [x] T015 [US2] Write unit tests in `tests/test_coronagraph.py` asserting on-axis starlight rejection $>99\%$, off-axis companion transmission $>80\%$, and identical direct PSF output when `enabled=False`

---

## Phase 5: User Story 3 - PyWFS Tip-Tilt Modulation & Slope Computation (Priority: P3)

**Goal**: Implement circular tip-tilt dynamic modulation and dual output telemetry (raw 4-pupil images and normalized slopes)

- [x] T016 [US3] Implement dynamic circular tip-tilt modulation in `src/tablasul/optics/pywfs.py` (integrating $N_{\text{mod}}$ points on circle of radius $r_{\text{mod}}$ in $\lambda_{\text{wfs}}/D$, with $r_{\text{mod}} = 0$ fast path)
- [x] T017 [US3] Implement quad-cell slope calculation $(S_x, S_y)$ in `src/tablasul/optics/pywfs.py` with subaperture threshold masking to prevent division-by-zero
- [x] T018 [US3] Write unit tests in `tests/test_pywfs.py` verifying linear slope response under small tip-tilt aberrations

---

## Phase 6: User Story 4 - Multi-Source & Spectral Filter Modeling (Priority: P4)

**Goal**: Support diverse astronomical sources, off-axis companions, internal calibration sources, and distinct WFS/Science spectral filter bands

- [x] T019 [US4] Implement `Source` hierarchy in `src/tablasul/core/source.py` (`PointSource`, `CompanionSource` with contrast/separation/position angle, and `CalibrationSource` for monochromatic lasers and flat lamps)
- [x] T020 [US4] Implement `FilterBand` in `src/tablasul/core/source.py` (top-hat, Gaussian, standard photometric bands, discrete channel discretization)
- [x] T021 [US4] Integrate multi-wavelength propagation and distinct WFS/Science arm bandpasses into `src/tablasul/engine/simulator.py`
- [x] T022 [US4] Write unit tests in `tests/test_source.py` verifying chromatic phase scaling ($\phi \propto 1/\lambda$), off-axis companion positioning, and calibration source clean single-wavelength execution

---

## Phase 7: User Story 5 - Synchronized Spatiotemporal Data Cubes (Priority: P5)

**Goal**: Implement dynamic atmospheric phase evolution and generate synchronized multi-dimensional cubes with temporal cadence controls

- [x] T023 [US5] Implement `KolmogorovScreen` and `FrozenFlowAtmosphere` in `src/tablasul/core/atmosphere.py` (Fried parameter $r_0$, outer scale $L_0$, wind vector $\vec{v}$, spatial interpolation drift)
- [x] T024 [US5] Implement `SynchronizedCubeEngine` in `src/tablasul/engine/sync.py` handling oversampled sub-step integration between WFS frame rate $f_{\text{wfs}}$ and science exposure time $\tau_{\text{psf}}$
- [x] T025 [US5] Implement cube assembly in `src/tablasul/engine/simulator.py` yielding `wfs_cube`, `psf_cube`, and `ground_truth_phase`
- [x] T026 [US5] Write unit tests in `tests/test_synchronization.py` verifying cube dimension contracts, 1-to-1 temporal alignment, and long-exposure PSF averaging

---

## Phase 8: User Story 6 - NCPA, Detector Noise & File Serialization (Priority: P6)

**Goal**: Add Non-Common Path Aberrations, Poisson/Gaussian detector noise models, and standard astronomical FITS / HDF5 / NPZ export

- [x] T027 [US6] Implement `ZernikeNCPA` and static phase injection in `src/tablasul/optics/ncpa.py` (applying differential aberrations exclusively to science or WFS path)
- [x] T028 [US6] Implement `DetectorModel` in `src/tablasul/detector/camera.py` (Poisson shot noise, Gaussian read noise $\sigma_{\text{RON}}$, dark current, bypass when `enabled=False`)
- [x] T029 [US6] Implement data serializers in `src/tablasul/io/export.py` (FITS with complete optical headers, HDF5, and `.npz`)
- [x] T030 [US6] Write unit tests in `tests/test_ncpa_and_noise.py` verifying NCPA presence in science PSF with flat WFS telemetry, and noise statistical distributions
- [x] T031 [US6] Write comprehensive toggle matrix test in `tests/test_toggles.py` executing all 64 toggle permutations asserting non-empty arrays, shape stability, and zero runtime crashes
