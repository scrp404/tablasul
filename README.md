# Tablasul

`tablasul` generates physically and mathematically grounded synchronized data cubes between a **Pyramid Wavefront Sensor (PyWFS)** and a **Science Focal Plane Camera (PSF)** behind an optional **Lyot Coronagraph**.

---

## Key Features

1. **First-Principles Wave Optics**: Paraxial Fraunhofer propagation via centered 2D FFTs with exact Parseval energy conservation and physical units (meters, radians, wavelengths $\lambda$, flux in photons).
2. **"Toy-to-Full-Scale" Toggle Philosophy**: Every single physical subsystem is independently toggleable (`enabled: bool`). With toggles off, it runs as an ultra-fast diffraction-limited Airy toy model. With toggles on, it models the Subaru 8.2m pupil, Kolmogorov turbulence, NCPA, modulated PyWFS, Lyot coronagraph, and detector noise.
3. **Synchronized Multidimensional Cubes**:
   - `wfs_cube`: `[T_wfs, 4, H, W]` raw 4-pupil intensity images, or `[T_wfs, 2, H, W]` normalized $(S_x, S_y)$ quad-cell slopes.
   - `psf_cube`: `[T_sci, H, W]` science camera PSF integrated over exposure time.
   - `ground_truth_phase`: `[T_wfs, H, W]` instantaneous entrance pupil wavefront phase screens.
   - Handles temporal oversampling where $M$ fast WFS frames integrate into a single long-exposure science image.
4. **Multi-Source & Off-Axis Companions**: Supports on-axis primary stars, off-axis exoplanet companions (contrast, separation in $\lambda/D$, position angle), and internal calibration sources (monochromatic lasers, flat lamps).
5. **NCPA & Detector Noise**: Non-Common Path Aberrations (Zernike polynomials) applied exclusively to the science path, plus Poisson photon shot noise, Gaussian read noise ($\sigma_{\text{RON}}$), and dark current.
6. **Data Export**: Direct serializers for multi-extension astronomical FITS and compressed NumPy `.npz` files.

---

## Documentation

* **[Optical Theory & Mathematical Principles](docs/optical_theory.md)**: Mathematical derivations of Fraunhofer propagation, 4-faceted pyramid phase mask, circular tip-tilt modulation, and Lyot coronagraph filtering.

---

## Quickstart

### 1. Minimal Frame Simulation ("Toy Model")

```python
from tablasul import OpticalSystem, CircularAperture, PyramidWFS

# Fast diffraction-limited toy model
system = OpticalSystem(
    aperture=CircularAperture(diameter=8.2, grid_size=256, samples_across_pupil=48),
    pywfs=PyramidWFS(modulation_radius=0.0, separation_ratio=2.0),
    wfs_wavelength=800e-9,                   # 800 nm visible
    science_wavelength=1650e-9,              # 1650 nm NIR H-band
)

result = system.simulate_frame()
print("Science PSF shape:", result.science_psf.shape)
print("PyWFS detector shape:", result.pywfs_detector.shape)
print("Slopes (Sx, Sy):", result.slopes[0].shape, result.slopes[1].shape)
```

---

### 2. High-Fidelity Simulation with Off-Axis Planet & Coronagraph

```python
from tablasul import (
    OpticalSystem,
    SubaruAperture,
    PyramidWFS,
    LyotCoronagraph,
    PointSource,
    CompanionSource,
    ZernikeNCPA,
    DetectorModel,
)

# On-axis star + 2x10^-3 off-axis companion at 6.0 lambda/D
star = PointSource(total_flux=1e7)
planet = CompanionSource(contrast=2e-3, separation_lod=6.0, position_angle_deg=45.0)

system = OpticalSystem(
    aperture=SubaruAperture(diameter=8.2, central_obscuration_ratio=0.29, spider_width=0.22),
    pywfs=PyramidWFS(modulation_radius=2.0, separation_ratio=2.2, n_modulation_points=16),
    coronagraph=LyotCoronagraph(fpm_radius_lod=2.5, lyot_stop_ratio=0.88, edge_taper_lod=0.5, enabled=True),
    sources=[star, planet],
    ncpa=ZernikeNCPA(coefficients={5: 0.10, 6: 0.15}), # Astigmatism in science arm
    science_detector=DetectorModel(read_noise=2.0, photon_noise=True),
    wfs_detector=DetectorModel(read_noise=1.0, photon_noise=True),
)

result = system.simulate_frame()
```

---

### 3. Generating Synchronized Cubes

```python
from tablasul import (
    OpticalSystem,
    SynchronizedCubeEngine,
    KolmogorovScreen,
    FrozenFlowAtmosphere,
    export_fits,
    export_npz,
)

atmosphere = FrozenFlowAtmosphere(
    screen_generator=KolmogorovScreen(r0=0.18, L0=25.0),
    wind_velocity=(10.0, 5.0),
)

system = OpticalSystem()

# 1000 Hz WFS, 10 ms science exposures (10:1 temporal oversampling ratio)
engine = SynchronizedCubeEngine(
    system=system,
    atmosphere=atmosphere,
    wfs_frame_rate=1000.0,
    science_exposure_time=0.010,
    wfs_output_mode="slopes", # or "quadrants"
)

dataset = engine.generate_cubes(n_science_frames=50)

print("WFS cube shape:", dataset.wfs_cube.shape)           # (500, 2, 48, 48)
print("PSF cube shape:", dataset.psf_cube.shape)           # (50, 256, 256)
print("Phase cube shape:", dataset.ground_truth_phase.shape) # (500, 256, 256)

export_fits(dataset, "output/cubes.fits")
export_npz(dataset, "output/cubes.npz")
```
