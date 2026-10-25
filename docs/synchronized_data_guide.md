# Simulating & Loading Synchronized WFS and PSF Data with Tablasul

This guide demonstrates how to simulate, export, and load synchronized multi-dimensional data cubes between a Pyramid Wavefront Sensor (PyWFS) and a Science Focal Plane Camera (PSF) using `tablasul`.

---

## 1. Synchronized Multi-Scale Optical Data Generation

In high-contrast adaptive optics systems, contemporaneous telemetry across multiple optical arms exhibits distinct spatiotemporal characteristics:

### Key Physical Features Modeled by Tablasul:
1. **Cadence Mismatch (Temporal Oversampling)**: WFS loops typically run at $1\text{--}3\,\text{kHz}$ ($\Delta t \approx 1\text{ ms}$), whereas science cameras integrate over longer exposures ($10\text{--}100\text{ ms}$ or longer). A single science exposure is the temporal integral over $M$ turbulent WFS states.
2. **Chromaticity**: WFS sensing often occurs in the visible (e.g. $800\,\text{nm}$) while science imaging occurs in the near-infrared (e.g. $1.65\,\mu\text{m}$). Optical phase scales inversely with wavelength ($\phi \propto 1/\lambda$).
3. **Non-Common Path Aberrations (NCPA)**: The science path contains quasi-static aberrations that the WFS does not sense.
4. **Data Formats**: Produces synchronized multi-extension FITS and compressed NumPy `.npz` cubes:
   - `wfs_cube`: `[T_wfs, 4, H, W]` raw 4-quadrant intensities or `[T_wfs, 2, H, W]` normalized $(S_x, S_y)$ slopes.
   - `psf_cube`: `[T_sci, H, W]` science camera integrated focal plane PSF.
   - `ground_truth_phase`: `[T_wfs, H, W]` instantaneous entrance pupil wavefront phase.

---

## 2. Generating Synchronized Cubes

```python
from tablasul import (
    OpticalSystem,
    SubaruAperture,
    PyramidWFS,
    LyotCoronagraph,
    PointSource,
    KolmogorovScreen,
    FrozenFlowAtmosphere,
    ZernikeNCPA,
    DetectorModel,
    SynchronizedCubeEngine,
    export_fits,
    export_npz,
)

# 1. Setup optical system
system = OpticalSystem(
    aperture=SubaruAperture(diameter=8.2, grid_size=256, samples_across_pupil=48),
    pywfs=PyramidWFS(modulation_radius=2.0, separation_ratio=2.0),
    coronagraph=LyotCoronagraph(fpm_radius_lod=2.5, enabled=True),
    sources=[PointSource(total_flux=5e6)],
    ncpa=ZernikeNCPA(coefficients={5: 0.15, 6: 0.20}, enabled=True),
    science_detector=DetectorModel(read_noise=2.0, photon_noise=True),
    wfs_detector=DetectorModel(read_noise=1.0, photon_noise=True),
    wfs_wavelength=800e-9,
    science_wavelength=1650e-9,
)

# 2. Dynamic frozen-flow turbulence
atmos = FrozenFlowAtmosphere(
    screen_generator=KolmogorovScreen(r0=0.18, L0=25.0),
    wind_velocity=(12.0, 3.0),
)

# 3. Synchronized engine with 10:1 oversampling
engine = SynchronizedCubeEngine(
    system=system,
    atmosphere=atmos,
    wfs_frame_rate=1000.0,       # 1 kHz WFS
    science_exposure_time=0.010, # 10 ms science exposure
    wfs_output_mode="slopes",    # or "quadrants"
)

# Generate 500 science frames (accompanied by 5000 WFS frames and phase screens)
dataset = engine.generate_cubes(n_science_frames=500)

# Save to disk
export_fits(dataset, "data/synchronized_cubes.fits")
export_npz(dataset, "data/synchronized_cubes.npz")
```

---

## 3. Loading Synchronized Data in PyTorch

```python
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader

class SynchronizedOpticalDataset(Dataset):
    """PyTorch Dataset for paired WFS telemetry and integrated science PSF frames."""

    def __init__(self, npz_path: str):
        data = np.load(npz_path, allow_pickle=True)
        self.wfs = torch.from_numpy(data["wfs_cube"]).float()        # (T_wfs, 2, H_pup, W_pup)
        self.psf = torch.from_numpy(data["psf_cube"]).float()        # (T_sci, H_focal, W_focal)
        self.phase = torch.from_numpy(data["ground_truth_phase"]).float()

        self.substeps = self.wfs.shape[0] // self.psf.shape[0]

    def __len__(self):
        return self.psf.shape[0]

    def __getitem__(self, idx):
        # M WFS frames corresponding to this 1 science exposure
        start = idx * self.substeps
        end = start + self.substeps
        wfs_chunk = self.wfs[start:end] # (M, 2, H, W)
        science_psf = self.psf[idx]     # (H, W)
        return wfs_chunk, science_psf
```
