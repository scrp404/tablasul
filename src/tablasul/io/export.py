"""Data cube export and serialization utilities (FITS and NPZ)."""

from __future__ import annotations
from pathlib import Path
from typing import TYPE_CHECKING
import numpy as np
from astropy.io import fits

if TYPE_CHECKING:
    from ..engine.sync import SynchronizedDataset


def export_fits(dataset: SynchronizedDataset, filepath: str | Path) -> Path:
    """Export synchronized dataset to a multi-extension astronomical FITS file."""
    path = Path(filepath)
    path.parent.mkdir(parents=True, exist_ok=True)

    # Primary HDU with metadata
    pri_hdr = fits.Header()
    pri_hdr["COMMENT"] = "Tablasul Synchronized WFS and PSF Simulation Dataset"
    for k, v in dataset.metadata.items():
        # FITS header keys must be <= 8 chars
        short_k = k[:8].upper()
        if isinstance(v, (int, float, str, bool)):
            pri_hdr[short_k] = v

    hdul = fits.HDUList([fits.PrimaryHDU(header=pri_hdr)])

    # HDU 1: WFS Telemetry Cube
    wfs_hdu = fits.ImageHDU(data=dataset.wfs_cube.astype(np.float32), name="WFS_CUBE")
    wfs_hdu.header["DESC"] = "Synchronized WFS telemetry"
    hdul.append(wfs_hdu)

    # HDU 2: Science PSF Cube
    psf_hdu = fits.ImageHDU(data=dataset.psf_cube.astype(np.float32), name="PSF_CUBE")
    psf_hdu.header["DESC"] = "Synchronized Science focal plane PSF"
    hdul.append(psf_hdu)

    # HDU 3: Ground Truth Wavefront Phase Cube
    phase_hdu = fits.ImageHDU(data=dataset.ground_truth_phase.astype(np.float32), name="PHASE_CUBE")
    phase_hdu.header["DESC"] = "Ground truth wavefront phase in radians"
    hdul.append(phase_hdu)

    hdul.writeto(path, overwrite=True)
    return path


def export_npz(dataset: SynchronizedDataset, filepath: str | Path) -> Path:
    """Export synchronized dataset to a compressed NumPy .npz archive."""
    path = Path(filepath)
    path.parent.mkdir(parents=True, exist_ok=True)

    np.savez_compressed(
        path,
        wfs_cube=dataset.wfs_cube,
        psf_cube=dataset.psf_cube,
        ground_truth_phase=dataset.ground_truth_phase,
        timestamps_wfs=dataset.timestamps_wfs,
        timestamps_psf=dataset.timestamps_psf,
        metadata=dataset.metadata,
    )
    return path
