"""Tests for NCPA, detector noise models, and data export."""

import tempfile
from pathlib import Path
import numpy as np
import pytest
from astropy.io import fits
from tablasul.core.aperture import CircularAperture
from tablasul.optics.ncpa import ZernikeNCPA
from tablasul.detector.camera import DetectorModel
from tablasul.engine.simulator import OpticalSystem
from tablasul.engine.sync import SynchronizedCubeEngine
from tablasul.io.export import export_fits, export_npz


def test_ncpa_isolation_between_arms():
    """Verify that NCPA degrades science PSF while WFS sees a flat wavefront."""
    # Add vertical astigmatism (Z6) to science arm
    ncpa = ZernikeNCPA(coefficients={6: 0.8}, enabled=True)
    system = OpticalSystem(
        aperture=CircularAperture(diameter=8.2, grid_size=64, samples_across_pupil=32),
        ncpa=ncpa,
    )
    res = system.simulate_frame()

    # WFS should be completely invariant to NCPA on the science arm
    sys_clean = OpticalSystem(
        aperture=CircularAperture(diameter=8.2, grid_size=64, samples_across_pupil=32),
        ncpa=ZernikeNCPA(coefficients={6: 0.8}, enabled=False),
    )
    res_clean = sys_clean.simulate_frame()
    assert np.allclose(res.slopes[0], res_clean.slopes[0], atol=1e-12)
    assert np.allclose(res.slopes[1], res_clean.slopes[1], atol=1e-12)

    # Science PSF peak should be degraded compared to unaberrated PSF
    strehl_ratio = np.max(res.science_psf) / np.max(res_clean.science_psf)
    assert strehl_ratio < 0.85, f"Expected degraded Strehl, got {strehl_ratio}"


def test_detector_noise_statistics():
    """Verify Poisson shot noise and Gaussian read noise statistics."""
    # Test flat uniform input
    det = DetectorModel(read_noise=2.0, dark_current=0.0, photon_noise=True, gain=1.0, bias=100.0)
    flat_flux = np.full((100, 100), 50.0, dtype=np.float64)
    noisy = det.apply(flat_flux, exposure_time=1.0, seed=42)

    # Expected mean = 50 + 100 (bias) = 150
    # Expected variance = Var(Poisson(50)) + Var(Gaussian(2)) = 50 + 4 = 54
    assert np.isclose(np.mean(noisy), 150.0, atol=2.0)
    assert np.isclose(np.var(noisy), 54.0, atol=6.0)


def test_fits_and_npz_export():
    """Verify FITS and NPZ export of synchronized cubes."""
    system = OpticalSystem(
        aperture=CircularAperture(diameter=8.2, grid_size=64, samples_across_pupil=32),
    )
    engine = SynchronizedCubeEngine(system=system, wfs_frame_rate=100.0, science_exposure_time=0.02)
    dataset = engine.generate_cubes(n_science_frames=2)

    with tempfile.TemporaryDirectory() as tmpdir:
        fits_path = Path(tmpdir) / "test_cubes.fits"
        npz_path = Path(tmpdir) / "test_cubes.npz"

        export_fits(dataset, fits_path)
        assert fits_path.exists()
        with fits.open(fits_path) as hdul:
            assert "WFS_CUBE" in hdul
            assert "PSF_CUBE" in hdul
            assert "PHASE_CUBE" in hdul
            assert hdul["WFS_CUBE"].data.shape == dataset.wfs_cube.shape

        export_npz(dataset, npz_path)
        assert npz_path.exists()
        with np.load(npz_path, allow_pickle=True) as loaded:
            assert "wfs_cube" in loaded
            assert "psf_cube" in loaded
            assert "ground_truth_phase" in loaded
