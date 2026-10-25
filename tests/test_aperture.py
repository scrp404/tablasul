"""Tests for aperture classes including analytical Subaru aperture and PupilImageAperture."""

import tempfile
from pathlib import Path
import numpy as np
import pytest
from astropy.io import fits
from tablasul.core.aperture import CircularAperture, SubaruAperture, PupilImageAperture
from tablasul.optics.pywfs import PyramidWFS
from tablasul.optics.coronagraph import LyotCoronagraph
from tablasul.engine.simulator import OpticalSystem


def test_subaru_aperture_structure():
    """Verify Subaru aperture geometry (primary diameter, central obscuration, 4 spiders)."""
    grid_size = 256
    ap = SubaruAperture(
        diameter=8.2,
        grid_size=grid_size,
        samples_across_pupil=64,
        central_obscuration_ratio=0.29,
        spider_width=0.22,
    )
    mask = ap.transmission_mask()
    assert mask.shape == (grid_size, grid_size)
    assert np.max(mask) == 1.0
    assert np.min(mask) == 0.0

    # Center must be obscured
    mid = grid_size // 2
    assert mask[mid, mid] == 0.0

    # Outer corner must be zero
    assert mask[0, 0] == 0.0


def test_pupil_image_aperture_from_array():
    """Verify loading PupilImageAperture from a 2D numpy array."""
    grid_size = 128
    dummy_pupil = np.zeros((grid_size, grid_size), dtype=np.float64)
    y, x = np.ogrid[:grid_size, :grid_size]
    r = np.hypot(x - grid_size // 2, y - grid_size // 2)
    dummy_pupil[(r <= 32) & (r > 10)] = 255.0  # Annular test pupil with raw pixel values

    ap = PupilImageAperture(
        image_or_path=dummy_pupil,
        diameter=8.2,
    )
    mask = ap.transmission_mask()
    assert mask.shape == (grid_size, grid_size)
    assert np.isclose(np.max(mask), 1.0)
    assert np.isclose(mask[grid_size // 2, grid_size // 2], 0.0)


def test_pupil_image_aperture_fits_and_npy(tmp_path: Path):
    """Verify loading from FITS and NPY files with resampling."""
    orig_size = 64
    y, x = np.ogrid[:orig_size, :orig_size]
    r = np.hypot(x - orig_size // 2, y - orig_size // 2)
    disk = (r <= 20).astype(np.float32)

    # 1. Test FITS
    fits_file = tmp_path / "subaru_pupil.fits"
    fits.writeto(fits_file, disk, overwrite=True)

    target_grid = 128
    ap_fits = PupilImageAperture(
        image_or_path=fits_file,
        diameter=8.2,
        grid_size=target_grid,
    )
    mask_fits = ap_fits.transmission_mask()
    assert mask_fits.shape == (target_grid, target_grid)
    assert np.max(mask_fits) > 0.9

    # 2. Test NPY
    npy_file = tmp_path / "subaru_pupil.npy"
    np.save(npy_file, disk)
    ap_npy = PupilImageAperture(image_or_path=npy_file, diameter=8.2)
    assert ap_npy.transmission_mask().shape == (orig_size, orig_size)


def test_pupil_image_aperture_in_optical_system(tmp_path: Path):
    """Verify that PupilImageAperture seamlessly works inside OpticalSystem end-to-end."""
    grid_size = 128
    y, x = np.ogrid[:grid_size, :grid_size]
    r = np.hypot(x - grid_size // 2, y - grid_size // 2)
    pupil_data = ((r <= 24) & (r >= 7)).astype(np.float64)

    fits_path = tmp_path / "custom_pupil.fits"
    fits.writeto(fits_path, pupil_data, overwrite=True)

    aperture = PupilImageAperture(
        image_or_path=fits_path,
        diameter=8.2,
        grid_size=grid_size,
    )

    optics = OpticalSystem(
        aperture=aperture,
        pywfs=PyramidWFS(),
        coronagraph=LyotCoronagraph(),
    )

    result = optics.simulate_frame()
    assert result.pywfs_detector is not None
    assert result.science_psf is not None
    assert result.pywfs_detector.shape == (grid_size, grid_size)
    assert result.science_psf.shape == (grid_size, grid_size)
    assert np.sum(result.pywfs_detector) > 0
    assert np.sum(result.science_psf) > 0


def test_pupil_image_aperture_invalid_inputs():
    """Verify clear error handling for non-existent files or invalid dimensions."""
    with pytest.raises(FileNotFoundError):
        PupilImageAperture("non_existent_file_path_12345.fits")

    with pytest.raises(ValueError):
        PupilImageAperture(np.ones((10, 10, 10)))  # 3D array
