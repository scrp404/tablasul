"""Tests for multi-source scenes, off-axis companions, and filter bands."""

import numpy as np
import pytest
from tablasul.core.aperture import CircularAperture
from tablasul.core.source import (
    PointSource,
    CompanionSource,
    CalibrationSource,
    FilterBand,
)
from tablasul.optics.coronagraph import LyotCoronagraph
from tablasul.engine.simulator import OpticalSystem


def test_off_axis_companion_detection():
    """Verify that an off-axis companion creates a localized PSF peak at expected separation."""
    n = 128
    samples_pupil = 32
    star = PointSource(total_flux=1.0)
    # Companion at 6 lambda/D along +x axis (position_angle=0)
    sep_lod = 6.0
    contrast = 1e-2  # 1% brightness
    planet = CompanionSource(contrast=contrast, separation_lod=sep_lod, position_angle_deg=0.0)

    # Use Lyot coronagraph with FPM radius 3.0 lambda/D to suppress on-axis star
    system = OpticalSystem(
        aperture=CircularAperture(diameter=8.2, grid_size=n, samples_across_pupil=samples_pupil),
        coronagraph=LyotCoronagraph(fpm_radius_lod=3.0, lyot_stop_ratio=0.90, enabled=True),
        sources=[star, planet],
    )
    res = system.simulate_frame()
    psf = res.science_psf

    # The companion should dominate outside the FPM
    mid = n // 2
    pixel_scale_lod = samples_pupil / n  # lod per pixel
    expected_x_pixel = mid + int(round(sep_lod / pixel_scale_lod))

    # Look for peak along the row through the center
    row = psf[mid, :]
    peak_col = np.argmax(row)
    assert np.isclose(peak_col, expected_x_pixel, atol=2)


def test_filter_band_sampling():
    """Verify FilterBand sampling channels and weights."""
    band_h = FilterBand.standard("H")
    assert band_h.central_wavelength == 1650e-9
    channels = band_h.sample_channels(n_channels=5)
    assert len(channels) == 5
    # Sum of weights must equal 1.0
    total_weight = sum(w for _, w in channels)
    assert np.isclose(total_weight, 1.0)
    # Middle channel should match central wavelength
    assert np.isclose(channels[2][0], 1650e-9)


def test_calibration_laser_source():
    """Verify CalibrationSource produces zero tilt."""
    laser = CalibrationSource(wavelength=635e-9, total_flux=1.0)
    tilt = laser.get_phase_tilt(grid_size=64, samples_across_pupil=32, wavelength=635e-9)
    assert np.all(tilt == 0.0)
