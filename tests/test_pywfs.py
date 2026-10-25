"""Tests for Pyramid WFS optical model and slope response."""

import numpy as np
import pytest
from tablasul.core.aperture import CircularAperture
from tablasul.optics.pywfs import PyramidWFS
from tablasul.engine.simulator import OpticalSystem


def test_unmodulated_pywfs_quadrant_symmetry():
    """Verify that unmodulated PyWFS with zero aberration produces balanced quadrants."""
    system = OpticalSystem(
        aperture=CircularAperture(diameter=8.2, grid_size=128, samples_across_pupil=48),
        pywfs=PyramidWFS(modulation_radius=0.0, separation_ratio=1.4),
        wavelength=800e-9,
    )
    res = system.simulate_frame()
    I1, I2, I3, I4 = res.pywfs_quadrants

    f1, f2, f3, f4 = np.sum(I1), np.sum(I2), np.sum(I3), np.sum(I4)
    fluxes = np.array([f1, f2, f3, f4])
    mean_flux = np.mean(fluxes)

    # All 4 quadrants must be balanced to within 2%
    relative_diff = np.abs(fluxes - mean_flux) / mean_flux
    assert np.all(relative_diff < 0.02), f"Quadrant fluxes unbalanced: {fluxes}"

    Sx, Sy, mask = res.slopes
    # Mean slopes across illuminated pupil under flat wavefront must be near zero
    assert np.abs(np.mean(Sx[mask])) < 0.01
    assert np.abs(np.mean(Sy[mask])) < 0.01


def test_pywfs_tip_tilt_slope_response():
    """Verify that small tip-tilt aberrations produce expected differential slope signals."""
    system = OpticalSystem(
        aperture=CircularAperture(diameter=8.2, grid_size=128, samples_across_pupil=32),
        pywfs=PyramidWFS(modulation_radius=1.5, n_modulation_points=16),
        wavelength=800e-9,
    )

    n = 128
    mid = n // 2
    y, x = np.ogrid[:n, :n]
    # Apply pure horizontal tilt (phase ramp along x)
    tilt_amp = 1.0  # radians of tilt across pupil
    phase_x = tilt_amp * (x - mid) / 16.0

    res_tilt = system.simulate_frame(phase_aberration=phase_x)
    Sx, Sy, mask = res_tilt.slopes

    mean_sx = np.mean(Sx[mask])
    mean_sy = np.mean(Sy[mask])

    # Horizontal slope should be significantly positive, vertical slope near zero
    assert mean_sx > 0.05, f"Expected positive Sx slope, got {mean_sx}"
    assert np.abs(mean_sy) < 0.015, f"Expected near-zero Sy slope, got {mean_sy}"
