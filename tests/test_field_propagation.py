"""Tests for complex wave field propagation and analytical Airy pattern verification."""

import numpy as np
import pytest
from scipy.special import j1
from tablasul.core.field import ComplexElectricField
from tablasul.core.aperture import CircularAperture
from tablasul.engine.simulator import OpticalSystem


def test_parseval_energy_conservation():
    """Verify that Fraunhofer FFT propagation strictly conserves total energy (Parseval)."""
    n = 128
    data = np.zeros((n, n), dtype=np.complex128)
    mid = n // 2
    y, x = np.ogrid[:n, :n]
    r = np.hypot(x - mid, y - mid)
    data[r <= 32] = 1.0

    pupil_field = ComplexElectricField(
        data=data,
        wavelength=800e-9,
        pixel_scale=8.2 / 64,
        plane="pupil",
    )

    pupil_energy = pupil_field.total_energy
    focal_field = pupil_field.propagate_to_focal()
    focal_energy = focal_field.total_energy

    assert np.isclose(focal_energy, pupil_energy, rtol=1e-10)

    # Re-propagate back to pupil
    reimaged_pupil = focal_field.propagate_to_pupil()
    reimaged_energy = reimaged_pupil.total_energy

    assert np.isclose(reimaged_energy, pupil_energy, rtol=1e-10)
    assert np.allclose(reimaged_pupil.data, pupil_field.data, atol=1e-10)


def test_analytical_airy_disk_verification():
    """Verify that direct unobstructed circular aperture PSF matches analytical J1 Airy disk."""
    n = 256
    samples_pupil = 64
    aperture = CircularAperture(diameter=8.2, grid_size=n, samples_across_pupil=samples_pupil)
    field = aperture.create_field(wavelength=800e-9, total_flux=1.0)
    focal = field.propagate_to_focal()
    psf = focal.intensity

    # Normalize PSF peak to 1
    psf_norm = psf / np.max(psf)

    # Coordinate system in units of lambda / D
    pixel_scale_lod = samples_pupil / n

    center = n // 2
    y, x = np.ogrid[:n, :n]
    r_lod = np.hypot(x - center, y - center) * pixel_scale_lod

    # Analytical Airy intensity: [2 * J1(pi * r) / (pi * r)]^2
    with np.errstate(divide="ignore", invalid="ignore"):
        arg = np.pi * r_lod
        analytical_airy = np.where(r_lod == 0, 1.0, (2.0 * j1(arg) / arg) ** 2)

    # Compare within central core and first 2 rings (r <= 2.5 lambda/D)
    mask = r_lod <= 2.5
    max_error = np.max(np.abs(psf_norm[mask] - analytical_airy[mask]))
    assert max_error < 0.005, f"Airy profile discrepancy too large: {max_error}"


def test_mft_toggleable_propagation():
    """Verify that Matrix Fourier Transform (MFT) can be toggled on demand and matches FFT peak."""
    system_fft = OpticalSystem(
        aperture=CircularAperture(diameter=8.2, grid_size=128, samples_across_pupil=32),
        science_propagator="fft",
    )
    system_mft = OpticalSystem(
        aperture=CircularAperture(diameter=8.2, grid_size=128, samples_across_pupil=32),
        science_propagator="mft",
        science_fov_lod=32.0,  # Match FFT FOV
        science_pixel_dims=128,
    )
    res_fft = system_fft.simulate_frame()
    res_mft = system_mft.simulate_frame()

    # Peak intensities must match to within 1%
    fft_peak = np.max(res_fft.direct_psf)
    mft_peak = np.max(res_mft.direct_psf)
    assert np.isclose(fft_peak, mft_peak, rtol=0.01)
    assert res_mft.science_psf.shape == (128, 128)
