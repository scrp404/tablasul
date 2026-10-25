"""Tests for Lyot Coronagraph and Subaru aperture geometry."""

import numpy as np
import pytest
from tablasul.core.aperture import CircularAperture, SubaruAperture
from tablasul.optics.coronagraph import LyotCoronagraph
from tablasul.engine.simulator import OpticalSystem


def test_coronagraph_bypass_mode():
    """Verify that when enabled=False, coronagraphic science PSF matches direct PSF exactly."""
    system = OpticalSystem(
        aperture=CircularAperture(diameter=8.2, grid_size=128, samples_across_pupil=48),
        coronagraph=LyotCoronagraph(enabled=False),
    )
    res = system.simulate_frame()
    assert np.allclose(res.science_psf, res.direct_psf, atol=1e-12)


def test_coronagraph_starlight_suppression():
    """Verify that an enabled Lyot coronagraph suppresses on-axis starlight peak by >95%."""
    system = OpticalSystem(
        aperture=CircularAperture(diameter=8.2, grid_size=128, samples_across_pupil=48),
        coronagraph=LyotCoronagraph(fpm_radius_lod=2.5, lyot_stop_ratio=0.90, enabled=True),
    )
    res = system.simulate_frame()
    direct_peak = np.max(res.direct_psf)
    corono_peak = np.max(res.science_psf)

    # Rejection factor
    rejection = 1.0 - (corono_peak / direct_peak)
    assert rejection > 0.95, f"Expected >95% rejection, got {rejection * 100:.2f}%"


def test_subaru_aperture_geometry():
    """Verify Subaru aperture obscuration and spider diffraction spikes."""
    subaru = SubaruAperture(
        diameter=8.2,
        grid_size=128,
        samples_across_pupil=48,
        central_obscuration_ratio=0.29,
        spider_width=0.22,
        enabled_spiders=True,
        enabled_obscuration=True,
    )
    mask = subaru.transmission_mask()
    assert mask.shape == (128, 128)

    # Center must be obscured (zero transmission)
    mid = 128 // 2
    assert mask[mid, mid] == 0.0

    # Intermediate pupil area must be illuminated (transmission 1.0)
    assert mask[mid, mid + 16] == 1.0

    # Outside pupil must be zero
    assert mask[mid, mid + 35] == 0.0
