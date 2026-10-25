"""Comprehensive toggle matrix test (all 6 optical subsystems)."""

import itertools
import numpy as np
import pytest
from tablasul.core.aperture import CircularAperture, SubaruAperture
from tablasul.core.atmosphere import KolmogorovScreen
from tablasul.optics.pywfs import PyramidWFS
from tablasul.optics.coronagraph import LyotCoronagraph
from tablasul.optics.ncpa import ZernikeNCPA
from tablasul.detector.camera import DetectorModel
from tablasul.engine.simulator import OpticalSystem


def test_all_toggle_permutations():
    """Verify that all combinations of the 6 optical toggles execute without crashes or shape mutations."""
    # 6 independent binary toggles:
    # 1. Obscuration (Subaru vs Clear)
    # 2. Coronagraph (enabled vs disabled)
    # 3. Modulation (modulated vs unmodulated)
    # 4. NCPA (enabled vs disabled)
    # 5. Science Detector Noise (enabled vs disabled)
    # 6. WFS Detector Noise (enabled vs disabled)

    toggles = [False, True]
    grid_size = 64
    samples = 32

    # Test all 64 permutations
    for obs_on, corono_on, mod_on, ncpa_on, sci_noise_on, wfs_noise_on in itertools.product(
        toggles, toggles, toggles, toggles, toggles, toggles
    ):
        aperture = (
            SubaruAperture(diameter=8.2, grid_size=grid_size, samples_across_pupil=samples)
            if obs_on
            else CircularAperture(diameter=8.2, grid_size=grid_size, samples_across_pupil=samples)
        )
        coronagraph = LyotCoronagraph(enabled=corono_on)
        pywfs = PyramidWFS(modulation_radius=2.0 if mod_on else 0.0, n_modulation_points=4)
        ncpa = ZernikeNCPA(coefficients={6: 0.5}, enabled=ncpa_on)
        sci_det = DetectorModel(enabled=sci_noise_on)
        wfs_det = DetectorModel(enabled=wfs_noise_on)

        system = OpticalSystem(
            aperture=aperture,
            pywfs=pywfs,
            coronagraph=coronagraph,
            ncpa=ncpa,
            science_detector=sci_det,
            wfs_detector=wfs_det,
        )

        res = system.simulate_frame()

        # Assert output contract
        assert res.science_psf.shape == (grid_size, grid_size)
        assert res.direct_psf.shape == (grid_size, grid_size)
        assert res.pywfs_detector.shape == (grid_size, grid_size)
        assert not np.isnan(res.science_psf).any()
        assert not np.isnan(res.pywfs_detector).any()
