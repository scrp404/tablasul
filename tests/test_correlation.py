"""Tests verifying physical and statistical coupling between synchronized WFS telemetry and PSF frames."""

import numpy as np
import pytest
from scipy.stats import pearsonr
from tablasul.core.aperture import CircularAperture, SubaruAperture
from tablasul.core.atmosphere import KolmogorovScreen
from tablasul.core.source import PointSource
from tablasul.optics.pywfs import PyramidWFS
from tablasul.optics.coronagraph import LyotCoronagraph
from tablasul.detector.camera import DetectorModel
from tablasul.engine.simulator import OpticalSystem


def test_tip_tilt_wfs_psf_centroid_correlation():
    """Verify that WFS slopes and PSF centroids have >0.99 correlation under random tip/tilt."""
    grid_size = 128
    samples_pupil = 32
    aperture = CircularAperture(diameter=8.2, grid_size=grid_size, samples_across_pupil=samples_pupil)
    pywfs = PyramidWFS(modulation_radius=2.0)
    system = OpticalSystem(aperture=aperture, pywfs=pywfs, sources=[PointSource(total_flux=1e7)])

    n_samples = 25
    rng = np.random.default_rng(42)
    tilts_x = rng.uniform(-0.5, 0.5, size=n_samples)
    tilts_y = rng.uniform(-0.5, 0.5, size=n_samples)

    y_grid, x_grid = np.indices((grid_size, grid_size))
    x_norm = (x_grid - grid_size // 2) / (samples_pupil / 2.0)
    y_norm = (y_grid - grid_size // 2) / (samples_pupil / 2.0)

    mean_sx = []
    mean_sy = []
    psf_centroid_x = []
    psf_centroid_y = []

    focal_y, focal_x = np.indices((grid_size, grid_size))

    for tx, ty in zip(tilts_x, tilts_y):
        phase = tx * x_norm + ty * y_norm
        res = system.simulate_frame(phase_aberration=phase)

        sx, sy, mask = res.slopes
        valid = mask > 0.5
        mean_sx.append(np.mean(sx[valid]))
        mean_sy.append(np.mean(sy[valid]))

        psf = res.direct_psf
        norm = np.sum(psf)
        cx = np.sum(focal_x * psf) / norm
        cy = np.sum(focal_y * psf) / norm
        psf_centroid_x.append(cx)
        psf_centroid_y.append(cy)

    corr_x, _ = pearsonr(mean_sx, psf_centroid_x)
    corr_y, _ = pearsonr(mean_sy, psf_centroid_y)

    assert abs(corr_x) > 0.99
    assert abs(corr_y) > 0.99


def test_independent_realizations_frame_matching():
    """Verify that contemporaneous WFS frames have distinctively higher correlation with their matching PSF than mismatched frames."""
    grid_size = 128
    samples_pupil = 32
    system = OpticalSystem(
        aperture=CircularAperture(diameter=8.2, grid_size=grid_size, samples_across_pupil=samples_pupil),
        pywfs=PyramidWFS(modulation_radius=2.0),
        sources=[PointSource(total_flux=1e7)],
    )
    screen_gen = KolmogorovScreen(r0=0.20, L0=25.0)

    n_realizations = 20
    wfs_frames = []
    psf_frames = []

    for i in range(n_realizations):
        phase = screen_gen.generate_screen(
            grid_size=grid_size,
            pixel_scale=8.2 / samples_pupil,
            target_wavelength=800e-9,
            seed=2000 + i,
        )
        res = system.simulate_frame(phase_aberration=phase)
        sx, sy, _ = res.slopes
        wfs_frames.append(np.stack([sx, sy], axis=0))
        psf_frames.append(res.direct_psf)

    wfs = np.array(wfs_frames)
    psf = np.array(psf_frames)

    wfs_flat = wfs.reshape(n_realizations, -1)
    psf_flat = psf.reshape(n_realizations, -1)

    dist_w = np.linalg.norm(wfs_flat[:, None, :] - wfs_flat[None, :, :], axis=2)
    dist_p = np.linalg.norm(psf_flat[:, None, :] - psf_flat[None, :, :], axis=2)

    sim_matrix = np.zeros((n_realizations, n_realizations))
    for i in range(n_realizations):
        for j in range(n_realizations):
            sim_matrix[i, j] = pearsonr(dist_w[i], dist_p[j])[0]

    diag_vals = np.diag(sim_matrix)
    off_diag = sim_matrix[~np.eye(n_realizations, dtype=bool)]

    # True matching pairs must have significantly higher correlation than mismatched pairs
    assert np.mean(diag_vals) > 0.50
    assert np.mean(off_diag) < 0.15

    # Top-1 accuracy must be >= 90%
    ranks = [np.where(np.argsort(sim_matrix[i])[::-1] == i)[0][0] for i in range(n_realizations)]
    top1 = np.mean(np.array(ranks) == 0)
    assert top1 >= 0.90
