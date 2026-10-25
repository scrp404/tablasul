"""Tests for spatiotemporal synchronization and multi-dimensional cube generation."""

import numpy as np
import pytest
from tablasul.core.aperture import CircularAperture
from tablasul.core.atmosphere import KolmogorovScreen, FrozenFlowAtmosphere
from tablasul.engine.simulator import OpticalSystem
from tablasul.engine.sync import SynchronizedCubeEngine


def test_synchronized_cubes_shapes_and_oversampling():
    """Verify that SynchronizedCubeEngine produces correct cube shapes and oversampling."""
    system = OpticalSystem(
        aperture=CircularAperture(diameter=8.2, grid_size=64, samples_across_pupil=32),
    )
    # 1000 Hz WFS, 4 ms science exposure => 4 WFS frames per 1 science frame
    engine = SynchronizedCubeEngine(
        system=system,
        atmosphere=None,
        wfs_frame_rate=1000.0,
        science_exposure_time=0.004,
        wfs_output_mode="quadrants",
    )
    n_sci = 5
    dataset = engine.generate_cubes(n_science_frames=n_sci)

    # 5 science frames * 4 = 20 WFS frames
    assert dataset.psf_cube.shape[0] == n_sci
    assert dataset.wfs_cube.shape[0] == 20
    assert dataset.ground_truth_phase.shape[0] == 20
    assert dataset.wfs_cube.shape[1] == 4  # 4 quadrants


def test_synchronized_cubes_slopes_mode_and_atmosphere():
    """Verify synchronized cubes in slopes mode with evolving frozen-flow atmosphere."""
    system = OpticalSystem(
        aperture=CircularAperture(diameter=8.2, grid_size=64, samples_across_pupil=32),
    )
    atmos = FrozenFlowAtmosphere(
        screen_generator=KolmogorovScreen(r0=0.20, L0=25.0),
        wind_velocity=(15.0, 0.0),
    )
    engine = SynchronizedCubeEngine(
        system=system,
        atmosphere=atmos,
        wfs_frame_rate=500.0,
        science_exposure_time=0.002,  # 1-to-1 ratio
        wfs_output_mode="slopes",
    )
    dataset = engine.generate_cubes(n_science_frames=3, seed=42)

    assert dataset.wfs_cube.shape[0] == 3
    assert dataset.wfs_cube.shape[1] == 2  # Sx, Sy
    assert dataset.psf_cube.shape[0] == 3
    # Ground truth phase should evolve over time (nonzero variance)
    p0 = dataset.ground_truth_phase[0]
    p1 = dataset.ground_truth_phase[1]
    assert not np.allclose(p0, p1)
