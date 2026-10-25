"""Synchronized multi-dimensional data cube generation engine."""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Literal
import numpy as np
from ..core.atmosphere import FrozenFlowAtmosphere
from .simulator import OpticalSystem


@dataclass
class SynchronizedDataset:
    """Container for synchronized multi-dimensional simulation data cubes."""

    wfs_cube: np.ndarray  # Shape: (T_wfs, 4, H, W) or (T_wfs, 2, H, W)
    psf_cube: np.ndarray  # Shape: (T_sci, H_focal, W_focal)
    ground_truth_phase: np.ndarray  # Shape: (T_wfs, H_pupil, W_pupil)
    timestamps_wfs: np.ndarray  # (T_wfs,) in seconds
    timestamps_psf: np.ndarray  # (T_sci,) in seconds
    metadata: Dict


class SynchronizedCubeEngine:
    """Generates synchronized data cubes between WFS and Science camera.

    Handles cadence oversampling where multiple WFS frames integrate into a
    single longer science camera exposure.
    """

    def __init__(
        self,
        system: OpticalSystem,
        atmosphere: FrozenFlowAtmosphere | None = None,
        wfs_frame_rate: float = 1000.0,  # Hz (1 ms per frame)
        science_exposure_time: float = 0.01,  # seconds (10 ms per frame)
        wfs_output_mode: Literal["quadrants", "slopes"] = "quadrants",
    ) -> None:
        self.system = system
        self.atmosphere = atmosphere
        self.wfs_frame_rate = wfs_frame_rate
        self.science_exposure_time = science_exposure_time
        self.wfs_output_mode = wfs_output_mode

    def generate_cubes(
        self,
        n_science_frames: int = 10,
        seed: int | None = None,
    ) -> SynchronizedDataset:
        """Run simulation over time and assemble synchronized data cubes."""
        dt_wfs = 1.0 / self.wfs_frame_rate
        substeps_per_science = max(1, int(round(self.science_exposure_time / dt_wfs)))
        total_wfs_frames = n_science_frames * substeps_per_science

        grid_size = self.system.aperture.grid_size
        samples_pupil = self.system.aperture.samples_across_pupil
        pixel_scale = self.system.aperture.pixel_scale

        if self.atmosphere is not None:
            self.atmosphere.initialize(
                pupil_grid_size=grid_size,
                pixel_scale=pixel_scale,
                target_wavelength=self.system.wfs_wavelength,
                seed=seed,
            )

        # Allocate arrays
        sample_res = self.system.simulate_frame()
        q_shape = sample_res.pywfs_quadrants[0].shape
        focal_shape = sample_res.science_psf.shape

        if self.wfs_output_mode == "quadrants":
            wfs_cube = np.zeros((total_wfs_frames, 4, q_shape[0], q_shape[1]), dtype=np.float32)
        else:
            wfs_cube = np.zeros((total_wfs_frames, 2, q_shape[0], q_shape[1]), dtype=np.float32)

        psf_cube = np.zeros((n_science_frames, focal_shape[0], focal_shape[1]), dtype=np.float32)
        phase_cube = np.zeros((total_wfs_frames, grid_size, grid_size), dtype=np.float32)
        timestamps_wfs = np.zeros(total_wfs_frames, dtype=np.float64)
        timestamps_psf = np.zeros(n_science_frames, dtype=np.float64)

        wfs_idx = 0
        current_time = 0.0

        for sci_idx in range(n_science_frames):
            science_accumulator = np.zeros(focal_shape, dtype=np.float64)
            timestamps_psf[sci_idx] = current_time

            for _ in range(substeps_per_science):
                if self.atmosphere is not None:
                    phase = self.atmosphere.get_phase_screen(current_time, grid_size)
                else:
                    phase = np.zeros((grid_size, grid_size), dtype=np.float64)

                res = self.system.simulate_frame(phase_aberration=phase)

                # Store WFS frame
                if self.wfs_output_mode == "quadrants":
                    wfs_cube[wfs_idx, 0] = res.pywfs_quadrants[0]
                    wfs_cube[wfs_idx, 1] = res.pywfs_quadrants[1]
                    wfs_cube[wfs_idx, 2] = res.pywfs_quadrants[2]
                    wfs_cube[wfs_idx, 3] = res.pywfs_quadrants[3]
                else:
                    wfs_cube[wfs_idx, 0] = res.slopes[0]
                    wfs_cube[wfs_idx, 1] = res.slopes[1]

                phase_cube[wfs_idx] = phase
                timestamps_wfs[wfs_idx] = current_time

                # Integrate science intensity
                science_accumulator += res.science_psf

                wfs_idx += 1
                current_time += dt_wfs

            psf_cube[sci_idx] = science_accumulator / substeps_per_science

        metadata = {
            "wfs_frame_rate_hz": self.wfs_frame_rate,
            "science_exposure_time_s": self.science_exposure_time,
            "substeps_per_science_frame": substeps_per_science,
            "n_science_frames": n_science_frames,
            "total_wfs_frames": total_wfs_frames,
            "wfs_wavelength_m": self.system.wfs_wavelength,
            "science_wavelength_m": self.system.science_wavelength,
            "wfs_output_mode": self.wfs_output_mode,
        }

        return SynchronizedDataset(
            wfs_cube=wfs_cube,
            psf_cube=psf_cube,
            ground_truth_phase=phase_cube,
            timestamps_wfs=timestamps_wfs,
            timestamps_psf=timestamps_psf,
            metadata=metadata,
        )
