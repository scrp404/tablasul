"""End-to-end optical simulation engine."""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Union
import numpy as np
from ..core.aperture import Aperture, CircularAperture
from ..core.field import ComplexElectricField
from ..core.source import Source, PointSource, FilterBand
from ..optics.pywfs import PyramidWFS
from ..optics.coronagraph import LyotCoronagraph
from ..optics.ncpa import StaticPhaseNCPA, ZernikeNCPA
from ..detector.camera import DetectorModel


@dataclass
class SimulationResult:
    """Container for a single simulated frame."""

    science_psf: np.ndarray
    direct_psf: np.ndarray
    pywfs_detector: np.ndarray
    pywfs_quadrants: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
    slopes: Tuple[np.ndarray, np.ndarray, np.ndarray]
    incident_phase: Optional[np.ndarray] = None


class OpticalSystem:
    """Master orchestrator for first-principles optical propagation.

    Attributes:
        aperture: The entrance pupil aperture model.
        pywfs: Pyramid wavefront sensor model.
        coronagraph: Lyot coronagraph model in the science arm.
        sources: List of mutually incoherent optical sources (stars, companions, lasers).
        ncpa: Non-Common Path Aberrations model applied to the science arm.
        science_detector: Noise model for the science camera.
        wfs_detector: Noise model for the WFS camera.
        wfs_wavelength: WFS sensing wavelength in meters (default: 800 nm, visible).
        science_wavelength: Science imaging wavelength in meters (default: 1650 nm, NIR H-band).
    """

    def __init__(
        self,
        aperture: Optional[Aperture] = None,
        pywfs: Optional[PyramidWFS] = None,
        coronagraph: Optional[LyotCoronagraph] = None,
        sources: Optional[List[Source]] = None,
        ncpa: Optional[Union[StaticPhaseNCPA, ZernikeNCPA]] = None,
        science_detector: Optional[DetectorModel] = None,
        wfs_detector: Optional[DetectorModel] = None,
        wavelength: Optional[float] = None,
        wfs_wavelength: float = 800e-9,
        science_wavelength: float = 1650e-9,
    ) -> None:
        self.aperture = aperture if aperture is not None else CircularAperture()
        self.pywfs = pywfs if pywfs is not None else PyramidWFS()
        self.coronagraph = coronagraph if coronagraph is not None else LyotCoronagraph(enabled=False)
        self.sources = sources if sources is not None else [PointSource(total_flux=1.0)]
        self.ncpa = ncpa
        self.science_detector = science_detector if science_detector is not None else DetectorModel(enabled=False)
        self.wfs_detector = wfs_detector if wfs_detector is not None else DetectorModel(enabled=False)

        if wavelength is not None:
            self.wfs_wavelength = wavelength
            self.science_wavelength = wavelength
        else:
            self.wfs_wavelength = wfs_wavelength
            self.science_wavelength = science_wavelength

    def simulate_frame(
        self,
        phase_aberration: Optional[np.ndarray] = None,
        ncpa_phase: Optional[np.ndarray] = None,
        exposure_time: float = 1.0,
    ) -> SimulationResult:
        """Simulate a single optical frame across all configured sources."""
        grid_size = self.aperture.grid_size
        samples_pupil = self.aperture.samples_across_pupil
        pupil_radius_pixels = samples_pupil / 2.0

        total_direct_psf = np.zeros((grid_size, grid_size), dtype=np.float64)
        total_science_psf = np.zeros((grid_size, grid_size), dtype=np.float64)
        total_wfs_detector = np.zeros((grid_size, grid_size), dtype=np.float64)

        # Evaluate static NCPA if configured
        active_ncpa = ncpa_phase
        if active_ncpa is None and self.ncpa is not None:
            if isinstance(self.ncpa, ZernikeNCPA):
                active_ncpa = self.ncpa.generate_phase(grid_size, pupil_radius_pixels)
            elif isinstance(self.ncpa, StaticPhaseNCPA):
                active_ncpa = self.ncpa.get_phase()

        I1_tot = None
        I2_tot = None
        I3_tot = None
        I4_tot = None

        for src in self.sources:
            source_tilt_wfs = src.get_phase_tilt(
                grid_size=grid_size,
                samples_across_pupil=samples_pupil,
                wavelength=self.wfs_wavelength,
            )
            source_tilt_sci = src.get_phase_tilt(
                grid_size=grid_size,
                samples_across_pupil=samples_pupil,
                wavelength=self.science_wavelength,
            )

            total_phase_wfs = source_tilt_wfs
            total_phase_sci = source_tilt_sci

            if phase_aberration is not None:
                phase_scale_sci = self.wfs_wavelength / self.science_wavelength
                total_phase_wfs = total_phase_wfs + phase_aberration
                total_phase_sci = total_phase_sci + (phase_aberration * phase_scale_sci)

            if active_ncpa is not None:
                total_phase_sci = total_phase_sci + active_ncpa

            # Science arm
            sci_pupil_field = self.aperture.create_field(
                wavelength=self.science_wavelength,
                total_flux=src.total_flux,
                phase_aberration=total_phase_sci,
            )
            sci_direct_focal = sci_pupil_field.propagate_to_focal()
            total_direct_psf += sci_direct_focal.intensity

            sci_corono_focal = self.coronagraph.propagate(
                focal_field=sci_direct_focal,
                samples_across_pupil=samples_pupil,
            )
            total_science_psf += sci_corono_focal.intensity

            # WFS arm
            wfs_pupil_field = self.aperture.create_field(
                wavelength=self.wfs_wavelength,
                total_flux=src.total_flux,
                phase_aberration=total_phase_wfs,
            )
            wfs_det, q1, q2, q3, q4 = self.pywfs.sense(
                pupil_field=wfs_pupil_field,
                pupil_radius_pixels=pupil_radius_pixels,
            )
            total_wfs_detector += wfs_det

            if I1_tot is None:
                I1_tot = q1.copy()
                I2_tot = q2.copy()
                I3_tot = q3.copy()
                I4_tot = q4.copy()
            else:
                I1_tot += q1
                I2_tot += q2
                I3_tot += q3
                I4_tot += q4

        # Apply detector noise if enabled
        final_science_psf = self.science_detector.apply(total_science_psf, exposure_time=exposure_time)
        final_wfs_detector = self.wfs_detector.apply(total_wfs_detector, exposure_time=exposure_time)

        if self.wfs_detector.enabled:
            I1_tot = self.wfs_detector.apply(I1_tot, exposure_time=exposure_time)
            I2_tot = self.wfs_detector.apply(I2_tot, exposure_time=exposure_time)
            I3_tot = self.wfs_detector.apply(I3_tot, exposure_time=exposure_time)
            I4_tot = self.wfs_detector.apply(I4_tot, exposure_time=exposure_time)

        Sx, Sy, mask = self.pywfs.compute_slopes(I1_tot, I2_tot, I3_tot, I4_tot)

        return SimulationResult(
            science_psf=final_science_psf,
            direct_psf=total_direct_psf,
            pywfs_detector=final_wfs_detector,
            pywfs_quadrants=(I1_tot, I2_tot, I3_tot, I4_tot),
            slopes=(Sx, Sy, mask),
            incident_phase=phase_aberration,
        )
