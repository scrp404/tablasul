"""Astronomical and calibration optical sources and filter bands."""

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import List, Tuple
import numpy as np


@dataclass
class FilterBand:
    """Spectral filter bandpass profile.

    Attributes:
        central_wavelength: Central wavelength in meters.
        bandwidth: Bandwidth in meters.
        name: Name of the filter band (e.g. 'R', 'I', 'J', 'H').
    """

    central_wavelength: float
    bandwidth: float
    name: str = "custom"

    @classmethod
    def standard(cls, name: str) -> FilterBand:
        """Create a standard astronomical photometric band."""
        name_upper = name.upper()
        standards = {
            "R": (650e-9, 100e-9),
            "I": (800e-9, 150e-9),
            "J": (1250e-9, 200e-9),
            "H": (1650e-9, 300e-9),
            "K": (2200e-9, 400e-9),
        }
        if name_upper not in standards:
            raise ValueError(f"Unknown filter band '{name}'. Standard bands: {list(standards.keys())}")
        cw, bw = standards[name_upper]
        return cls(central_wavelength=cw, bandwidth=bw, name=name_upper)

    def sample_channels(self, n_channels: int = 1) -> List[Tuple[float, float]]:
        """Return list of (wavelength, weight) pairs across the band."""
        if n_channels <= 1:
            return [(self.central_wavelength, 1.0)]
        half_bw = self.bandwidth / 2.0
        wavelengths = np.linspace(
            self.central_wavelength - half_bw,
            self.central_wavelength + half_bw,
            n_channels,
        )
        # Uniform top-hat weighting normalized to 1
        weight = 1.0 / n_channels
        return [(float(wl), weight) for wl in wavelengths]


@dataclass
class Source(ABC):
    """Abstract base class for optical sources."""

    total_flux: float = 1.0

    @abstractmethod
    def get_phase_tilt(
        self,
        grid_size: int,
        samples_across_pupil: int,
        wavelength: float,
    ) -> np.ndarray:
        """Compute the 2D linear phase tilt across the entrance pupil in radians."""
        pass


@dataclass
class PointSource(Source):
    """On-axis or off-axis astronomical point source (e.g. primary star).

    Attributes:
        total_flux: Total photon flux.
        offset_lod: Offset on sky in units of lambda/D as (offset_x, offset_y).
    """

    total_flux: float = 1.0
    offset_lod: Tuple[float, float] = (0.0, 0.0)

    def get_phase_tilt(
        self,
        grid_size: int,
        samples_across_pupil: int,
        wavelength: float,
    ) -> np.ndarray:
        ox, oy = self.offset_lod
        if ox == 0.0 and oy == 0.0:
            return np.zeros((grid_size, grid_size), dtype=np.float64)

        mid = grid_size // 2
        y, x = np.ogrid[:grid_size, :grid_size]
        # In normalized FFT coordinates, 1 lambda/D across pupil diameter D corresponds
        # to a 2*pi phase ramp over the pupil diameter (samples_across_pupil pixels).
        # tilt = 2*pi * (ox * x_pupil + oy * y_pupil) / samples_across_pupil
        tilt = 2.0 * np.pi * (ox * (x - mid) + oy * (y - mid)) / samples_across_pupil
        return tilt


@dataclass
class CompanionSource(Source):
    """Off-axis faint companion (e.g. exoplanet or binary companion).

    Attributes:
        contrast: Flux contrast ratio relative to primary star (e.g. 1e-4).
        separation_lod: Radial separation in units of lambda/D.
        position_angle_deg: Position angle on sky in degrees (0 = +x, 90 = +y).
    """

    contrast: float = 1e-4
    separation_lod: float = 4.0
    position_angle_deg: float = 0.0
    total_flux: float = field(init=False)

    def __post_init__(self) -> None:
        self.total_flux = self.contrast

    def get_phase_tilt(
        self,
        grid_size: int,
        samples_across_pupil: int,
        wavelength: float,
    ) -> np.ndarray:
        rad = np.radians(self.position_angle_deg)
        ox = self.separation_lod * np.cos(rad)
        oy = self.separation_lod * np.sin(rad)
        mid = grid_size // 2
        y, x = np.ogrid[:grid_size, :grid_size]
        tilt = 2.0 * np.pi * (ox * (x - mid) + oy * (y - mid)) / samples_across_pupil
        return tilt


@dataclass
class CalibrationSource(Source):
    """Internal laboratory calibration source (monochromatic laser or flat lamp).

    Attributes:
        wavelength: Exact laser emission wavelength in meters.
        total_flux: Calibration flux.
    """

    wavelength: float = 635e-9
    total_flux: float = 1.0

    def get_phase_tilt(
        self,
        grid_size: int,
        samples_across_pupil: int,
        wavelength: float,
    ) -> np.ndarray:
        return np.zeros((grid_size, grid_size), dtype=np.float64)
