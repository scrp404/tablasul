"""Telescope entrance pupil and aperture masks."""

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass
import numpy as np
from .field import ComplexElectricField


@dataclass
class Aperture(ABC):
    """Abstract base class for pupil apertures."""

    diameter: float
    grid_size: int
    samples_across_pupil: int

    @property
    def pixel_scale(self) -> float:
        """Physical pupil grid spacing (meters per pixel)."""
        return self.diameter / self.samples_across_pupil

    @abstractmethod
    def transmission_mask(self) -> np.ndarray:
        """Return 2D float array in [0, 1] representing pupil transmission."""
        pass

    def create_field(
        self,
        wavelength: float,
        total_flux: float = 1.0,
        phase_aberration: np.ndarray | None = None,
    ) -> ComplexElectricField:
        """Create an incident complex electric field initialized with the aperture mask."""
        mask = self.transmission_mask()
        pupil_flux = np.sum(mask**2)
        amplitude = mask * np.sqrt(total_flux / max(pupil_flux, 1e-12))

        if phase_aberration is not None:
            data = amplitude * np.exp(1j * phase_aberration)
        else:
            data = amplitude.astype(np.complex128)

        return ComplexElectricField(
            data=data,
            wavelength=wavelength,
            pixel_scale=self.pixel_scale,
            plane="pupil",
        )


@dataclass
class CircularAperture(Aperture):
    """Unobstructed circular aperture ("toy model" baseline)."""

    diameter: float = 8.2
    grid_size: int = 128
    samples_across_pupil: int = 64

    def transmission_mask(self) -> np.ndarray:
        n = self.grid_size
        center = n // 2
        y, x = np.ogrid[:n, :n]
        r = np.hypot(x - center, y - center)
        radius_pixels = self.samples_across_pupil / 2.0
        return (r <= radius_pixels).astype(np.float64)


@dataclass
class SubaruAperture(Aperture):
    """Subaru Telescope 8.2m aperture with central obscuration and 4 spiders."""

    diameter: float = 8.2
    grid_size: int = 128
    samples_across_pupil: int = 64
    central_obscuration_ratio: float = 0.29
    spider_width: float = 0.22  # meters
    enabled_spiders: bool = True
    enabled_obscuration: bool = True

    def transmission_mask(self) -> np.ndarray:
        n = self.grid_size
        center = n // 2
        y, x = np.ogrid[:n, :n]
        coords_x = (x - center) * self.pixel_scale
        coords_y = (y - center) * self.pixel_scale
        r = np.hypot(coords_x, coords_y)
        radius = self.diameter / 2.0

        mask = (r <= radius).astype(np.float64)

        if self.enabled_obscuration and self.central_obscuration_ratio > 0:
            obs_radius = radius * self.central_obscuration_ratio
            mask[r < obs_radius] = 0.0

        if self.enabled_spiders and self.spider_width > 0:
            half_w = self.spider_width / 2.0
            d1 = np.abs(coords_y - coords_x) / np.sqrt(2)
            d2 = np.abs(coords_y + coords_x) / np.sqrt(2)
            spider_mask = (d1 < half_w) | (d2 < half_w)
            mask[spider_mask] = 0.0

        return mask
