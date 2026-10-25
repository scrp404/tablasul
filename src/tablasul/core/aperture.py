"""Telescope entrance pupil and aperture masks."""

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Union
import numpy as np
from scipy.ndimage import zoom
from astropy.io import fits
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


class PupilImageAperture(Aperture):
    """Entrance pupil loaded from a file (FITS, NPY, PNG, etc.) or preloaded 2D array.

    Enables using experimentally measured Subaru pupil files, interferometer maps,
    or CAD aperture masks.

    Attributes:
        image_or_path: Path to file (.fits, .npy, .png, .tif) or 2D numpy array.
        diameter: Physical aperture diameter in meters (default: 8.2 m).
        grid_size: Simulation grid size N x N (default: matches loaded image).
        samples_across_pupil: Number of pixels across pupil diameter (if None, auto-estimated).
    """

    def __init__(
        self,
        image_or_path: Union[str, Path, np.ndarray],
        diameter: float = 8.2,
        grid_size: Optional[int] = None,
        samples_across_pupil: Optional[int] = None,
    ) -> None:
        raw_image = self._load_data(image_or_path)

        # Normalize transmission to [0, 1]
        max_val = np.max(raw_image)
        if max_val > 0:
            raw_image = raw_image / max_val
        raw_image = np.clip(raw_image, 0.0, 1.0)

        # Determine grid size
        if grid_size is None:
            grid_size = raw_image.shape[0]

        # Resample if needed
        if raw_image.shape != (grid_size, grid_size):
            zoom_factor = grid_size / raw_image.shape[0]
            self._mask = zoom(raw_image, zoom_factor, order=1)
            self._mask = np.clip(self._mask, 0.0, 1.0)
        else:
            self._mask = raw_image.astype(np.float64)

        # Estimate samples across pupil if not provided
        if samples_across_pupil is None:
            nonzeros = np.where(self._mask > 0.05)
            if len(nonzeros[0]) > 0:
                y_span = np.max(nonzeros[0]) - np.min(nonzeros[0])
                x_span = np.max(nonzeros[1]) - np.min(nonzeros[1])
                samples_across_pupil = max(int(max(y_span, x_span)), 1)
            else:
                samples_across_pupil = grid_size // 2

        super().__init__(
            diameter=diameter,
            grid_size=grid_size,
            samples_across_pupil=samples_across_pupil,
        )

    def _load_data(self, source: Union[str, Path, np.ndarray]) -> np.ndarray:
        if isinstance(source, np.ndarray):
            if source.ndim != 2:
                raise ValueError(f"Pupil image array must be 2D, got shape {source.shape}")
            return source.astype(np.float64)

        path = Path(source)
        if not path.exists():
            raise FileNotFoundError(f"Pupil file not found: {path}")

        suffix = path.suffix.lower()
        if suffix in (".fits", ".fit", ".fts"):
            data = fits.getdata(path)
            # If 3D cube, take first slice
            if data.ndim == 3:
                data = data[0]
            return data.astype(np.float64)
        elif suffix == ".npy":
            return np.load(path).astype(np.float64)
        else:
            # Attempt image load with PIL
            from PIL import Image
            img = Image.open(path).convert("L")
            return np.array(img, dtype=np.float64)

    def transmission_mask(self) -> np.ndarray:
        return self._mask.copy()
