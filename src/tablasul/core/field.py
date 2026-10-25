"""Complex electric field representations and wave propagation."""

from __future__ import annotations
from dataclasses import dataclass
from typing import Literal
import numpy as np


@dataclass
class ComplexElectricField:
    """Represents a 2D complex wave field E(x, y) = A * exp(i * phi).

    Attributes:
        data: 2D complex numpy array.
        wavelength: Radiation wavelength in meters.
        pixel_scale: Grid spacing (in meters for pupil planes, radians for focal planes).
        plane: 'pupil' or 'focal'.
        fov_lod: Field of view in lambda/D (if focal plane generated via MFT).
    """

    data: np.ndarray
    wavelength: float
    pixel_scale: float
    plane: Literal["pupil", "focal"] = "pupil"
    fov_lod: float | None = None

    def __post_init__(self) -> None:
        if self.data.ndim != 2:
            raise ValueError(f"Electric field data must be 2D, got shape {self.data.shape}")
        if self.data.shape[0] != self.data.shape[1]:
            raise ValueError(f"Electric field grid must be square, got {self.data.shape}")
        if self.data.dtype not in (np.complex64, np.complex128):
            self.data = self.data.astype(np.complex128)

    @property
    def grid_size(self) -> int:
        """Grid dimension along one axis (N for an N x N grid)."""
        return self.data.shape[0]

    @property
    def intensity(self) -> np.ndarray:
        """Instantaneous intensity |E|^2."""
        return np.abs(self.data) ** 2

    @property
    def phase(self) -> np.ndarray:
        """Phase angle in radians [-pi, pi]."""
        return np.angle(self.data)

    @property
    def total_energy(self) -> float:
        """Total integrated flux sum(|E|^2)."""
        return float(np.sum(self.intensity))

    def apply_phase(self, phase_map: np.ndarray) -> ComplexElectricField:
        """Apply an additive phase delay phi(x, y) in radians."""
        if phase_map.shape != self.data.shape:
            raise ValueError(
                f"Phase map shape {phase_map.shape} does not match field shape {self.data.shape}"
            )
        new_data = self.data * np.exp(1j * phase_map)
        return ComplexElectricField(
            data=new_data,
            wavelength=self.wavelength,
            pixel_scale=self.pixel_scale,
            plane=self.plane,
            fov_lod=self.fov_lod,
        )

    def apply_mask(self, mask: np.ndarray) -> ComplexElectricField:
        """Apply an amplitude transmission mask A(x, y) in [0, 1]."""
        if mask.shape != self.data.shape:
            raise ValueError(
                f"Mask shape {mask.shape} does not match field shape {self.data.shape}"
            )
        new_data = self.data * mask
        return ComplexElectricField(
            data=new_data,
            wavelength=self.wavelength,
            pixel_scale=self.pixel_scale,
            plane=self.plane,
            fov_lod=self.fov_lod,
        )

    def propagate_to_focal(self) -> ComplexElectricField:
        """Propagate pupil field to the focal plane via standard paraxial Fraunhofer 2D FFT.

        Preserves total energy (Parseval's theorem): sum(|E_focal|^2) == sum(|E_pupil|^2).
        """
        if self.plane != "pupil":
            raise ValueError("Can only propagate to focal plane from a pupil plane.")

        n = self.grid_size
        shifted = np.fft.ifftshift(self.data)
        transformed = np.fft.fft2(shifted)
        focal_data = np.fft.fftshift(transformed) / n

        focal_pixel_scale = self.wavelength / (n * self.pixel_scale)

        return ComplexElectricField(
            data=focal_data,
            wavelength=self.wavelength,
            pixel_scale=focal_pixel_scale,
            plane="focal",
        )

    def propagate_to_pupil(self) -> ComplexElectricField:
        """Propagate focal plane field back to the pupil plane via inverse Fraunhofer FFT.

        Preserves total energy (Parseval's theorem).
        """
        if self.plane != "focal":
            raise ValueError("Can only propagate to pupil plane from a focal plane.")

        n = self.grid_size
        shifted = np.fft.ifftshift(self.data)
        transformed = np.fft.ifft2(shifted)
        pupil_data = np.fft.fftshift(transformed) * n

        pupil_pixel_scale = self.wavelength / (n * self.pixel_scale)

        return ComplexElectricField(
            data=pupil_data,
            wavelength=self.wavelength,
            pixel_scale=pupil_pixel_scale,
            plane="pupil",
        )

    def propagate_to_focal_mft(
        self,
        fov_lod: float,
        n_focal: int,
        samples_across_pupil: int,
    ) -> ComplexElectricField:
        """Propagate pupil field to an arbitrary focal plane via Matrix Fourier Transform (MFT).

        Allows decoupling the focal plane field of view (in lambda/D) and pixel count
        from the pupil grid resolution.
        """
        if self.plane != "pupil":
            raise ValueError("Can only propagate to focal plane from a pupil plane.")

        n_pupil = self.grid_size
        mid_p = n_pupil // 2
        mid_f = n_focal // 2

        # Normalized pupil coordinates
        u = (np.arange(n_pupil) - mid_p) / samples_across_pupil
        # Focal coordinates in lambda/D
        xf = (np.arange(n_focal) - mid_f) * (fov_lod / n_focal)

        # 1D Fourier kernels
        kx = np.exp(-2j * np.pi * np.outer(xf, u))
        ky = np.exp(-2j * np.pi * np.outer(xf, u))

        # Direct matrix multiplication
        norm_factor = fov_lod / (samples_across_pupil * n_focal)
        focal_data = norm_factor * (ky @ self.data @ kx.T)

        focal_pixel_scale = (fov_lod / n_focal) * (self.wavelength / (samples_across_pupil * self.pixel_scale))

        return ComplexElectricField(
            data=focal_data,
            wavelength=self.wavelength,
            pixel_scale=focal_pixel_scale,
            plane="focal",
            fov_lod=fov_lod,
        )
