"""Pyramid Wavefront Sensor (PyWFS) optical model and slope computation."""

from __future__ import annotations
from dataclasses import dataclass
from typing import Tuple
import numpy as np
from ..core.field import ComplexElectricField


@dataclass
class PyramidWFS:
    """Optical model of a 4-faceted Pyramid Wavefront Sensor.

    Attributes:
        separation_ratio: Separation distance between pupil centers in pupil diameters (default: 2.0).
        modulation_radius: Dynamic modulation radius in units of lambda/D (0 for unmodulated).
        n_modulation_points: Number of discrete integration points along the modulation circle.
    """

    separation_ratio: float = 2.0
    modulation_radius: float = 0.0
    n_modulation_points: int = 16

    def generate_prism_phase_mask(self, grid_size: int, pupil_radius_pixels: float) -> np.ndarray:
        """Create the 4-faceted pyramid phase mask in the focal plane."""
        n = grid_size
        mid = n // 2
        y, x = np.ogrid[:n, :n]
        u = x - mid
        v = y - mid

        shift_pixels = pupil_radius_pixels * self.separation_ratio
        tilt_factor = 2.0 * np.pi * shift_pixels / n

        # Four-faceted pyramid phase profile
        phase_mask = tilt_factor * (np.abs(u) + np.abs(v))
        return phase_mask

    def sense(
        self,
        pupil_field: ComplexElectricField,
        pupil_radius_pixels: float,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Propagate pupil field through the PyWFS and re-image the 4 pupils.

        Returns:
            Tuple of:
                full_detector_image: 2D array of the detector containing all 4 pupils
                I1, I2, I3, I4: Extracted 2D pupil quadrant images
        """
        n = pupil_field.grid_size
        prism_phase = self.generate_prism_phase_mask(n, pupil_radius_pixels)

        if self.modulation_radius == 0.0 or self.n_modulation_points <= 1:
            focal_field = pupil_field.propagate_to_focal()
            modulated_focal = focal_field.apply_phase(prism_phase)
            reimaged = modulated_focal.propagate_to_pupil()
            full_detector = reimaged.intensity
        else:
            full_detector = np.zeros((n, n), dtype=np.float64)
            thetas = np.linspace(0, 2 * np.pi, self.n_modulation_points, endpoint=False)
            mid = n // 2
            y, x = np.ogrid[:n, :n]
            xp = (x - mid) / pupil_radius_pixels
            yp = (y - mid) / pupil_radius_pixels

            for theta in thetas:
                tilt = np.pi * self.modulation_radius * (xp * np.cos(theta) + yp * np.sin(theta))
                tilted_pupil = pupil_field.apply_phase(tilt)
                focal_field = tilted_pupil.propagate_to_focal()
                modulated_focal = focal_field.apply_phase(prism_phase)
                reimaged = modulated_focal.propagate_to_pupil()
                full_detector += reimaged.intensity

            full_detector /= self.n_modulation_points

        shift = int(round(pupil_radius_pixels * self.separation_ratio))
        mid = n // 2
        r = int(round(pupil_radius_pixels * 1.15))

        def extract_quad(cy: int, cx: int) -> np.ndarray:
            y0, y1 = max(0, cy - r), min(n, cy + r)
            x0, x1 = max(0, cx - r), min(n, cx + r)
            return full_detector[y0:y1, x0:x1]

        # In numpy ifft2, positive focal phase tilt shifts conjugate coordinates towards negative indices:
        # Focal Q1 (u>0, v>0) lands at row = mid - shift, col = mid - shift
        # Focal Q2 (u<0, v>0) lands at row = mid - shift, col = mid + shift
        # Focal Q3 (u<0, v<0) lands at row = mid + shift, col = mid + shift
        # Focal Q4 (u>0, v<0) lands at row = mid + shift, col = mid - shift
        I1 = extract_quad(mid - shift, mid - shift)  # Focal Q1 (u>0, v>0)
        I2 = extract_quad(mid - shift, mid + shift)  # Focal Q2 (u<0, v>0)
        I3 = extract_quad(mid + shift, mid + shift)  # Focal Q3 (u<0, v<0)
        I4 = extract_quad(mid + shift, mid - shift)  # Focal Q4 (u>0, v<0)

        return full_detector, I1, I2, I3, I4

    def compute_slopes(
        self,
        I1: np.ndarray,
        I2: np.ndarray,
        I3: np.ndarray,
        I4: np.ndarray,
        flux_threshold_ratio: float = 0.05,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Compute normalized quad-cell slopes (Sx, Sy).

        Returns:
            Tuple of:
                Sx: Horizontal slope map
                Sy: Vertical slope map
                valid_mask: Boolean mask of valid illuminated subapertures
        """
        I_tot = I1 + I2 + I3 + I4
        max_flux = np.max(I_tot) if np.max(I_tot) > 0 else 1.0
        valid_mask = I_tot >= (flux_threshold_ratio * max_flux)

        Sx = np.zeros_like(I_tot)
        Sy = np.zeros_like(I_tot)

        Sx[valid_mask] = ((I1[valid_mask] + I4[valid_mask]) - (I2[valid_mask] + I3[valid_mask])) / I_tot[valid_mask]
        Sy[valid_mask] = ((I1[valid_mask] + I2[valid_mask]) - (I3[valid_mask] + I4[valid_mask])) / I_tot[valid_mask]

        return Sx, Sy, valid_mask
