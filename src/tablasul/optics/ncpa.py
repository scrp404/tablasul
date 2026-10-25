"""Non-Common Path Aberrations (NCPA) modeling."""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass
class StaticPhaseNCPA:
    """Static 2D Non-Common Path Aberration phase map.

    Attributes:
        phase_map: 2D array of phase aberration in radians.
        enabled: Boolean toggle.
    """

    phase_map: np.ndarray
    enabled: bool = True

    def get_phase(self) -> np.ndarray:
        if not self.enabled:
            return np.zeros_like(self.phase_map)
        return self.phase_map


@dataclass
class ZernikeNCPA:
    """Generates NCPA from low-order Zernike polynomials.

    Attributes:
        coefficients: Dictionary mapping Noll index (or name) to RMS amplitude in radians.
        enabled: Boolean toggle.
    """

    coefficients: dict[int, float]
    enabled: bool = True

    def generate_phase(self, grid_size: int, pupil_radius_pixels: float) -> np.ndarray:
        """Evaluate Zernike phase screen across the pupil."""
        if not self.enabled:
            return np.zeros((grid_size, grid_size), dtype=np.float64)

        mid = grid_size // 2
        y, x = np.ogrid[:grid_size, :grid_size]
        rho = np.hypot(x - mid, y - mid) / pupil_radius_pixels
        theta = np.arctan2(y - mid, x - mid)
        pupil_mask = rho <= 1.0

        phase = np.zeros((grid_size, grid_size), dtype=np.float64)

        # Standard Noll indices:
        # Z4: Defocus = sqrt(3) * (2*rho^2 - 1)
        # Z5: Oblique Astigmatism = sqrt(6) * rho^2 * sin(2*theta)
        # Z6: Vertical Astigmatism = sqrt(6) * rho^2 * cos(2*theta)
        # Z7: Vertical Coma = sqrt(8) * (3*rho^3 - 2*rho) * sin(theta)
        # Z8: Horizontal Coma = sqrt(8) * (3*rho^3 - 2*rho) * cos(theta)
        # Z11: Spherical = sqrt(5) * (6*rho^4 - 6*rho^2 + 1)
        for j, amp in self.coefficients.items():
            if amp == 0.0:
                continue
            if j == 4:
                z = np.sqrt(3) * (2 * rho**2 - 1)
            elif j == 5:
                z = np.sqrt(6) * (rho**2) * np.sin(2 * theta)
            elif j == 6:
                z = np.sqrt(6) * (rho**2) * np.cos(2 * theta)
            elif j == 7:
                z = np.sqrt(8) * (3 * rho**3 - 2 * rho) * np.sin(theta)
            elif j == 8:
                z = np.sqrt(8) * (3 * rho**3 - 2 * rho) * np.cos(theta)
            elif j == 11:
                z = np.sqrt(5) * (6 * rho**4 - 6 * rho**2 + 1)
            else:
                z = np.zeros_like(rho)

            phase += amp * z

        phase[~pupil_mask] = 0.0
        return phase
