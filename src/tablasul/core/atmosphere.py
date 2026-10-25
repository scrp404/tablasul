"""Atmospheric turbulence generation and frozen-flow dynamics."""

from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np


@dataclass
class KolmogorovScreen:
    """Von Kármán / Kolmogorov atmospheric turbulence phase screen generator.

    Attributes:
        r0: Fried parameter at 500 nm in meters (e.g. 0.15 m for good seeing).
        L0: Outer scale in meters (default: 25.0 m).
        reference_wavelength: Reference wavelength for r0 in meters (default: 500 nm).
    """

    r0: float = 0.15
    L0: float = 25.0
    reference_wavelength: float = 500e-9

    def r0_at(self, wavelength: float) -> float:
        """Scale Fried parameter to target wavelength (r0 ~ lambda^(6/5))."""
        return self.r0 * (wavelength / self.reference_wavelength) ** 1.2

    def generate_screen(
        self,
        grid_size: int,
        pixel_scale: float,
        target_wavelength: float,
        seed: int | None = None,
    ) -> np.ndarray:
        """Generate a 2D phase screen (in radians at target_wavelength)."""
        rng = np.random.default_rng(seed)
        n = grid_size
        mid = n // 2
        r0_target = self.r0_at(target_wavelength)

        # Spatial frequencies in cycles per meter
        fx = (np.arange(n) - mid) / (n * pixel_scale)
        fy = (np.arange(n) - mid) / (n * pixel_scale)
        FX, FY = np.meshgrid(fx, fy)
        k = np.hypot(FX, FY)

        # Modified von Kármán power spectral density
        k0 = 1.0 / self.L0
        psd = 0.0229 * (r0_target ** (-5.0 / 3.0)) * ((k**2 + k0**2) ** (-11.0 / 6.0))
        psd[mid, mid] = 0.0

        # Complex Gaussian random field
        white_noise = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2.0)
        field_fft = np.fft.ifftshift(white_noise * np.sqrt(psd) * (1.0 / pixel_scale))
        phase_screen = np.real(np.fft.ifft2(field_fft)) * n
        return phase_screen


@dataclass
class FrozenFlowAtmosphere:
    """Dynamic atmosphere advancing according to Taylor's frozen-flow hypothesis.

    Uses continuous sub-pixel Fourier phase-shifting for exact, artifact-free translation.

    Attributes:
        screen_generator: KolmogorovScreen generator.
        wind_velocity: Wind vector (v_x, v_y) in meters/second.
    """

    screen_generator: KolmogorovScreen = field(default_factory=KolmogorovScreen)
    wind_velocity: tuple[float, float] = (10.0, 0.0)
    _base_screen_fft: np.ndarray | None = field(default=None, init=False, repr=False)
    _pixel_scale: float = field(default=8.2 / 64, init=False, repr=False)
    _target_wavelength: float = field(default=800e-9, init=False, repr=False)
    _grid_size: int = field(default=64, init=False, repr=False)

    def initialize(
        self,
        pupil_grid_size: int,
        pixel_scale: float,
        target_wavelength: float,
        seed: int | None = None,
    ) -> None:
        """Generate the base phase screen and store its Fourier representation."""
        self._grid_size = pupil_grid_size
        self._pixel_scale = pixel_scale
        self._target_wavelength = target_wavelength

        screen = self.screen_generator.generate_screen(
            grid_size=pupil_grid_size,
            pixel_scale=pixel_scale,
            target_wavelength=target_wavelength,
            seed=seed,
        )
        self._base_screen_fft = np.fft.fft2(screen)

    def get_phase_screen(
        self,
        time_seconds: float,
        pupil_grid_size: int,
    ) -> np.ndarray:
        """Extract pupil phase screen at time_seconds via continuous Fourier phase shift."""
        if self._base_screen_fft is None or self._grid_size != pupil_grid_size:
            self.initialize(pupil_grid_size, self._pixel_scale, self._target_wavelength)

        vx, vy = self.wind_velocity
        dx = vx * time_seconds
        dy = vy * time_seconds

        n = pupil_grid_size
        # Frequencies for standard unshifted 2D FFT
        fx = np.fft.fftfreq(n, d=self._pixel_scale)
        fy = np.fft.fftfreq(n, d=self._pixel_scale)
        FX, FY = np.meshgrid(fx, fy)

        # Fourier shift theorem: f(x - dx, y - dy) <=> F(fx, fy) * exp(-2*pi*i*(fx*dx + fy*dy))
        shift_phase = -2.0 * np.pi * (FX * dx + FY * dy)
        shifted_fft = self._base_screen_fft * np.exp(1j * shift_phase)
        return np.real(np.fft.ifft2(shifted_fft))
