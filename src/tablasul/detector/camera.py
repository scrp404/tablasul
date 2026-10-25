"""Science and WFS detector noise models."""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass
class DetectorModel:
    """Scientific detector noise model (EMCCD / sCMOS / InGaAs).

    Attributes:
        read_noise: Readout noise standard deviation in electrons (RMS).
        dark_current: Dark current rate in electrons/second/pixel.
        photon_noise: Whether to simulate Poisson photon shot noise.
        gain: Detector gain in electrons per Digital Number (e-/ADU).
        bias: Constant electronic bias offset added to output in ADU.
        enabled: Boolean toggle. When False, outputs pure floating point intensities.
    """

    read_noise: float = 1.0
    dark_current: float = 0.01
    photon_noise: bool = True
    gain: float = 1.0
    bias: float = 100.0
    enabled: bool = True

    def apply(
        self,
        intensity: np.ndarray,
        exposure_time: float = 1.0,
        seed: int | None = None,
    ) -> np.ndarray:
        """Apply detector noise to an incident photon flux array."""
        if not self.enabled:
            return intensity.copy()

        rng = np.random.default_rng(seed)

        # 1. Expected electron signal = incident photons + dark current
        expected_electrons = np.maximum(intensity + (self.dark_current * exposure_time), 0.0)

        # 2. Poisson photon shot noise
        if self.photon_noise:
            electrons = rng.poisson(expected_electrons).astype(np.float64)
        else:
            electrons = expected_electrons.copy()

        # 3. Gaussian readout noise
        if self.read_noise > 0:
            ron = rng.normal(loc=0.0, scale=self.read_noise, size=intensity.shape)
            electrons += ron

        # 4. Conversion to ADU (counts) and bias addition
        counts = (electrons / self.gain) + self.bias
        return np.maximum(counts, 0.0)
