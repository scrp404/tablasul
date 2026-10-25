"""Classical Lyot Coronagraph optical subsystem."""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ..core.field import ComplexElectricField


@dataclass
class LyotCoronagraph:
    """Classical Lyot Coronagraph with opaque focal plane mask and undersized Lyot stop.

    Attributes:
        fpm_radius_lod: Focal plane mask radius in units of lambda / D (default: 2.5).
        lyot_stop_ratio: Undersizing factor for the Lyot stop relative to pupil diameter (default: 0.90).
        edge_taper_lod: Smooth taper width at FPM edge in lambda/D to suppress Gibbs ringing (0 for hard-edge).
        enabled: Boolean toggle. When False, light bypasses the coronagraph completely.
    """

    fpm_radius_lod: float = 2.5
    lyot_stop_ratio: float = 0.90
    edge_taper_lod: float = 0.4
    enabled: bool = True

    def propagate(
        self,
        focal_field: ComplexElectricField,
        samples_across_pupil: int,
    ) -> ComplexElectricField:
        """Propagate a focal-plane field through the coronagraph.

        If enabled=False, returns the input field unchanged.
        """
        if not self.enabled:
            return focal_field

        n = focal_field.grid_size
        mid = n // 2
        y, x = np.ogrid[:n, :n]

        # 1. Focal Plane Mask (FPM)
        pixel_scale_lod = samples_across_pupil / n
        r_lod = np.hypot(x - mid, y - mid) * pixel_scale_lod

        if self.edge_taper_lod <= 0.0:
            # Hard-edge opaque mask
            fpm_mask = (r_lod >= self.fpm_radius_lod).astype(np.float64)
        else:
            # Smooth cosine apodized edge transition
            r_in = max(0.0, self.fpm_radius_lod - self.edge_taper_lod / 2.0)
            r_out = self.fpm_radius_lod + self.edge_taper_lod / 2.0
            taper_range = r_out - r_in

            fpm_mask = np.zeros_like(r_lod)
            # Fully transmitting outside r_out
            fpm_mask[r_lod >= r_out] = 1.0
            # Smooth transition in [r_in, r_out]
            transition = (r_lod >= r_in) & (r_lod < r_out)
            fpm_mask[transition] = 0.5 * (1.0 - np.cos(np.pi * (r_lod[transition] - r_in) / taper_range))

        masked_focal = focal_field.apply_mask(fpm_mask)

        # 2. Re-image to Lyot pupil plane
        lyot_field = masked_focal.propagate_to_pupil()

        # 3. Apply undersized Lyot stop
        pupil_radius_pixels = samples_across_pupil / 2.0
        lyot_radius_pixels = pupil_radius_pixels * self.lyot_stop_ratio
        r_pupil = np.hypot(x - mid, y - mid)
        lyot_stop_mask = (r_pupil <= lyot_radius_pixels).astype(np.float64)

        masked_lyot = lyot_field.apply_mask(lyot_stop_mask)

        # 4. Propagate to final science focal plane
        science_focal = masked_lyot.propagate_to_focal()
        return science_focal
