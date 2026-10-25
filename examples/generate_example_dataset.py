"""Script to generate realistic example data and visualization for Tablasul with wide pupil spacing and astronomical contrast."""

from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

from tablasul import (
    OpticalSystem,
    CircularAperture,
    SubaruAperture,
    PyramidWFS,
    LyotCoronagraph,
    PointSource,
    CompanionSource,
    KolmogorovScreen,
    FrozenFlowAtmosphere,
    ZernikeNCPA,
    DetectorModel,
    SynchronizedCubeEngine,
    export_fits,
    export_npz,
)

def run():
    out_dir = Path("/Users/seamus/Research/tablasul/examples/output")
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Optical Train Configuration
    # Wide 2x2 pupil grid allocation on 256x256 detector
    grid_size = 256
    samples_pupil = 36  # Pupil diameter in pixels

    # On-axis star + 2x10^-3 off-axis companion at 6.0 lambda/D (45 deg)
    star = PointSource(total_flux=1e7)
    planet = CompanionSource(contrast=2e-3, separation_lod=6.0, position_angle_deg=45.0)

    pupil = CircularAperture(
        diameter=8.2,
        grid_size=grid_size,
        samples_across_pupil=samples_pupil,
    )

    # Wide separation ratio (2.2 pupil diameters)
    pywfs = PyramidWFS(
        separation_ratio=2.2,
        modulation_radius=2.0,
        n_modulation_points=16,
    )

    # Lyot Coronagraph with smooth apodized FPM edge
    coronagraph = LyotCoronagraph(
        fpm_radius_lod=2.5,
        lyot_stop_ratio=0.88,
        edge_taper_lod=0.5,
        enabled=True,
    )

    ncpa = ZernikeNCPA(coefficients={5: 0.10, 6: 0.15}, enabled=True)
    detector = DetectorModel(read_noise=1.0, photon_noise=True, bias=0.0, enabled=False)

    system = OpticalSystem(
        aperture=pupil,
        pywfs=pywfs,
        coronagraph=coronagraph,
        sources=[star, planet],
        ncpa=ncpa,
        science_detector=detector,
        wfs_detector=detector,
        wfs_wavelength=800e-9,
        science_wavelength=1650e-9,
    )

    # 2. Turbulence Screen
    atmos_gen = KolmogorovScreen(r0=0.18, L0=25.0)
    phase_screen = atmos_gen.generate_screen(
        grid_size=grid_size,
        pixel_scale=pupil.pixel_scale,
        target_wavelength=800e-9,
        seed=123,
    )

    # 3. Simulate Frame
    res = system.simulate_frame(phase_aberration=phase_screen)

    # 4. Generate Synchronized Cubes
    atmos = FrozenFlowAtmosphere(screen_generator=atmos_gen, wind_velocity=(10.0, 2.0))
    engine = SynchronizedCubeEngine(
        system=system,
        atmosphere=atmos,
        wfs_frame_rate=1000.0,
        science_exposure_time=0.010,
        wfs_output_mode="quadrants",
    )
    dataset = engine.generate_cubes(n_science_frames=10, seed=42)

    # Export cubes
    fits_path = out_dir / "scexao_example_cubes.fits"
    npz_path = out_dir / "scexao_example_cubes.npz"
    export_fits(dataset, fits_path)
    export_npz(dataset, npz_path)
    print(f"Exported FITS to: {fits_path}")
    print(f"Exported NPZ to:  {npz_path}")

    # 5. Publication-grade Diagnostic Visualization
    fig, axes = plt.subplots(1, 4, figsize=(22, 5.5), facecolor="#0e1117")

    # Plot 1: Entrance Pupil Phase Screen (cropped to aperture)
    mid = grid_size // 2
    r_box = int(samples_pupil * 0.9)
    crop_pupil = phase_screen[mid - r_box : mid + r_box, mid - r_box : mid + r_box]
    im0 = axes[0].imshow(crop_pupil, cmap="RdBu_r", origin="lower")
    axes[0].set_title("Entrance Pupil Phase $\\phi(x, y)$ [rad]", color="white", fontsize=12, pad=10)
    axes[0].tick_params(colors="white")
    plt.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04).ax.tick_params(colors="white")

    # Plot 2: PyWFS Detector (Wide 2x2 Grid of 4 Pupils)
    im1 = axes[1].imshow(res.pywfs_detector, cmap="inferno", origin="lower")
    axes[1].set_title("PyWFS Detector (Wide 2x2 Pupils)", color="white", fontsize=12, pad=10)
    axes[1].tick_params(colors="white")
    plt.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04).ax.tick_params(colors="white")

    # Science Focal Plane Units: Normalized to Unocculted Stellar Peak (Contrast Units)
    direct_peak = np.max(res.direct_psf)
    contrast_direct = res.direct_psf / direct_peak
    contrast_corono = res.science_psf / direct_peak

    # Extents in lambda/D
    pixel_scale_lod = samples_pupil / grid_size  # lambda/D per pixel
    crop_radius_lod = 12.0
    crop_pixels = int(round(crop_radius_lod / pixel_scale_lod))
    extent = [-crop_radius_lod, crop_radius_lod, -crop_radius_lod, crop_radius_lod]

    c_dir_crop = contrast_direct[mid - crop_pixels : mid + crop_pixels, mid - crop_pixels : mid + crop_pixels]
    c_cor_crop = contrast_corono[mid - crop_pixels : mid + crop_pixels, mid - crop_pixels : mid + crop_pixels]

    # Plot 3: Direct Science PSF (Contrast Units)
    im2 = axes[2].imshow(
        c_dir_crop,
        cmap="magma",
        origin="lower",
        extent=extent,
        norm=LogNorm(vmin=1e-5, vmax=1.0),
    )
    axes[2].set_title("Direct Science PSF (Contrast)", color="white", fontsize=12, pad=10)
    axes[2].set_xlabel("$\\Delta x$ [$\\lambda/D$]", color="white")
    axes[2].set_ylabel("$\\Delta y$ [$\\lambda/D$]", color="white")
    axes[2].tick_params(colors="white")
    cb2 = plt.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04)
    cb2.ax.tick_params(colors="white")
    cb2.set_label("Normalized Contrast", color="white")

    # Plot 4: Coronagraphic Science PSF (Contrast Units)
    im3 = axes[3].imshow(
        c_cor_crop,
        cmap="magma",
        origin="lower",
        extent=extent,
        norm=LogNorm(vmin=1e-5, vmax=1e-2),
    )
    axes[3].set_title("Coronagraphic PSF (Companion at 6 $\\lambda/D$)", color="white", fontsize=12, pad=10)
    axes[3].set_xlabel("$\\Delta x$ [$\\lambda/D$]", color="white")
    axes[3].set_ylabel("$\\Delta y$ [$\\lambda/D$]", color="white")
    axes[3].tick_params(colors="white")
    cb3 = plt.colorbar(im3, ax=axes[3], fraction=0.046, pad=0.04)
    cb3.ax.tick_params(colors="white")
    cb3.set_label("Normalized Contrast", color="white")

    # Circle indicating FPM radius (2.5 lambda/D)
    fpm_circle = plt.Circle((0, 0), 2.5, color="cyan", fill=False, linestyle="--", linewidth=1.2, label="FPM (2.5 $\\lambda/D$)")
    axes[3].add_patch(fpm_circle)
    axes[3].legend(loc="upper right", facecolor="#161b22", edgecolor="#30363d", labelcolor="white")

    for ax in axes:
        for spine in ax.spines.values():
            spine.set_color("#30363d")

    plt.tight_layout()
    img_path = out_dir / "scexao_simulation_example.png"
    plt.savefig(img_path, dpi=200, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close()
    print(f"Saved diagnostic figure to: {img_path}")

    # Copy to artifacts directory
    artifact_dir = Path("/Users/seamus/.gemini/antigravity/brain/bd82fcab-1201-4c84-a32e-014e18a799e1")
    artifact_dir.mkdir(parents=True, exist_ok=True)
    import shutil
    shutil.copy(img_path, artifact_dir / "scexao_simulation_example.png")
    print("Copied figure to artifacts directory.")

if __name__ == "__main__":
    run()
