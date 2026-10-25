# Feature Specification: Refine PyWFS Pupil Separation, Coronagraph Contrast, and Documentation

**Feature Branch**: `002-optics-refinement`  
**Created**: 2026-10-25  
**Status**: Draft  
**Project**: Tablasul (`/Users/seamus/Research/tablasul`)  
**Input**: Increase PyWFS 4-pupil spatial separation to match real lab detectors, refine coronagraphic diffraction physics and contrast-normalized visualization, and package into a shareable Git repository with comprehensive documentation.

---

## User Scenarios & Testing

### User Story 1 - Configurable Wide PyWFS Pupil Separation (Priority: P1)

As an adaptive optics researcher viewing PyWFS detector frames,  
I want the 4 re-imaged circular pupil footprints to have wide, realistic separation on the detector ($2.0\times$ pupil diameter or larger),  
So that the 4 quadrants resemble real-world camera layouts with clear dark spacing between pupils and no crowding.

**Acceptance Scenarios**:
1. **Given** a PyWFS with default `separation_ratio = 2.0`, **When** simulated on a circular or Subaru aperture, **Then** the distance between pupil centers is at least twice the pupil diameter, leaving a clear gap $> 1.0\times D_{\text{pup}}$ between the inner edges of adjacent pupils.
2. **Given** quadrant extraction logic, **When** `separation_ratio` is modified by the user (e.g. from 1.5 to 2.5), **Then** the quadrant extraction windows dynamically track the re-imaged pupil centers without cropping pupil edges.

---

### User Story 2 - Coronagraph High-Contrast Physics & Display Refinement (Priority: P2)

As a high-contrast imaging scientist,  
I want the coronagraphic PSF to accurately reflect astronomical contrast ratios ($10^{-6}$ to $10^{-1}$ relative to the unocculted stellar peak) with optional smooth-edge mask apodization,  
So that the on-axis star suppression and off-axis companion detection are displayed in standard normalized contrast units rather than arbitrary counts.

**Acceptance Scenarios**:
1. **Given** a Lyot coronagraph with an on-axis star, **When** normalized to the peak intensity of the unocculted direct PSF, **Then** the coronagraphic core exhibits a contrast depth of $< 10^{-3}$ while preserving an off-axis companion at its true contrast (e.g. $10^{-3}$ or $10^{-4}$).
2. **Given** the option for smooth/tapered mask transitions (`fpm_type='hard'` or `'smooth'`), **When** smooth is selected, **Then** edge-diffraction ringing at the mask boundary is suppressed.

---

### User Story 3 - Repository Packaging & Documentation Suite (Priority: P3)

As a researcher sharing Tablasul with colleagues and collaborators,  
I want a clean Git repository with comprehensive documentation (optical theory, math, usage guides, API reference, and ready-to-run scripts),  
So that any team member can clone, install via `uv`, and immediately understand the physics and code.

**Acceptance Scenarios**:
1. **Given** the repository root, **When** `git status` or `git log` is run, **Then** the project is under version control with a clean `.gitignore` excluding caches and virtual environments.
2. **Given** documentation in `README.md` and `docs/`, **When** inspected, **Then** the mathematical wave-optics principles, PyWFS geometry, coronagraph design, and synchronization engine are clearly documented with code examples.
