# Tasks: Refine PyWFS Pupil Separation, Coronagraph Contrast, and Documentation

**Branch**: `002-optics-refinement` | **Spec**: [spec.md](file:///Users/seamus/Research/tablasul/specs/002-optics-refinement/spec.md) | **Plan**: [plan.md](file:///Users/seamus/Research/tablasul/specs/002-optics-refinement/plan.md)

---

## Phase 1: PyWFS Separation Refinement (User Story 1)

- [x] T032 Update `PyramidWFS` default `separation_ratio` to `2.0` in `src/tablasul/optics/pywfs.py` and ensure quadrant extraction window radius scales appropriately
- [x] T033 Verify quadrant separation and slope response tests pass with wide separation in `tests/test_pywfs.py`

---

## Phase 2: Coronagraph Contrast & Smooth FPM Refinement (User Story 2)

- [x] T034 Add optional `edge_taper_lod` to `LyotCoronagraph` in `src/tablasul/optics/coronagraph.py`
- [x] T035 Update `examples/generate_example_dataset.py` with wide pupil separation, smooth FPM, and true normalized astronomical contrast plots ($10^{-5}$ to $1$)

---

## Phase 3: Git Repository & Comprehensive Documentation (User Story 3)

- [x] T036 Create `.gitignore` in repository root
- [x] T037 Initialize Git repository with `git init` and create initial commit
- [x] T038 Create `docs/optical_theory.md` with complete mathematical and wave-optics derivations
- [x] T039 Create `docs/psf_reconstruction_guide.md` with guide on training ML models with the synchronized cubes
- [x] T040 Update `README.md` with full documentation index, architecture diagram, and GitHub sharing instructions

---

## Phase 4: Decoupled MFT Science Camera & Measured Pupil Support

- [x] T041 Implement Matrix Fourier Transform (`propagate_to_focal_mft`) in `src/tablasul/core/field.py` for decoupled science FOV and detector pixel scaling
- [x] T042 Implement `PupilImageAperture` in `src/tablasul/core/aperture.py` supporting custom measured Subaru pupil files (FITS, NPY, PNG)
- [x] T043 Add unit tests in `tests/test_aperture.py` covering array, FITS, NPY loading, and end-to-end optical system execution
