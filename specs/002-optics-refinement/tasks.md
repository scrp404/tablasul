# Tasks: Refine PyWFS Pupil Separation, Coronagraph Contrast, and Documentation

**Branch**: `002-optics-refinement` | **Spec**: [spec.md](file:///Users/seamus/Research/tablasul/specs/002-optics-refinement/spec.md) | **Plan**: [plan.md](file:///Users/seamus/Research/tablasul/specs/002-optics-refinement/plan.md)

---

## Phase 1: PyWFS Separation Refinement (User Story 1)

- [ ] T032 Update `PyramidWFS` default `separation_ratio` to `2.0` in `src/tablasul/optics/pywfs.py` and ensure quadrant extraction window radius scales appropriately
- [ ] T033 Verify quadrant separation and slope response tests pass with wide separation in `tests/test_pywfs.py`

---

## Phase 2: Coronagraph Contrast & Smooth FPM Refinement (User Story 2)

- [ ] T034 Add optional `edge_taper_lod` to `LyotCoronagraph` in `src/tablasul/optics/coronagraph.py`
- [ ] T035 Update `examples/generate_example_dataset.py` with wide pupil separation, smooth FPM, and true normalized astronomical contrast plots ($10^{-5}$ to $1$)

---

## Phase 3: Git Repository & Comprehensive Documentation (User Story 3)

- [ ] T036 Create `.gitignore` in repository root
- [ ] T037 Initialize Git repository with `git init` and create initial commit
- [ ] T038 Create `docs/optical_theory.md` with complete mathematical and wave-optics derivations
- [ ] T039 Create `docs/psf_reconstruction_guide.md` with guide on training ML models with the synchronized cubes
- [ ] T040 Update `README.md` with full documentation index, architecture diagram, and GitHub sharing instructions
