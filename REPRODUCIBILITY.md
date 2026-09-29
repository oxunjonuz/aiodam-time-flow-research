# Release packaging and environment notes

This release preserves the completed research scripts, reports, measurement
artifacts and PDFs. It omits the local Python virtual environment, Python caches,
compiled Lean `.olean` files and full copies of external papers/video transcripts.
No research claims or numerical artifacts were rewritten during packaging.

## Python

See `publication/README.md` for the analysis commands and their recorded outcomes.
The original environment is Linux / Python 3.12 inside the AIODAM container.
Scripts require NumPy, SciPy, SymPy and h5py; waveform and sky-map steps also use
phenomxpy, healpy and astropy-healpix. These dependencies are not bundled.

Some scripts use the original absolute project path
`/work/time_flow_research_20260927`. For an unchanged reproduction, mount this
repository at that path in a Linux environment. Some steps refer to
`work/env/venv`; recreate an environment there when using those steps.

## Lean

The original verification scripts expect Lean at `/opt/lean4/bin/lean` and
Mathlib at `/work/Shopify/audit-work/mathlib`. These are environment assumptions,
not portable installation instructions. The `.lean` sources are included;
Lean and the Mathlib cache must be installed separately. Do not assume a latest
Mathlib release is compatible with the original environment.

Recorded from the original container on 2026-09-29:

- Lean **4.19.0**, toolchain `leanprover/lean4:v4.19.0`;
- Mathlib Git commit **`c44e0c8ee63ca166450922a373c7409c5d26b00b`**.

## Verification boundary

For this publication upload, the supplied PDF and publication-consistency
checkers were run successfully (77 PDF checks and 212 consistency checks).
These checks compare the package with its recorded artifacts. They are not an
independent rerun of the full physics analysis or an external peer review.

`publication/MANIFEST.md` covers the publication subdirectory. Root
`SHA256SUMS` covers the complete release, excluding Git metadata and itself.
Checksums establish file integrity, not the truth of scientific claims.
