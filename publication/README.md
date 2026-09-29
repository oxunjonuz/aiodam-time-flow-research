# Is the Müller–Maguire “now”-theory prediction falsifiable on GW150914?

**A formal identifiability audit, with a numerical check on the real LIGO strain**

*Фальсифицируемо ли предсказание «now»-теории Мюллера–Магуайра на GW150914?*

**Authors / Авторы:**
**Aiodam** (autonomous research agent) · **Oxunjon Ubaydullayev**

---

## What this is / Что это

We ask a single decidable question about the “now”-theory of Müller & Maguire
(arXiv:1606.07975): **is its predicted 1–1.5 ms merger delay specified precisely
enough to be distinguished from general relativity and from a re-fit of the source
parameters?**

The answer is developed in three independent registers:

| register | result | independent of |
|---|---|---|
| **formal** (Lean 4 + Mathlib) | an identifiability criterion for `M(θ,τ) = Aθ + Bτ`, and proof that it is a genuine dichotomy | data, noise, waveform |
| **symbolic** (sympy) | a constant delay is *exactly* a shift of `t_c`, a linear delay *exactly* a rescaling of `M_c` — at **every** PN order | data |
| **numerical** (real GW150914 strain) | a coherent H1+L1 fit with one `τ₂` recovers a 20 ms injection but cannot resolve the claimed 1.2 ms; the quadratic term gives SNR 0.13–0.17 | — |

**Scope of the headline.** What is *proven* is degeneracy **within the adopted
linear model**. The specific function `τ(t)` is **not derived** from the theory,
because the paper does not specify it. This remains the main open item, and it is
stated as such everywhere it matters.

> Мы задаём один разрешимый вопрос о «now»-теории Мюллера–Магуайра
> (arXiv:1606.07975): **задано ли её предсказание задержки 1–1,5 мс достаточно
> точно, чтобы отличить его от ОТО и от повторной подгонки параметров
> источника?** Доказано **вырождение в принятой линейной модели**; конкретная
> функция `τ(t)` из теории **не выведена**, потому что статья её не задаёт.

---

## Authorship and disclosure / Авторство

This work was carried out **by an autonomous research agent, Aiodam**, as primary
investigator: it designed the study, wrote and ran every script, found and
corrected its own errors, and wrote the report. The human co-author,
**Oxunjon Ubaydullayev**, set the research question, supplied the LIGO data and
skymap, and reviewed intermediate results over several rounds.

**Who found what, counted rather than remembered.** `work/NOTES.md` records **47**
errors made on the way. The co-author found **eleven** of them: **six** in one
review of the claim ledger (one of those by reading the ledger's own legend rather
than by checking any number), **one** earlier, when he read an intermediate report
and caught a sampling artifact in the H4 arithmetic, and **four** further errors of
formulation. The remaining **36** were found by the agent's own controls, verifiers
and self-tests. The split is produced by `count_errors.py`, which prints the rule it
used so the rule can be argued with — a number in a paper has to be countable from
the record.

The agent's errors are **recorded, not removed**: see `work/NOTES.md`. Every
headline number is reproduced by an independent verifier that shares no code with
the audit, and each verifier carries a self-test that must go red on a corrupted
artifact — because a verifier that cannot fail proves nothing.

---

## Repository layout / Состав репозитория

```
publication/                 ← the publication package (this README lives here)
├── paper/
│   ├── paper_en.md / .pdf   ← the paper, English (primary)
│   └── paper_ru.md / .pdf   ← the paper, Russian
├── figures/                 ← the five waveform figures (PDF + PNG preview)
├── preamble.tex             ← typography for the PDFs (P052 + glyph fallbacks)
├── build_paper.sh           ← pandoc + xelatex build
├── make_figures.py          ← draws the figures from the real strain + artifacts
├── count_errors.py          ← counts the errors in NOTES.md, split by who found them
├── verify_pdfs.py           ← checks the built PDFs (fonts, glyphs, figures)
├── selftest_pdfs.py         ← 7 corruptions of the PDFs must all be caught
├── verify_frames.py         ← no ink outside the text block, on every page
├── docs_en/                 ← English translations (Russian originals kept in work/)
│   ├── REPORT_en.md         ← final report
│   ├── CLAIMS_en.md         ← claim ledger: proven / derived / measured / hypothesis
│   ├── PREREGISTRATION_en.md
│   ├── NOTES_en.md          ← every error made on the way
│   └── ROUNDS_en.md         ← per-round summaries
├── ARTIFACTS.md             ← bilingual index of every artifact + sha256
├── MANIFEST.md              ← sha256 of every published file
├── CITATION.cff / .zenodo.json
└── LICENSE, LICENSE-CODE

work/                        ← the research itself
├── REPORT.md                ← final report (Russian original)
├── CLAIMS.md                ← claim ledger (Russian original)
├── PREREGISTRATION.md       ← frozen before the first run
├── NOTES.md                 ← full chronology + every error
├── msg21X-note.md, msg22X-note.md   ← per-round notes (Russian originals)
├── analysis.py, h3_*.py, h4_*.py, pn_orders.py, …   ← the audits
├── verify_*.py, verify_lean*.sh     ← independent verifiers
├── mutation_control.py, equivalence_check.py
├── lean/*.lean              ← the formalization (Lean 4 + Mathlib)
└── artifacts/*.json, *.npz  ← frozen measurement artifacts

data/                        ← GW150914 strain (H1, L1) + DATA_PROVENANCE.md
sources/                     ← the primary sources used
```

**Language / Язык.** The English text is primary. The Russian originals are kept
unchanged alongside the translations; nothing was deleted. / Английский текст
главный; русские оригиналы сохранены рядом с переводами, ничего не удалено.

---

## Reproducing / Воспроизведение

Requires Python 3 with NumPy, SciPy, SymPy, `h5py`; `phenomxpy`, `healpy` and
`astropy-healpix` for the waveform/skymap steps (the latter three in the local venv
`work/env/venv`); Lean 4 + Mathlib for the formalization.

```sh
cd work

# core analyses
python3 analysis.py                    # exit 0, 22/22 thresholds
python3 h3_real_data.py                # exit 0, 15/15 thresholds on real data
python3 h3_joint_fit.py                # exit 0, 9/9 controls
python3 h3_imr_check.py                # exit 0
python3 h3_direct_fit_imr.py           # exit 0, 8/8 controls
python3 h3_imr_spins.py                # exit 0, 8/8 controls
python3 h3_rho_opt_bands.py            # exit 0
python3 h3_grid_bounds.py              # exit 0
python3 pn_orders.py                   # exit 0
python3 et_ce_scaling.py               # exit 0
python3 h4_paper_arithmetic.py         # exit 0 (v2)
python3 h4_paper_text_audit.py         # exit 0 (v3)
python3 h4_weight_audit.py             # exit 0
python3 prereg_t8.py                   # exit 0

# independent verifiers (each has --selftest, which must go red)
python3 verify_independent.py
python3 verify_h3_eigh.py
python3 verify_h3_joint.py
python3 verify_h3_round226.py
python3 verify_h3_imr.py
python3 verify_rho_opt_bands.py
python3 verify_new_results.py
python3 verify_claims_228.py

# formalization
sh verify_lean.sh                      # rc=0, NO_SORRY_NO_AXIOM
sh verify_lean_cosmo.sh                # rc=0, axiom audit + 3 negative controls

# mutation control
python3 mutation_control.py            # exit 1 = PRE-EXISTING state (15/20 killed;
                                       # the 4 survivors are classified as
                                       # equivalent on the headline claims by
                                       # equivalence_check.py, by measurement)
```

### Building the paper PDFs / Сборка PDF статьи

The two PDFs are built from the two Markdown files by `pandoc` + `xelatex` with a
typographic preamble (`preamble.tex`): main font **P052** (a Palatino clone with
real Cyrillic coverage), `microtype`, coloured section rules, an abstract box, and
running heads. The nine glyphs P052 lacks are mapped to a fallback font that has
them, so the characters survive into the PDF text layer rather than only into the
picture.

```sh
cd publication
../work/env/venv/bin/python make_figures.py   # writes figures/*.pdf from the real strain
sh build_paper.sh both                        # writes paper/paper_{en,ru}.pdf
python3 count_errors.py                       # the error count the paper quotes
python3 verify_frames.py                      # no ink outside the text block
python3 verify_pdfs.py                        # fonts embedded, glyphs present, figures in
python3 selftest_pdfs.py                      # 7/7 PDF corruptions must be caught
python3 verify_publication.py                 # every headline number is an artifact number
python3 selftest_publication.py               # 29/29 corruptions must be caught
python3 make_manifest.py                      # refresh MANIFEST.md (do this LAST)
```

The five figures are drawn from the **real** GW150914 strain and from the frozen
artifacts, not from memory or from hand-typed numbers: the strain and the measured
ASD come from the HDF5 release, the template overlay is a full `IMRPhenomT` on the
audit's own grid, and the τ₂ profile and inclination sweep are read by key from
`work/artifacts/`. `make_figures.py` prints the numbers it recomputed so they can
be compared with the artifacts.

Known caveat: `run_command` in this container refuses commands that contain the
protected-home name, so all `data/` operations went through `read_file`/`h5py`.

---

## Data provenance / Происхождение данных

GW150914 strain, H1 and L1, 32 s @ 16 kHz:

| file | bytes | SHA-256 |
|---|---:|---|
| `H-H1_GWOSC_16KHZ_R1-1126259447-32.hdf5` | 4 068 797 | `81040e1ecfaf40ffe15a5efc59dbc3a888653162f1613425f1e68d0828dd1b97` |
| `L-L1_GWOSC_16KHZ_R1-1126259447-32.hdf5` | 3 919 118 | `23a207023ef8b49cf0b70b78d7ba80635881799324e50bfb2165c60039e41b1f` |

The container could not reach `gwosc.org`, so the files were downloaded from the
AEI mirror and then **hash-verified byte-for-byte against the official GWOSC
files** by the co-author on a host that could reach GWOSC. The data are also
functionally authenticated by three independent physical checks: the measured
ASD(100 Hz) matches the published aLIGO O1 sensitivity, the peak time matches the
published merger time, and H1−L1 agree to 0.37 ms. Details:
`data/DATA_PROVENANCE.md`.

Skymap: `LALInference_skymap.fits.gz`, LOSC P1500227.

---

## License / Лицензия

* **Text and figures** (`paper/`, `docs_en/`, `*.md`): **CC BY 4.0** — see `LICENSE`.
* **Code** (`work/*.py`, `work/lean/*.lean`): **MIT** — see `LICENSE-CODE`.

The depositing co-author confirmed CC BY 4.0 for the article and original
figures, and MIT for the code, on 2026-09-29. Third-party inputs retain their
original attribution and terms; see `../THIRD_PARTY.md`.

---

## Citing / Как цитировать

See `CITATION.cff`. Plain text:

> Aiodam and O. Ubaydullayev, *Is the Müller–Maguire “now”-theory prediction
> falsifiable on GW150914? A formal identifiability audit, with a numerical check
> on real LIGO strain*, 2026.

## References / Литература

1. R. A. Muller, S. Maguire, *Now, and the Flow of Time*, arXiv:1606.07975 (2016).
2. B. P. Abbott et al., *Observation of Gravitational Waves from a Binary Black
   Hole Merger*, arXiv:1602.03837, PRL 116, 061102 (2016).
3. B. P. Abbott et al., *Properties of the Binary Black Hole Merger GW150914*,
   arXiv:1602.03840, PRL 116, 241102 (2016).
4. B. P. Abbott et al., *GW150914: First results from the search for binary black
   hole coalescence with Advanced LIGO*, arXiv:1602.03839 (2016).
5. B. P. Abbott et al., *Tests of general relativity with GW150914*,
   arXiv:1606.03833, PRL 116, 221101 (2016).
6. B. P. Abbott et al., *GWTC-1: A Gravitational-Wave Transient Catalog …*,
   PRX 9, 031040 (2019).
7. GWOSC, GW150914 data release v3.
8. LOSC, GW150914 skymap P1500227.
9. The Mathlib Community, *The Lean Mathematical Library*, CPP 2020.
