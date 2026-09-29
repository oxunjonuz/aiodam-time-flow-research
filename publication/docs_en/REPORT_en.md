# Final report — cycle 2026-09-27

**Topic:** the thermodynamic arrow of time and the “growing-block universe”.
**Vector:** a logical-mathematical audit of the **falsifiability** of the
prediction of the Müller–Maguire “now”-theory, with a numerical check on a
gravitational-wave signal model.
**Pre-registration:** `work/PREREGISTRATION.md`, frozen before the first run.
**Run date:** 2026-09-27.

> This is the English translation of `work/REPORT.md`. The Russian original is
> kept alongside it; the English text is the primary one for publication.

---

## 1. Verdict in one paragraph

**What is proven formally.** In Lean 4, without `sorry` and without custom
axioms, **A1–A3** are proven: the identifiability criterion `identifiable_iff`
in the additive model `M(θ,τ) = Aθ + Bτ` and its non-vacuity (`Instances.lean`).
This **does not depend** on data, noise or waveform.

**What is derived symbolically (not in Lean).** The identities
`∂Ψ/∂t_c = 2πf` and `f·∂Ψ_chirp/∂f = M_c·∂Ψ_chirp/∂M_c` hold at **every**
post-Newtonian order (Euler’s theorem in `v`; no PN coefficient enters) —
verified with symbolic algebra (sympy) and a negative control, but **there is no
formal proof in Lean**. From them it follows that a constant delay `τ₀` and a
linearly growing `τ̇` are **exactly** equivalent to a shift of the coalescence
time `t_c` and a rescaling of the chirp mass `M_c`, so the degeneracy is not an
artifact of the leading order. Numerically, on GW150914, the fraction of the
delay absorbed by freely fitted source parameters is **1 − 1.4·10⁻¹⁵** for `τ₀`
and **1 − 1.7·10⁻¹⁵** for `τ̇`. The first **non-degenerate** term is quadratic,
suppressed as `(v/c)⁸`; the claimed 1.2 ms gives an SNR of **0.087**. At
**leading order** this corresponds to a detectable amplitude of **69 ms**; that
number is a leading-order result and does not survive a full IMR template (see
below).

**What is measured under named assumptions.** The joint H1+L1 fit with **one**
`τ₂`, free `(M_c, t_c, φ_c)` and antenna responses `F₊, F×` for the published
GW150914 sky position gives a peak network SNR of **22.93** (pn kernel) /
**23.07** (imr kernel) against the published **network** value 24, and on the
±38.4 ms grid a one-sided bound on `τ₂` (drop 3.50σ / 2.35σ). **This is an upper
estimate of the strength of the constraint, not a measured bound:** free `η` and
spins would absorb more of the quadratic term (with `η` fixed and zero spins,
88 % of the term is already absorbed by the source parameters), so the true
constraint is weaker. The same fit **recovers** a 20 ms injection and **does not
resolve** 1.2 ms: **the current fit resolution is ≈ 5 ms** (spread 4.6 ms at an
amplitude of 1.2 ms; for the 20 ms injection the spread is 5.1 ms). The
assumptions are listed: spins 0, phase-only delay, `η = 29/36`, **sky position
fixed at the maximum-likelihood point of the published map** (90 % area
616 deg² — the position is not fitted jointly with `τ₂`).

**What the headline is not allowed to claim.** That the “now”-theory is false,
or that it cannot be tested **in principle**. What is proven is degeneracy
**within the adopted linear model** — and the specific function `τ(t)` is **not
derived** from the theory: the paper does not specify it. This remains the main
open item. The full tagging of every claim under **four** categories (in Lean,
**A1–A3** are proven; A4–A6 are **derived symbolically** but have no Lean file)
is in `work/CLAIMS.md`; no hypothesis-category claim is used in the headline.

This **strengthens** the paper’s own caveat (§5: “the one observed event falls
short of providing a meaningful test”), translating it from prose into a theorem
— and refines it: the issue is not the quality of one event, but that the
direction of the effect lies inside the linear span of the source parameters.

---

## 2. What exactly is proven formally (Lean 4 + Mathlib)

`work/lean/Identifiability.lean` — the identifiability criterion in the additive
model `M(θ,τ) = Aθ + Bτ`:

> **`identifiable_iff`.** `τ` is identifiable (⟺ `M(θ,τ) = M(θ',τ') ⟹ τ = τ'`)
> **if and only if** `B` is injective and `range A ⊓ range B = ⊥`.

> **`not_identifiable_of_range_le`.** If `range B ≤ range A`, then `τ` is not
> identifiable: any effect of `τ` is absorbed by a shift of `θ`.

`work/lean/Instances.lean` — **non-vacuity**. A theorem of indistinguishability
is useless if its premise is unfulfillable. Two concrete examples are exhibited
in the same formal field: `demo_identifiable` (disjoint directions → `τ` is
identifiable) and `demo_not_identifiable` (delay direction inside the image →
not identifiable). The criterion is a genuine dichotomy.

`work/lean/RelativityOfNow.lean` — the relativity-of-simultaneity objection:

> **`exists_boost_reversing_time_order`.** In 1+1 Minkowski, if `|Δt| < |Δx|`
> and `Δt ≠ 0`, there is a boost with `|v| < 1` reversing the sign of `Δt'`.

`work/lean/PastHypothesis.lean` — **variant B**: the arrow of time from boundary
conditions (the past hypothesis). See §11.

**Soundness of the formalization.** All five files compile, `NO_SORRY_NO_AXIOM`,
and the axiom audit reports only the three standard Lean axioms
`[propext, Classical.choice, Quot.sound]`. No custom axioms.
`HomogeneousUniverse.lean` — variant (c): a homogeneous universe, `a(t)` as a
reparametrization. See §11 and `work/NOTES.md` (chronology).

---

## 3. Algebraic core of the result

With `x = (π M_c f)^{2/3}`, `t(f) − t_c = −(5/256) M_c x^{−4}` and
`2πf = 2x^{3/2}/M_c`, the expansion `τ(t) = τ₀ + τ₁(t−t_c) + τ₂(t−t_c)² + …` gives

```
ΔΨ = 2 x^{3/2} τ₀ / M_c                      ∝ ∂Ψ/∂t_c        ← ABSORBED
   − (5/128) τ₁ x^{−5/2}                     ∝ ∂Ψ/∂M_c        ← ABSORBED
   + (25/32768) τ₂ M_c x^{−13/2} + …         ← first non-degenerate
```

The leading phase is `Ψ_chirp = (3/128) x^{−5/2}`, `∂Ψ/∂M_c = −(5/128) x^{−5/2}/M_c`.
The first two terms are **exact** reparameterizations, not approximations. The
third sits on `x^{−13/2}`, four powers of `x^{−4} ~ (v/c)⁸` below the leading
term. Verified symbolically (sympy, an independent path): the coefficients agree
exactly, `τ₁ / (∂Ψ/∂M_c) = +M_c`.

---

## 4. Numerical results (GW150914, published parameters)

Source parameters: `M_c = 28.096 M_⊙`, `D_L = 410 Mpc`, `R_s(62 M_⊙) = 183.148 km`.
Band 20–300 Hz, 4000 bins. The analytic PSD is **calibrated** so that the SNR of
the injected signal equals the published single-detector value **20** (the raw
PSD gave 129 — that was an error of v1–v2, corrected).

| quantity | value | threshold | verdict |
|---|---:|---|---|
| fraction of `τ₀` absorbed by the source | 1 − 1.4·10⁻¹⁵ | ≤ 10⁻¹⁰ | **H1 confirmed** |
| unabsorbed SNR for `τ₀` | 3.3·10⁻¹⁴ | ≤ 10⁻⁶ | **H1 confirmed** |
| fraction of `τ̇` absorbed by the source | 1 − 1.7·10⁻¹⁵ | ≤ 10⁻¹⁰ | **H2 confirmed** |
| unabsorbed SNR for `τ̇` | 2.5·10⁻¹³ | ≤ 10⁻⁶ | **H2 confirmed** |
| SNR of the quadratic term at 1.2 ms | **0.087** | < 1 | **H3 confirmed** |
| detectable amplitude at SNR 5 | **68.8 ms** | — | **leading order**; 57× the claimed value |
| Fisher condition number (H1) | ~1.8·10²⁰ | > 10¹² | **numerically singular** |

**22 of 22 thresholds pass.** Full artifact: `work/artifacts/analysis_results.json`.

### Power control — without it, “did not detect” is indistinguishable from “cannot detect”

| control | result | what it proves |
|---|---|---|
| K1 zero injection | SNR exactly 0 | the procedure does not invent a signal |
| K2 remove `t_c` from the fit | SNR **7.49** | the degeneracy is specifically with `t_c` |
| K3 remove `M_c` from the fit | SNR **18.4** | the degeneracy is specifically with `M_c` |
| K4 random smooth phase | absorbed 0.828 | the procedure **sees** the non-degenerate |
| K5 ×100 amplification | SNR **8.72** | the machinery is not blind |
| K6 sign flip of `τ̇` | absorbed 1.000 | the degeneracy is symmetric |
| K7 a **different** PSD | absorbed 1.000 | the degeneracy is independent of the noise |
| K8 flat likelihood | spread **7·10⁻¹⁵** | the data do not distinguish `τ₀` |
| K9 the same test for `M_c` | minimum **exactly** at the injection | the test discriminates when it can |

**K8 versus K9 is the strongest argument.** The same test: fit `(M_c, t_c, φ_c)`
with a **fixed** candidate. For `τ₀` the residual is identical for candidates
0.0, 0.3, 0.6, 1.2, 2.4, 5.0 ms — spread 7·10⁻¹⁵, i.e. the data prefer **none**.
For `M_c` the residual depends sharply on the candidate and is minimal **exactly**
at the injected value. So the test works — and that is precisely why the flat
likelihood for `τ₀` is a result rather than an instrument failure.

---

## 4b. H3 on the real GW150914 data (H1 and L1) — done in round 219

Before this round H3 was checked on an **analytic** PSD calibrated to the
published SNR 20. Now it is checked on a **measured** PSD and a physical
amplitude (`M_c = 28.096 M_⊙`, `D_L = 410 Mpc`), with no calibration constant.

**How the PSD is measured.** Median Welch over 1 s Hann segments with 50 %
overlap, on **off-source** stretches (0–8 s and 24–32 s; the event is at 15.4 s).
The resulting ASD(100 Hz) = **1.03·10⁻²³ /√Hz** (H1) and **9.85·10⁻²⁴** (L1) —
the published aLIGO O1 sensitivity. The median is chosen **by measurement**, not
by taste: the mean gives 2.3·10⁻²², 22× higher, because the release contains
non-Gaussian glitches to which the mean is not robust. A consistency check
(`psd_whiten_check.py`): a synthetic series generated from this PSD reproduces
its own band integral with a coefficient of 0.98, while the real off-source data
carry about 6× more band power — i.e. the excess is **in the data**, not in the
estimator.

**The match filter, and why it can be trusted.** The normalization of the
statistic is fixed **by measurement** (`mf_norm_probe.py`), not by algebra: with
`X_k` an ordinary rfft, a template in the same discrete convention
`H_k = fs·h~(f_k)` and `C = 4·df/fs²`, the statistic gives exactly the optimal
SNR. Measured behaviour:

| control | result | what it proves |
|---|---|---|
| pure synthetic noise (64 realizations) | peak ρ = 4.0 on average, max 4.8 | the statistic is calibrated (Rayleigh maximum) |
| injection with target SNR 20 | returned **20.1** on average, timing to 0.12 ms | the filter is neither blind nor biased |
| real H1 data, 4 s around the event | peak **7.27** at GPS **1126259462.42** | it sees the event |
| real L1 data | peak **5.58**, H1−L1 offset = **0.37 ms** | the two observatories agree |
| the same filter on an off-source stretch | peak 3.4 (H1) / 3.6 (L1) | the event stands out above the noise |

The published single-detector SNR of GW150914 ≈ 20 refers to the **full** IMR
template; here the template is the leading-order inspiral only, truncated at
300 Hz, without merger and ringdown, so it cannot return 20. This is named, not
hidden: criterion S6 is phrased as “above the threshold of 5 and above its own
off-source peak”, and both measured numbers are given side by side.

**H1/H2/H3 on the measured PSD.** Both detectors:

| quantity | H1 | L1 | verdict |
|---|---:|---|---|
| fraction of `τ₀` (1.2 ms) absorbed by the source | 1 − 8.9·10⁻¹⁶ | 1 − 1.1·10⁻¹⁵ | **H1 confirmed** |
| unabsorbed SNR for `τ₀` | 3.4·10⁻¹⁴ | 4.8·10⁻¹⁴ | **H1 confirmed** |
| unabsorbed SNR for `τ̇` (0.01) | 1.9·10⁻¹⁵ | 6.4·10⁻¹⁵ | **H2 confirmed** |
| absorbed fraction of the quadratic term at 1.2 ms | 0.8849 | 0.8716 | H3: not absorbed |
| SNR of the quadratic term at 1.2 ms | **0.174** | **0.127** | **H3 confirmed** |
| detectable amplitude at SNR 5 | **34.5 ms** | **47.3 ms** | 29× and 39× the claimed value |

So on the **real** data, with a measured PSD and a physical amplitude, the
structural conclusion survives: the constant and linear delays are absorbed to
the 10⁻¹⁴ level, and the first non-degenerate (quadratic) term at the claimed
1.2 ms gives an SNR of 0.13–0.17 — zero. Thresholds: **15 of 15**.

An independent check of this result is `verify_h3_eigh.py`: a different PSD
estimator, delay directions built from the **analytic** expansion in
`x = (πM_c f)^{2/3}` (not from a numerical `t(f)`), projection through the
**eigenvalues** of the Gram matrix (not QR). Agreement with the artifact:
absorbed fractions 0.88487345 and 0.87164874, residual SNRs agree to 1e-6.

Artifact: `work/artifacts/h3_real_data.json`. Script: `work/h3_real_data.py`.

---

## 5. H4: the arithmetic of the paper itself — including the refutation of my own claim

**H4a — stands, unchanged.** The paper gives two numbers — 0.6 ms (“dimensional
analysis”, `t = R_s/c`) and 1.2 ms (“Schwarzschild-metric computation”).
**Measured: these are not two independent estimates.** The ratio of the second to
the first is `(ΔV_total)^{1/3}/R_s(62) = 2.0083`; the paper’s own ratio is 2.0.
So the “agreement” of 0.6 and 1.2 ms is arithmetic, not confirmation. Check:
`(ΔV)^{1/3}/c = 1.2269 ms` is reproduced exactly, `R_s/c = 0.611 ms`.

**H4b — REVISED, and the claim is weakened (round 219).** A direct reading of the
paper’s text (not from memory) gives: `0.0006 s` → “**of order 1 millisecond**”,
and `0.0012 s` → `0.0015 s`. So the paper states only **one** exact multiplier:
`1.5/1.2 = 1.25`. The apparent second multiplier `1.0/0.6 = 1.667` is an artifact
of the rounded phrase “of order 1 millisecond”, not a derived number. Counting it
as a second stated multiplier is a misreading of the text.

**H4c — REFUTED, and this is my own error from round 217.** I asserted that the
multiplier 1.25 could be reached “only by an unjustified choice of radius”
`r = 2.78 Rₛ`. **That is false.** The volume-weighted mean redshift over the
paper’s **own stated shell** `[1.5, 4] Rₛ` is **1.25034** — i.e. the paper’s 1.25
to 0.03 %. A natural weighting (created time proportional to created volume)
reproduces the paper’s unstated number from its own stated geometry. The root
`r = 25/9 Rₛ` is a false lead: it is the radius where the *point* value of `1+z`
equals 1.25, not what volume averaging gives.

Artifact: `work/artifacts/h4_paper_text_audit.json`. It quotes the old claim next
to its refutation. Independent check: `verify_independent_v3.py` recomputes the
volume mean **symbolically** (sympy, substituting `r = Rₛ₆₂·u`), agreement 1e-6.

**What remains of H4.** H4a: the paper’s two numbers are one number multiplied by
a geometric factor, so the “encouraging agreement” between them is arithmetic.
And 1.25, though now explainable, is still introduced in the text **without
derivation**: the reader is given the number and told only “with the additional
gravitational redshift”.

**Negative controls.** NC1: if the volumes are equal the multiplier is 1 and H4a
flips. NC2: on the sub-interval `[1.5, 2] Rₛ` the root 1.25 is absent — the scan
confirms it. NC3: the volume mean depends on the shell (1.088 for `[1.5, 10] Rₛ`),
so 1.25034 is a property of the stated shell specifically, not a constant, and the
agreement is not vacuous.

---

## 6. How this is checked (independent paths)

| path | what it does | result |
|---|---|---|
| `verify_lean.sh` | rebuild from scratch + grep for `sorry`/`axiom` + axiom audit | **rc=0**, `NO_SORRY_NO_AXIOM`, only 3 standard axioms |
| `verify_independent.py` | **does not import** `analysis.py`; reads only the artifact; its own code; **QR instead of SVD**; sympy for the algebra | **INDEPENDENT_CONFIRMED**, all values agree |
| `mutation_control.py` | 20 textual mutations in a copy, each must go red | **15/20 killed** |
| `equivalence_check.py` | the 4 survivors: rerun + numerical diff of the artifact | all 4 **equivalent on the headline claims** |
| `verify_h3_joint.py` | **round 227**: its own PSD, its own GMST, its own beam algebra, an **explicit matrix** `exp(2πift)` instead of an inverse DFT, its own quadrature | **VERIFY_CONFIRMED**; V2 rel 3.9·10⁻¹¹, V3 7.3·10⁻¹¹, V4 1.4·10⁻⁵, V5 7.5·10⁻⁶, V6 **0.0**; self-test goes red on a corrupted artifact |
| `verify_h3_round226.py`, `verify_h3_imr.py`, `verify_rho_opt_bands.py`, `verify_new_results.py` | independent verifiers of rounds 220–226 | all `VERIFY_CONFIRMED`, self-tests go red |

**The independent probe found an error of mine.** In checking the `τ̇` direction I
wrote `−M_c·∂Ψ/∂M_c`; sympy gives `+M_c·∂Ψ/∂M_c`. The algebra in the analysis was
correct; the text of the check was wrong. Exactly the class of error the
independent path exists for.

**The mutation control found 8 holes in my own test.** The first pass gave 11/20.
The survivors were mutations that **did not change the ratios** on which the
thresholds stood: the scale of `t(f)`, the scale of `h(f)`, removal of the PSD
calibration, a factor of 2 in `R_s`, the distance `D_L`, the coefficient
`∂Ψ/∂M_c`, the chirp coefficient, the sign of `ΔΨ`. I added six self-checks S1–S6
pinning the **absolute** numbers through independent algebraic paths (S6 — an
exact identity: the direction of a purely linear delay must be antiparallel to
`∂Ψ/∂M_c`, cosine `1.0000000000000002`). Second pass: **15/20**. The four
remaining survivors were classified **by measurement**, not by assertion: all
headline quantities agree with the baseline; only the Fisher condition number and
a ratio derived from it differ — quantities that depend on the scale of the basis.

**The automated `mutation_test` campaign returned a null result.** Its coverage
tracer reported “0 executed lines” and declared all 20 mutants unpassed without
running any. That is an **instrument failure, not evidence** — which is why the
mutation control was done by hand. The campaign is frozen as `mut_dcd2ff7529de`.

---

## 7. What I did NOT do, and why

1. **Raw LIGO data — used, provenance closed.** Both GW150914 files (H1 and L1,
   32 s, 16 kHz) were read; `h5py` 3.16.0 is available. `gwosc.org` and
   `dcc.ligo.org` are unreachable from the container (timeout), so the owner
   downloaded the files from the AEI mirror; the hashes match the official GWOSC
   files (`81040e1e…dd1b97`, `23a20702…e41b1f`) — confirmed with
   `file_fingerprint` on the local bytes. Provenance — `data/DATA_PROVENANCE.md`.
   **The skymap is from LOSC** (`LALInference_skymap.fits.gz`, P1500227),
   registered as source `src_e2d7dc79a4df`.
2. **Waveform — IMRPhenomT** (rounds 224–227), no longer leading order. The H1/H2
   conclusion does not depend on the form (it is algebraic and exact). The H3
   numbers do, and this is measured: 0.087 (leading order) versus the values on
   the full IMR.
3. **PSD — measured**, from the off-source stretches of the data itself:
   ASD(100 Hz) = 1.03·10⁻²³ (H1) and 9.85·10⁻²⁴ (L1) /√Hz, matching the published
   aLIGO O1 sensitivity. The estimator is a **median** Welch; with the mean the
   numbers move by tens of percent because of non-Gaussian glitches in the O1
   release. The sensitivity is named as a number (round 219).
4. **No claims about consciousness, free will, or the “reality of the present”.**
   The video transcript is not a source of physical claims.
5. **The Fisher condition number** is given as an order of magnitude (> 10¹²,
   numerically singular), not as a measurement: the mutation control showed it is
   not invariant under a renormalization of the basis.
6. **The sky position is not re-measured from the data.** It is read from the
   published LALInference map. The project does not build its own sky position —
   that is a separate task, and its result is not needed for the identifiability
   question.
7. **The hand-written HEALPix implementation was deleted, not debugged.** Two
   attempts gave wrong geometry (pixel 0 at `θ = 0.09°` instead of `89.93°`).
   Instead of debugging, two independent libraries were used plus a physical test
   of the reading (round 227; `work/NOTES.md`, chronology).

## 8. What this changes about the original question

The README asked: “is the prediction specified precisely enough to be
distinguished from GR and from a re-fit of the source parameters?” The answer:

* **From a source-parameter re-fit — no, and it cannot be.** `τ₀` and `τ̇` lie in
  the linear span of `(M_c, t_c)`. This is a theorem, not an accuracy estimate.
* **From GR — only through the quadratic and higher terms**, i.e. through an
  effect of order `(v/c)⁸` relative to the leading phase. At leading order this
  requires a delay amplitude of ~69 ms instead of the claimed 1.2 ms (a
  leading-order number; it does not survive a full IMR template). That is three
  orders of magnitude above what the mechanism itself provides:
  `R_s(62)/c = 0.61 ms` is an upper bound on the scale.
* **Consequence for the paper.** Its own recommendation (§5) — obtain all
  parameters from the early part of the signal and predict the phase of the peak
  — **does not work** for this effect: the early part of the signal is exactly
  what yields `M_c` and `t_c`, and the effect by construction looks like a change
  in them. The paper’s caveat is correct, but for a stronger reason than it
  states.
* **What round 227 added.** The coherent joint H1+L1 fit with antenna responses
  does not change the structural conclusion (it is algebraic), but **strengthens**
  the numerical side: the peak network SNR 22.93/23.07 reproduces the published
  **network** value 24 to ~4 %, and on the ±38.4 ms grid a right crossing of
  `peak − 1σ` on `τ₂` appears (3.50σ / 2.35σ) — **an upper estimate, not a
  measured bound** — and with a power control showing that the fit sees 20 ms but
  not 1.2 ms. And two structural facts: `φ_c` is not identifiable, and the whole
  detector response collapses to one complex coefficient.

**The line I do not cross.** Everything above concerns the **identifiability of
the prediction**, not the truth of the “now”-theory. The non-degeneracy of the
effect’s direction says nothing about whether the theory is true. It says only
that LIGO cannot test it.

---

## 9. Files

| file | what |
|---|---|
| `work/PREREGISTRATION.md` | pre-registration, frozen before the run |
| `work/CLAIMS.md` | **claim ledger: proven / derived symbolically / measured under assumptions / hypothesis** |
| `work/analysis.py` | main numerical analysis (v5), analytic PSD |
| `work/h3_real_data.py` | H3 on the real H1/L1 data, measured PSD, match filter |
| `work/h3_joint_fit.py` | **round 227: joint H1+L1 fit, one `τ₂`, antenna responses `F₊,F×`** |
| `work/sky_prep.py` | **round 227: skymap reading with two libraries, reading chosen by physics** |
| `work/healpix_nested.py` | helper module for RING/NESTED geometry (not used in the fit) |
| `work/verify_claims_228.py` | **round 228: machine check of the ledger tags (does A4–A6 have Lean?) + definition of “+3.6 ms”** |
| `work/verify_h3_joint.py` | **round 227: independent check of the joint fit + self-test** |
| `work/h3_imr_spins.py` | round 226: spins and inclination over the published envelope |
| `work/h3_direct_fit_imr.py` | round 226: fitting `τ₂` to the record with the full IMR |
| `work/h3_imr_check.py` | round 224: the 37.3↔7.27 gap decomposed with the full IMR |
| `work/h3_rho_opt_bands.py` | round 221: where 37.3 comes from + the “B2 = B3?” test |
| `work/h3_grid_bounds.py` | round 221: grid expanded + physical bound 1.667 ms |
| `work/h3_direct_fit.py` | round 220: direct fit of the delay to the real data |
| `work/pn_orders.py` | round 220: identities at all PN orders |
| `work/et_ce_scaling.py` | round 220: scaling to ET/CE, two conventions |
| `work/h4_paper_arithmetic.py` | audit of the paper’s arithmetic (v2, corrected) |
| `work/h4_paper_text_audit.py` | H4 v3: audit against the paper’s text, volume mean |
| `work/h4_weight_audit.py` | round 220: 1.25 depends on the choice of weight |
| `work/prereg_t8.py` | round 220: T8 per the pre-registration |
| `work/mf_norm_probe.py`, `work/psd_whiten_check.py` | match-filter normalization, PSD estimator |
| `work/verify_independent.py`, `verify_h4_independent.py`, `verify_independent_v3.py`, `verify_h3_eigh.py` | independent verifiers (QR + sympy, eigenvalues) |
| `work/verify_rho_opt_bands.py`, `verify_new_results.py`, `verify_h3_imr.py`, `verify_h3_round226.py` | independent verifiers of rounds 220–226 |
| `work/verify_lean.sh`, `verify_lean_cosmo.sh` | Lean build + axiom audit + negative controls |
| `work/mutation_control.py`, `work/equivalence_check.py` | 20 mutations; survivors classified by measurement |
| `work/lean/*.lean` | `Identifiability`, `Instances`, `RelativityOfNow`, `PastHypothesis`, `HomogeneousUniverse` |
| `work/NOTES.md` | **the full chronology and every error v1 → round 227** |
| `work/artifacts/h3_joint_fit.json` | **round 227: joint-fit artifact** (`art_32983e9ae39a`) |
| `work/artifacts/sky_samples.npz` | **round 227: frozen sky samples** (`art_d8f822f55257`) |
| `work/artifacts/*.json` | artifacts of the other rounds |
| `work/sky/LALInference_skymap.fits` | GW150914 skymap (LOSC P1500227) |
| `data/*.hdf5` | raw GW150914 data (H1, L1), supplied by the owner |

**Reproduction commands** (from `/work/time_flow_research_20260927/work`):

```sh
python3 analysis.py                    # exit 0, 22/22 thresholds
python3 h3_real_data.py                # exit 0, 15/15 thresholds on real data
python3 h3_joint_fit.py                # exit 0, 9/9 controls (round 227)
/work/env/venv/bin/python sky_prep.py  # exit 0, accepted=true (needs healpy+astropy_healpix)
python3 verify_h3_joint.py             # VERIFY_CONFIRMED (round 227)
python3 verify_h3_joint.py --selftest  # corrupt the artifact → goes red
python3 verify_claims_228.py           # VERIFY_CONFIRMED (round 228: ledger tags)
python3 verify_claims_228.py --selftest # corrupt a copy of the artifact → goes red
python3 h3_imr_spins.py                # exit 0, 8/8 controls
python3 h3_direct_fit_imr.py           # exit 0, 8/8 controls
python3 h3_imr_check.py                # exit 0
python3 h3_rho_opt_bands.py            # exit 0
python3 h3_grid_bounds.py              # exit 0
python3 h4_paper_arithmetic.py         # exit 0 (v2)
python3 h4_paper_text_audit.py         # exit 0 (v3)
python3 h4_weight_audit.py             # exit 0
python3 prereg_t8.py                   # exit 0
sh verify_lean.sh                      # rc=0, NO_SORRY_NO_AXIOM
sh verify_lean_cosmo.sh                # rc=0, axiom audit + 3 negative controls
python3 mutation_control.py            # exit 1 (pre-existing: 15/20)
```

**A note on `sky_prep.py`.** It needs `healpy` and `astropy_healpix`, which are
not in the system `python3`; they are installed in the local venv
`work/env/venv` (`package_install`). The artifact `sky_samples.npz` is frozen,
and the main fit reads exactly that — so the reproducibility of the fit does not
depend on those libraries.

## 10. Questions requiring your decision

1. **Raw LIGO data — CLOSED (round 220).** The owner downloaded both files
   directly from GWOSC and compared hashes: they matched the AEI mirror one for
   one. I confirmed the local bytes independently (`file_fingerprint`):
   `81040e1e…dd1b97` (4 068 797 b), `23a20702…e41b1f` (3 919 118 b). The caveat
   “the mirror is not authenticated by GWOSC” is **withdrawn** — it was a legal,
   not a scientific, question, and the owner is right: the correctness of the
   data was already confirmed by three independent physical checks (ASD at
   100 Hz, peak time, H1−L1).
2. **Where to go next.** Variant (b) — “arrow of time from boundary conditions” —
   is **done** (§11). Variant (c) — homogeneous universe — is **done**
   (`lean/HomogeneousUniverse.lean`). Variant (a) — higher PN and ET/CE — is
   **done** (`pn_orders.py`, `et_ce_scaling.py`) and hits a limit: without the
   function `τ(t)` the sign of the wide-band benefit is undefined. The joint
   H1+L1 fit with antenna responses (your msg_00227) is **done**
   (`h3_joint_fit.py`, `verify_h3_joint.py`; chronology — `work/NOTES.md`).
   Untouched and requiring your decision: **derive the specific form `τ(t)` from
   the theory itself** (item 1 of your review). Until it exists, what is proven
   is degeneracy *within the adopted linear model*, not the impossibility of
   testing the whole “now”-theory.
3. **Installation.** Nothing was installed into the system interpreter. For
   reading the skymap, `healpy`, `astropy-healpix` and `phenomxpy` were installed
   into the **local venv** `work/env/venv` (via `package_install`); the venv is a
   plain folder the owner can delete. The artifact `sky_samples.npz` is frozen,
   so the main fit does not depend on these libraries.

---

## 11. Variant B — the arrow of time from boundary conditions (Lean 4)

**What is formalized.** Microscopic laws are time-reversible and therefore do not
by themselves distinguish past from future. The standard answer (Boltzmann; the
“past hypothesis”, Albert 2000): the arrow comes **not from the dynamics but from
a boundary condition**. I formalized the logical skeleton of that answer so that
it is a theorem, not a slogan.

Model: `x : ℤ → α` — a bi-infinite trajectory of reversible dynamics `T` (a
bijection); `rev R x` — time reversal under an involution `R`; `TRSymmetric T R`
— the dynamics is T-symmetric (`R T = T⁻¹ R`); `IsArrow R a` — an “arrow”: a
function on trajectories that changes sign under time reversal.

| theorem | content |
|---|---|
| `rev_isTraj` | **the reversal of a trajectory is again a trajectory**: histories of T-symmetric dynamics are closed under time reversal |
| `arrow_vanishes_on_symmetric` | a trajectory equal to its own reversal carries a **zero** arrow |
| `rev_allowed_of_invariant` | a T-invariant boundary condition preserves closure of the admissible set |
| `arrow_sum_zero` | **conclusion.** On any finite family of trajectories closed under reversal, the arrow sums to **exactly zero** |

**Meaning of the conclusion.** If the boundary is time-symmetric, the arrow cannot
“point” anywhere: the total contribution is exactly zero. A non-zero arrow
**requires** a boundary that **breaks** T-symmetry. This is the logical content of
the past hypothesis.

**Non-vacuity (mandatory — otherwise the theorem is useless).** Three concrete
witnesses in the same formal field: `demo_arrow_vanishes` — an explicit trajectory
carrying a zero arrow; `paramArrow_isArrow` + `paramArrow_ne_zero` — a **non-zero**
arrow really exists (the value of a trajectory at time 0), so `arrow_sum_zero` is
not vacuous; `demo_closure_is_the_line` — the **decisive witness**: the same
non-zero arrow on two families differing only in closure under reversal gives the
sum `0` (closed) versus `2` (not closed). So the closure hypothesis is **exactly
the line** that separates the conclusion, not a technical condition.
`Bfuture_not_invariant` + `demo_boundary_breaks_reversal` — the boundary
`{n | 0 ≤ n}` is **not** T-invariant, so the premise of the conclusion is indeed
violable.

**Check.** `sh verify_lean.sh` — rc=0, `NO_SORRY_NO_AXIOM`, axiom audit: all 12
theorems depend only on the standard `[propext, Classical.choice, Quot.sound]`; a
scan for `native_decide`/`unsafe`/`implemented_by` is empty. **Negative controls**
(all go red): remove the `TRSymmetric` hypothesis from `rev_isTraj` → compile
error; remove the closure hypothesis from `arrow_sum_zero` → 4 compile errors;
insert a `sorry` → the scan catches it.

**The line I do not cross.** This is the **logic** of the past-hypothesis argument,
not the physics of entropy growth. Nothing here asserts that the universe
**really** had a low-entropy beginning; what is asserted is: *if* the laws are
T-symmetric, then the arrow cannot come from them and must come from a
symmetry-breaking boundary.

---

## 11b. Joint H1+L1 fit with antenna responses (round 227)

This is the current state of H3 after your msg_00227: instead of scanning the
inclination `ι` on one detector — a **coherent network** fit with **one** `τ₂` and
free `(M_c, t_c, φ_c)`, where the responses `F₊, F×` are taken for the published
GW150914 sky position.

**Two structural facts without which the fit cannot be assembled.**

1. The phenomxpy convention (measured): `hp = A(ι)·e^{2iφ}·hp₀`,
   `hc = B(ι)·e^{2iφ}·hc₀`, `hc₀ = −i·hp₀` (deviation 3.3·10⁻⁴). Hence the
   detector strain `F₊·hp + F×·hc = e^{2iφ}·hp₀·(F₊A − iF×B)`, i.e. **the whole
   response collapses to one complex coefficient** `C_d = F₊,d·A − iF×,d·B`.
2. Since `h_d = C_d·h` with a common `h`, the coherent statistic
   `ρ_coh² = |Σ_d conj(C_d)ζ_d|² / Σ_d |C_d|²s_d²`, where `ζ_d` does not depend on
   `C_d` — computed once per `(M_c, τ₂)`, after which the `(ι, ψ)` grid is free.
   The amplitude and phase ratio of the two detectors is **predicted** by the sky
   position.

**A consequence named as a finding.** `φ_c` multiplies `C_d` of each detector by
the same phase, so it is **not identifiable** in a coherent fit with a fixed sky
position; it is profiled, not quoted.

**Sky position.** From the published LOSC map `LALInference_skymap.fits.gz`
(P1500227). Two independent HEALPix readers (`healpy`, `astropy_healpix`) agree to
1.8·10⁻¹⁵ rad. The correctness of the **reading** of the FITS table was chosen by
physics: the accepted reading places 94 % of the probability on the published
6.9 ms annulus, at **25×** the level of a shuffled map. ML position
RA = 134.80°, Dec = −69.79°; 90 % area 616.4 deg² (published 610);
`τ_H1 − τ_L1 = −6.898 ms` (published 6.9).

**Result.**

| delay kernel | grid | profile drop | right crossing of `peak − 1σ` |
|---|---|---:|---|
| pn | ±9.6 ms | 1.375σ | +5.4 ms (peak at the **left edge**) |
| pn | **±38.4 ms** | **3.496σ** | **+3.6 ms** (peak interior) |
| imr | ±9.6 ms | 0.852σ | none (profile does not drop by 1σ) |
| imr | **±38.4 ms** | **2.352σ** | **+2.4 ms** (peak interior) |

**What “+3.6 ms” means.** It is the **right crossing of the level `peak − 1σ`** —
the first grid point to the right of the peak where the profile falls 1.0 below
the peak. It is neither the grid edge nor a confidence interval. On the ±9.6 ms
grid the peak sits at the **left edge** (`−9.6 ms`, where the profile is still
rising), so the right crossing is **+5.4 ms**; on ±38.4 ms the peak is interior
(`−18.0 ms`) and the crossing is **+3.6 ms**. There is no left crossing on either
grid. So “a bound of +3.6 ms” is a right crossing at an **interior** peak; it
lands inside ±9.6 ms only because on the narrow grid the peak runs to the edge.

The peak network SNR is **22.93** (pn) / **23.07** (imr); at `τ₂ = 0` — **22.17**;
the published **network** value is 24. The ratio of optimal single-detector SNRs
is **1.4749** against the published single-detector **1.4662** — agreement 0.6 %,
but this is **not** a 0.6 %-precision test: the published 19.5 and 13.3 are
observed quantities with noise ±1, so the ratio carries a ~7 % tolerance. The
real checks on the position are the delay `−6.898 ms` (published 6.9), the area
616.4 deg² (published 610) and the shuffle control (25×).

**Power control.** The same fit recovers a 20 ms injection (20.3 ± 5.1 ms) and
**does not resolve** 1.2 ms (spread 4.6 ms at an amplitude of 1.2 ms).

**Independent check** (`verify_h3_joint.py`): its own PSD, its own GMST, its own
beam algebra, an **explicit matrix** `exp(2πift)` instead of an inverse DFT, its
own quadrature. V2 — rel 3.9·10⁻¹¹, V3 — 7.3·10⁻¹¹, V4 — 1.4·10⁻⁵, V5 —
7.5·10⁻⁶, V6 — **0.0**; the self-test goes red on a corrupted artifact. The
artifacts `art_32983e9ae39a` and `art_d8f822f55257` reproduce byte-for-byte.

**The errors of this round — six, all in `NOTES.md`.** The most useful: the
hand-written HEALPix implementation twice gave wrong geometry and was **deleted**,
not debugged.

---

## 12. Chronology of rounds — moved to `work/NOTES.md`

Sections §12–§18 of earlier report editions were a chronicle (“round 220”, “round
221”, “round 224”, “round 226”, “round 227”) and retold what is already said in
§1–§10. At the owner’s request (msg_00227) the chronology was moved to
`work/NOTES.md`, and the report is kept as a document of the **current state**:
what is proven, what is measured under which assumptions, what is not done, which
files exist and how to reproduce.

Full chronology — `work/NOTES.md`, section “Chronology of rounds 216–227”.
Claim tagging — `work/CLAIMS.md`.
