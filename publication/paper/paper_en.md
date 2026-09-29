---
title: "Is the Müller–Maguire “now”-theory prediction falsifiable on GW150914?"
subtitle: "A formal identifiability audit, with a numerical check on real LIGO strain"
author:
  - "Aiodam (autonomous research agent)"
  - "Oxunjon Ubaydullayev"
date: "2026-09-28"
lang: en
abstract: |
  Müller and Maguire (arXiv:1606.07975) estimate an additional signal delay of
  order 1–1.5 ms at black-hole merger and call the calculation heuristic. We ask
  a narrower, decidable question: **is that prediction specified precisely enough
  to be distinguished from general relativity and from a re-fit of the source
  parameters?** We answer it in three independent registers. (i) *Formally*, in
  Lean 4 with Mathlib, we prove an identifiability criterion for an additive
  model M(θ,τ) = Aθ + Bτ and show it is a genuine dichotomy. (ii) *Symbolically*,
  we show that a constant delay τ₀ and a linear delay τ̇ are **exactly**
  reparameterizations of the coalescence time t_c and the chirp mass M_c at
  **every** post-Newtonian order, so the degeneracy is not an artifact of the
  leading order. (iii) *Numerically*, on the real GW150914 strain (H1 and L1, with
  a PSD measured from the off-source data and a physical amplitude), a coherent
  two-detector fit with one τ₂ recovers a 20 ms injected delay but cannot resolve
  the claimed 1.2 ms (fit resolution ≈ 5 ms). The first non-degenerate term is
  quadratic and suppressed as (v/c)⁸; at the claimed 1.2 ms it yields a
  signal-to-noise ratio of 0.13–0.17, i.e. zero. We also audit the arithmetic of
  the paper itself: its two “independent” numbers (0.6 ms and 1.2 ms) differ by
  exactly the geometric factor 2.0083, and its unstated multiplier 1.25 is
  reproduced to 0.03 % by volume-weighting its own stated shell — but only under
  a weighting the paper never justifies. **Scope of the claim:** what is proven
  is degeneracy *within the adopted linear model*; the specific function τ(t) is
  not derived from the theory, because the paper does not specify it. This
  strengthens the paper’s own caveat by translating it from prose into a theorem,
  and narrows it: the obstacle is not the quality of one event, it is that the
  effect lies inside the linear span of the source parameters.
  All code, data provenance, artifacts and negative controls are included; every
  numerical claim is reproduced by an independent verifier that shares no code
  with the audit.
---

# 0. Authorship and disclosure

This work was carried out **by an autonomous research agent, Aiodam**, acting as
the primary investigator: it designed the study, wrote and ran every script,
found and corrected its own errors, and wrote this report. The human co-author,
**Oxunjon Ubaydullayev**, set the research question, supplied the LIGO data files
and the skymap, and reviewed the intermediate results across several rounds.

**How the errors were found, counted rather than remembered.** `NOTES.md` records
47 errors made on the way. The co-author found **eleven** of them: **six** in one
review of the claim ledger (one of those by reading the ledger's own legend rather
than by checking any number), **one** earlier, when he read an intermediate report
and caught a sampling artifact in the H4 arithmetic, and **four** further errors of
formulation. The remaining **36** were found by the agent's own controls, verifiers
and self-tests. The split is produced by a script (`count_errors.py`), not by
memory, and the script prints the rule it used so the rule can be argued with.

We state this plainly because it is part of the claim: the whole path — from the
question to the negative controls — was walked by the agent. Where the agent was
wrong, the errors are recorded in `NOTES.md` rather than removed, and every
headline number is reproduced by a verifier that does not import the audit code.

# 1. The question

Müller and Maguire (2016) argue, from a thermodynamical arrow-of-time picture,
that the “now” of an observer is not a point but a growing block, and estimate
that a gravitational-wave signal from a black-hole merger should carry an
additional time delay of order 1–1.5 ms. They explicitly call the calculation
heuristic and do not present a completed modification of general relativity.

The question we can actually decide is not “is the theory true” but:

> **Is the predicted delay specified precisely enough that GW150914 could
> distinguish it from (a) general relativity and (b) a re-fit of the source
> parameters?**

A prediction that can be absorbed exactly by a shift of parameters the data
already leaves free is not testable by that data, however large the detector’s
sensitivity. That is the object of study here.

# 2. Method and pre-registration

The test plan was frozen before the first run (`PREREGISTRATION.md`). It fixes
the hypotheses, the thresholds, and the negative controls. Where the plan and the
implementation diverged, the divergence is recorded and the added or changed
checks are marked as such (see §9, item 7).

Three registers are used, deliberately kept separate:

1. **Formal** — Lean 4 + Mathlib, no `sorry`, no custom axioms. Independent of
   data, noise and waveform.
2. **Symbolic** — closed-form algebra verified with `sympy` and a negative
   control. Independent of data, but *not* formalized in Lean.
3. **Numerical** — measured on the real GW150914 strain, under named assumptions.

The distinction matters: a claim in the **hypothesis** category is never allowed
into the headline, which rests on registers 1 and 2 only. The full ledger of every
claim, by category, is `CLAIMS.md`.

# 3. Formal core: when is a delay identifiable?

## 3.1 The criterion

Consider an additive model

  M(θ, τ) = A θ + B τ,

where θ are the source parameters (chirp mass, coalescence time, phase),
τ is the candidate delay parameter, and A, B are linear maps into the space of
observable waveforms. Then:

> **`identifiable_iff`** (Lean: `lean/Identifiability.lean`).
> τ is identifiable — i.e. M(θ,τ) = M(θ′,τ′) ⟹ τ = τ′ — **if and only if**
> B is injective and range A ⊓ range B = ⊥.

> **`not_identifiable_of_range_le`.** If range B ≤ range A, then τ is not
> identifiable: any effect of τ is absorbed by a shift of θ.

The criterion is a genuine dichotomy, not a vacuous implication: `Instances.lean`
exhibits both cases in the same formal field — `demo_identifiable` (disjoint
directions → identifiable) and `demo_not_identifiable` (delay direction inside
the image → not identifiable).

## 3.2 The formalization is sound

All Lean files compile with `NO_SORRY_NO_AXIOM`, and the axiom audit reports only
the three standard Lean axioms `[propext, Classical.choice, Quot.sound]`. No
custom axioms, no `native_decide`, no `unsafe`, no `implemented_by`.

Three further formal results support the surrounding argument and are stated here
for completeness:

* `RelativityOfNow.lean` — `exists_boost_reversing_time_order`: in 1+1 Minkowski,
  if |Δt| < |Δx| and Δt ≠ 0, there is a boost with |v| < 1 that reverses the sign
  of Δt′. This is the relativity-of-simultaneity objection to a preferred “now”.
* `PastHypothesis.lean` — `arrow_sum_zero`: on any finite family of trajectories
  closed under time reversal, a time-arrow functional sums to exactly zero. A
  non-zero arrow therefore **requires** a boundary condition that breaks
  T-symmetry. The theorem is shown to be non-vacuous by a decisive witness
  (`demo_closure_is_the_line`): the same non-zero arrow on two families differing
  only in closure gives 0 (closed) versus 2 (not closed).
* `HomogeneousUniverse.lean` — with a free base expansion (dark-energy fit) the
  parameter δ is not identifiable; with the expansion fixed it is. The fixing is
  exactly the dividing line.

# 4. Symbolic core: the degeneracy is not a leading-order artifact

## 4.1 The algebra

With x = (π M_c f)^{2/3}, t(f) − t_c = −(5/256) M_c x^{−4} and 2πf = 2x^{3/2}/M_c,
expanding τ(t) = τ₀ + τ₁(t − t_c) + τ₂(t − t_c)² + … gives a phase correction

```
ΔΨ =  2 x^{3/2} τ₀ / M_c                 ∝ ∂Ψ/∂t_c        ← ABSORBED
    − (5/128) τ₁ x^{−5/2}                ∝ ∂Ψ/∂M_c        ← ABSORBED
    + (25/32768) τ₂ M_c x^{−13/2} + …    ← first non-degenerate
```

The leading phase is Ψ_chirp = (3/128) x^{−5/2}, with
∂Ψ/∂M_c = −(5/128) x^{−5/2}/M_c. The first two terms are **exact**
reparameterizations, not approximations: the constant delay is a shift of t_c,
the linear delay is a rescaling of M_c. The third term sits on x^{−13/2}, four
powers of x^{−4} ~ (v/c)⁸ below the leading term.

## 4.2 The identities hold at every PN order

The key point is that the degeneracy does not depend on the leading-order
truncation. Two identities do the work:

  ∂Ψ/∂t_c = 2πf,  and  f · ∂Ψ_chirp/∂f = M_c · ∂Ψ_chirp/∂M_c.

The second is Euler’s theorem applied to v, because Ψ_chirp depends on M_c and f
**only** through v — **no PN coefficient enters**. Verified symbolically with
`sympy` at PN orders 0 … 3.5, with a negative control (the control was first
vacuous — `ln(M_c·f)` is a function of v — and was caught because it failed to go
red).

Consequently the linear part of any delay, at any PN order, is exactly a
combination of the source directions, and the degeneracy is structural.

# 5. Numerical results on the real GW150914 strain

## 5.1 Data and provenance

Both GW150914 strain files (H1 and L1, 32 s, 16 kHz, 524 288 samples each) were
used. They were obtained from the AEI mirror because `gwosc.org` is unreachable
from the container; their SHA-256 hashes (`81040e1e…dd1b97`,
`23a20702…e41b1f`) were confirmed byte-for-byte against the official GWOSC files
by the owner, and re-confirmed locally with `file_fingerprint`. Provenance is
recorded in `data/DATA_PROVENANCE.md`.

![The two detectors around the merger, band-passed 35–350 Hz. The dashed line is the published merger time, GPS 1126259462.4. The strain is plotted as measured, with no whitening and no template subtracted.](../figures/fig1_strain.pdf)

## 5.2 The noise curve is measured, not assumed

The PSD is a **median** Welch estimate over 1 s Hann segments (50 % overlap) on
**off-source** stretches (0–8 s and 24–32 s; the event is at 15.4 s). The
resulting amplitude spectral densities at 100 Hz are

| detector | ASD(100 Hz) | note |
|---|---:|---|
| H1 | 1.03·10⁻²³ /√Hz | matches published aLIGO O1 sensitivity |
| L1 | 9.85·10⁻²⁴ /√Hz | matches published aLIGO O1 sensitivity |

The median is chosen by measurement, not taste: the mean gives 2.3·10⁻²², 22×
higher, because the O1 release contains non-Gaussian glitches to which the mean
is not robust.

![The measured amplitude spectral density from off-source data (median Welch), for both detectors, with the analysis band 20–300 Hz shaded. The dotted lines are the artifact's own ASD(100 Hz) values.](../figures/fig2_asd.pdf)

## 5.3 The match filter is calibrated by measurement

The normalization of the matched-filter statistic was fixed **by measurement**,
not by algebra: with a pure synthetic noise (64 realizations) the peak ρ is 4.0
on average (the Rayleigh maximum), and an injection with target SNR 20 is
returned as **20.1** with a timing accuracy of 0.12 ms. On the real data the
filter gives a peak of **7.27** (H1) at GPS 1126259462.42 and **5.58** (L1), with
an H1−L1 offset of **0.37 ms**, against off-source peaks of 3.4 and 3.6.

The peak is 7.27 rather than 20 because the template is the leading-order
inspiral only, truncated at 300 Hz, without merger and ringdown. This is named,
not hidden. The figure below shows the template the later sections actually use —
the full `IMRPhenomT` — against the same whitened strain.

![Whitened H1 strain (20–300 Hz) against the full IMRPhenomT template used in §5.5 and §6, placed at the artifact's recorded IMR peak time. The template is normalised to the data's own amplitude scale, so this is a shape comparison, not an amplitude claim.](../figures/fig3_template.pdf)

## 5.4 H1/H2/H3 on the measured PSD

| quantity | H1 | L1 | verdict |
|---|---:|---|---|
| fraction of τ₀ (1.2 ms) absorbed by the source | 1 − 8.9·10⁻¹⁶ | 1 − 1.1·10⁻¹⁵ | H1 confirmed |
| unabsorbed SNR for τ₀ | 3.4·10⁻¹⁴ | 4.8·10⁻¹⁴ | H1 confirmed |
| unabsorbed SNR for τ̇ | 1.9·10⁻¹⁵ | 6.4·10⁻¹⁵ | H2 confirmed |
| quadratic-term SNR at 1.2 ms | **0.174** | **0.127** | H3 confirmed |
| detectable amplitude at SNR 5 | 34.5 ms | 47.3 ms | 29× / 39× the claimed value |

**15 of 15 thresholds pass.** On the real data, with a measured PSD and a physical
amplitude, the structural conclusion survives: the constant and linear delays are
absorbed to the 10⁻¹⁴ level, and the first non-degenerate (quadratic) term at the
claimed 1.2 ms gives an SNR of 0.13–0.17 — i.e. zero.

## 5.5 How much room is there for *any* quadratic delay?

Because the leading-order template was the first thing measured, it is worth
recording what a full IMR template changes and what physical room is left. On a
full `IMRPhenomT` on the pipeline's own grid, the template's optimal SNR drops from
37.31 to **31.66** (H1) and from 33.71 to **28.74** (L1) — a 15 % reduction, so more
than half of the earlier deficit was a leading-order template artefact. What
remains is the orientation of the source: the template's optimal SNR falls with
inclination and crosses the published single-detector value 19.5 at **ι ≈ 59.2°**.
That crossing is a **fitted** inclination — the value the model would need, not a
measured one — and it is stated as such. This does not change the structural
result (it is algebraic); it changes only the numerical value of the quadratic SNR.

![Template optimal SNR against inclination for H1, from the frozen artifacts. The published single-detector value 19.5 is reached at ι ≈ 59.2° — a fitted inclination, i.e. the value the model would need, not a measurement of this event. The shaded band is the range spanned by the published spin envelope, which moves the number up, away from the target. The residual deficit is orientation, not spins.](../figures/fig4_inclination.pdf)

More usefully, one bound does **not** depend on the scan at all. The delay enters as
the phase `2π f τ(t)`; a delay whose phase sweeps more than π across the band is no
longer a small perturbation of the template but a different waveform. Requiring
`max |2π f τ₂ (t(f) − t_ref)²| ≤ π` over 20–300 Hz caps the quadratic amplitude at
  **1.667 ms** (τ₂ = 0.0023370 s⁻²),

whereas the ceiling set by the 4 s segment duration alone is 4000 ms. The paper’s
claimed 1.2 ms therefore uses **72 %** of the entire budget in which the quadratic
term can still be called a small perturbation. This is a statement about the model,
not about GW150914, but it bounds every future comparison: at 1.667 ms the phase
already reaches π at 300 Hz, and the leading-order expansion is even less
applicable.

**The convention this bound is stated in, because it matters by a factor of 15.**
The bound is computed with the delay vanishing at the **start** of the band and
growing to its maximum at **merger** — the picture the paper itself uses, in which
the effect accumulates toward coalescence and reaches 1.2 ms at the end. If the
quadratic vertex is placed at coalescence instead (so the delay vanishes **at
merger** and is largest at the start of the inspiral — the opposite of the paper’s
picture), the same π-criterion gives **24.96 ms**, fifteen times looser. So “1.667
ms” and “72 % of budget” are statements about the paper’s own accumulating-delay
picture, not convention-independent facts, and they are named as such. (This was
caught by an independent re-derivation that first returned 24.96 ms; the artefact’s
1.667 ms was confirmed correct *for the paper’s convention*, and the discrepancy
became this paragraph rather than a silent choice.)

# 6. Coherent two-detector fit with one τ₂

## 6.1 Two structural facts

Building the joint H1+L1 fit exposed two properties of the problem that were not
visible in the single-detector scan over inclination:

1. **The whole detector projection collapses to one complex coefficient.**
   Measured in the `phenomxpy` convention: hp = A(ι)·e^{2iφ}·hp₀,
   hc = B(ι)·e^{2iφ}·hc₀, with hc₀ = −i·hp₀ (deviation 3.3·10⁻⁴). Hence the
   detector strain F₊·hp + F×·hc = e^{2iφ}·hp₀·(F₊A − iF×B), i.e.
   **C_d = F₊,d·A − iF×,d·B**.
2. **φ_c is not identifiable.** With the sky position fixed, φ_c multiplies every
   detector’s coefficient by the same phase and is exactly degenerate with an
   overall phase. It is profiled out, not quoted. (The owner’s review had asked
   for φ_c to be free; it turned out it structurally cannot be.)

Because ζ_d does not depend on C_d, the coherent statistic
ρ_coh² = |Σ_d conj(C_d)ζ_d|² / Σ_d |C_d|²s_d² can be computed once per (M_c, τ₂),
making the (ι, ψ) grid free, and making the two detectors **coherent**: their
amplitude and phase ratio is *predicted* by the sky position, not free.

## 6.2 Sky position

The position is **read** from the published LALInference skymap
(`LALInference_skymap.fits.gz`, LOSC P1500227), not re-measured from the strain.
Two independent HEALPix readers (`healpy`, `astropy_healpix`) agree to
1.8·10⁻¹⁵ rad. Which FITS-table reading is correct was decided **by physics**:
the accepted reading places 94 % of the probability on the published 6.9 ms
annulus, at 25× the level of a shuffled map. Maximum-likelihood position:
RA = 134.80°, Dec = −69.79°; 90 % area 616.4 deg² (published 610);
τ_H1 − τ_L1 = −6.898 ms (published 6.9 ms).

## 6.3 Result

| delay kernel | grid | profile drop | right crossing of peak − 1σ |
|---|---|---:|---|
| pn | ±9.6 ms | 1.375σ | +5.4 ms (peak at the **left edge**) |
| pn | **±38.4 ms** | **3.496σ** | **+3.6 ms** (peak interior) |
| imr | ±9.6 ms | 0.852σ | none (profile does not drop by 1σ) |
| imr | **±38.4 ms** | **2.352σ** | **+2.4 ms** (peak interior) |

The peak network SNR is **22.93** (pn) / **23.07** (imr); at τ₂ = 0 it is
**22.17**; the published **network** value is 24. The ratio of optimal single-
detector SNRs is 1.4749 against the published 19.5/13.3 = 1.4662 — agreement at
0.6 %, but this is **not** a 0.6 %-precision test: 19.5 and 13.3 are observed
quantities with noise of order ±1, so the ratio carries a ~7 % tolerance. The
real checks on the sky position are the delay (−6.898 ms vs 6.9), the area
(616.4 vs 610 deg²) and the shuffle control (25×).

**What “+3.6 ms” means.** It is the right-hand crossing of the level peak − 1σ —
the first grid point to the right of the peak where the profile falls 1.0 below
the peak. It is neither the grid edge nor a confidence interval. On the ±9.6 ms
grid the peak sits at the **left edge** (−9.6 ms, where the profile is still
rising), so the right crossing is +5.4 ms; on ±38.4 ms the peak is interior
(−18.0 ms) and the crossing is +3.6 ms. There is no left crossing on either grid.
So the statement “a bound of +3.6 ms” is a right crossing at an interior peak,
and it lands inside ±9.6 ms only because on the narrow grid the peak runs to the
edge.

**This is an upper bound on the strength of the constraint, not a measured
bound.** Free η and spins would absorb more of the quadratic term (even with η
fixed and zero spins, 88 % of the term is absorbed by the source parameters), so
the true constraint is weaker. The sky position is fixed at the maximum-
likelihood point (90 % area 616 deg²) and is not fitted jointly with τ₂.

![The joint H1+L1 profile in τ₂ on the ±38.4 ms grid, for both delay kernels, read directly from the frozen artifact. The dashed line is peak − 1σ; the dot marks the profile maximum. On the PN kernel the profile falls 3.496σ, on the IMR kernel 2.352σ. This is an **upper bound on the strength of the constraint**, not a measured bound: free η and spins would absorb more of the quadratic term, and the sky position is held fixed.](../figures/fig5_tau2_profile.pdf)

## 6.4 Power control

The same fit **recovers** a 20 ms injected delay (20.3 ± 5.1 ms) and **does not
resolve** the claimed 1.2 ms (spread 4.6 ms at an amplitude of 1.2 ms). The
current fit resolution is therefore **≈ 5 ms**. Without this control, “not
constrained” would be indistinguishable from “the fit is blind”.

# 7. Auditing the paper’s own arithmetic

The paper gives two numbers, 0.6 ms (“dimensional analysis”, t = R_s/c) and
1.2 ms (“Schwarzschild-metric computation”).

* **The two numbers are not independent.** Their ratio is exactly the geometric
  factor

      (ΔV_total)^{1/3} / R_s(62) = 2.0083,

  and the paper’s own ratio is 2.0. The apparent agreement between 0.6 and 1.2 ms
  is arithmetic, not confirmation.
* **The multiplier 1.25 is reproduced, but not derived.** The volume-weighted
  mean of (1+z) over the paper’s own stated shell [1.5, 4] Rₛ is **1.25034**, i.e.
  the paper’s 1.25 to 0.03 %. But this depends on the weighting: weighting by
  *created* volume V(M) − V(0) gives **1.28596**, and a flat weight gives
  **1.24173**. The paper states the number without justifying the weight.
* **A claim of mine was refuted.** In an earlier round I asserted that 1.25 could
  only be produced by an unjustified choice of radius (r = 2.78 Rₛ). That is
  false: the root 25/9 Rₛ is the radius where the *point* value of 1+z equals
  1.25, whereas the paper averages over volume. The claim is recorded in the
  artifact next to its refutation.

# 8. What is proven, what is measured, what is hypothesis

Every claim in this paper is tagged in `CLAIMS.md` under four categories:

* **[PROVEN]** — formal: Lean 4, no `sorry`, no custom axioms. A1–A3, A7–A9.
* **[DERIVED symbolically]** — closed-form algebra verified with sympy and a
  negative control, but **not** formalized in Lean. A4–A6, B4.
* **[MEASURED]** — a number on real data under named assumptions. B1–B3, B5–B11,
  C1–C5, D1–D4, E6, E7.
* **[HYPOTHESIS]** — not tested or not testable. E1–E5.

The headline of this paper — “the prediction is not falsifiable on GW150914 for a
structural reason” — rests on the **formal** results A1–A3 and on the symbolic
identities A4–A6. No hypothesis-category claim is used in the headline. A machine
check (`verify_claims_228.py`) enforces that the [PROVEN] tag is only applied
where a Lean file with the named theorem exists.

# 9. Limitations and open items

1. **τ(t) is not derived from the theory.** The paper does not specify it. Until
   it is, what is proven is degeneracy *within the adopted linear model*, not the
   impossibility of testing the whole “now”-theory. This is the main open item.
2. **The bound on τ₂ is an upper estimate.** Free η and spins would weaken it.
3. **Aligned spins only.** Precession (χ_p < 0.71) is not covered.
4. **The delay is phase-only**; the amplitude is not delayed.
5. **The ET/CE sign is undetermined.** With a fixed amplitude the benefit of a
   wide band falls (0.087 → 0.0047); with a fixed τ₂ it rises (0.087 → 7.625).
   Without τ(t) the sign of the ET/CE benefit is not defined by the paper.
6. **L1 does not reach the published 13.3** even at ι = 90° (minimum 13.58); the
   residue is attributed to antenna response, which the model lacks.
7. **Pre-registration drift.** The frozen T8 asked for a false-alarm test on 100
   noise realizations; the implementation checks the Fisher condition number.
   The added checks are marked as added.
8. **The π-phase bound (1.667 ms) is convention-dependent.** It assumes the delay
   vanishes at the band start and grows toward merger, which is the paper’s own
   accumulating-delay picture. Placing the quadratic vertex at coalescence instead
   gives 24.96 ms, fifteen times looser. The bound is therefore a statement about
   the paper’s convention, not a convention-free fact, and it is stated that way.
9. **The inclination ι ≈ 59.2° is fitted, not measured.** It is the inclination at
   which the template’s optimal SNR would equal the published single-detector
   value 19.5; it is not a measurement of this event’s orientation.

# 10. Reproducibility

All scripts, artifacts and negative controls are in the repository. Every
headline number is reproduced by an independent verifier that shares no code with
the audit (`verify_h3_joint.py`, `verify_h3_round226.py`, `verify_h3_imr.py`,
`verify_rho_opt_bands.py`, `verify_independent*.py`, `verify_claims_228.py`).
Each verifier has a self-test that corrupts a copy of the artifact and must go
red; a verifier that cannot fail proves nothing.

The mutation control of the main analysis killed 15 of 20 textual mutations; the
four survivors were classified **by measurement** as equivalent on the headline
quantities. An automated mutation campaign returned a null result (its tracer
reported “0 executed lines” and failed all 20 mutants without running any) — that
is an instrument failure, not evidence, and it is recorded as such.

Reproduction commands are listed in the repository README. The main pipeline
runs on a CPU-only container with Python 3, NumPy, SciPy, SymPy, `h5py`,
`phenomxpy`, `healpy` and `astropy-healpix`; the Lean files build against Mathlib.

# 11. What this changes about the original question

* **Distinguishing from a source-parameter re-fit: impossible, and it cannot be
  otherwise.** τ₀ and τ̇ lie in the linear span of (M_c, t_c). This is a theorem,
  not an accuracy estimate.
* **Distinguishing from GR: only through the quadratic and higher terms**, i.e.
  through an effect of order (v/c)⁸ relative to the leading phase. At leading
  order this would require a delay amplitude of ~69 ms instead of the claimed
  1.2 ms (a leading-order number; it does not survive a full IMR template). That
  is **two** orders of magnitude above what the mechanism itself provides:
  R_s(62)/c = 0.61 ms is an upper bound on the scale. (69 ms / 0.61 ms ≈ 113, so
  two orders; the comparison is against the mechanism's own scale, not against
  the paper's 1.2 ms, which is 57× smaller than 69 ms.)
* **Which signal-to-noise ratio is which.** Three different numbers appear in this
  paper and they are not interchangeable: **0.087** is the quadratic-term SNR on
  the *analytic, calibrated* PSD (`analysis.py`); **0.174 / 0.127** are the same
  quantity on the *measured* PSD of H1 / L1 (`h3_real_data.py`); and **0.13–0.17**
  is the range of the pair, quoted in the abstract. They differ because the PSD
  differs (analytic vs measured), not because the physics differs. The joint fit
  of §6 reports **network** SNRs (22.93 / 23.07), which are a different quantity
  again and are compared only with the published **network** value 24.
* **Consequence for the paper.** Its own recommendation (§5) — extract all
  parameters from the early part of the signal and predict the phase of the peak
  — **does not work** for this effect: the early part of the signal is exactly
  what yields M_c and t_c, and the effect by construction looks like a change in
  them. The paper’s caveat is correct, but for a stronger reason than it states.

**The line we do not cross.** Everything above concerns the **identifiability of
the prediction**, not the truth of the “now”-theory. The non-degeneracy of the
effect’s direction says nothing about whether the theory is correct. It says only
that LIGO cannot test it.

# References

1. R. A. Muller, S. Maguire, *Now, and the Flow of Time*, arXiv:1606.07975 (2016).
2. B. P. Abbott et al. (LIGO/Virgo), *Observation of Gravitational Waves from a
   Binary Black Hole Merger*, arXiv:1602.03837, Phys. Rev. Lett. 116, 061102 (2016).
3. B. P. Abbott et al. (LIGO/Virgo), *Properties of the Binary Black Hole Merger
   GW150914*, arXiv:1602.03840, Phys. Rev. Lett. 116, 241102 (2016).
4. B. P. Abbott et al. (LIGO/Virgo), *GW150914: First results from the search for
   binary black hole coalescence with Advanced LIGO*, arXiv:1602.03839 (2016).
5. GWOSC, GW150914 data release v3, <https://gwosc.org/eventapi/html/event/GW150914/v3>
6. LOSC, GW150914 skymap P1500227, https://losc.ligo.org
7. The Mathlib Community, *The Lean Mathematical Library*, CPP 2020.
