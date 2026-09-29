# Round-by-round summaries (rounds 217–228)

*English translation of the Russian round notes `work/msg21X-note.md` and
`work/msg22X-note.md`. The Russian originals are kept alongside; the English text
is the primary one for publication. The full chronology and the full error list
live in `NOTES.md`; this file is the compact per-round record.*

---

## Round 217 — H4 corrected, and variant B

**Date:** 2026-09-27.

### 1. H4 corrected (found by the owner)

**The defect.** Report §5/H4 asserted that the multiplier 1.25 “is not reached at
any choice of radius” from the range 1.5–4.0 Rₛ. **That was false.** The code
`h4_paper_arithmetic.py` (lines 99–100) tested `1+z` at only **four** points
`{1.5, 2, 3, 4}`, where `1+z` = {1.732, 1.414, 1.225, 1.155}, and declared
unreachability from those four samples. The equation `1+z(r) = 1.25` has the root
`r/Rₛ = 1/(1 − 1/1.25²) = 2.7778`, and it lies **inside** `[1.5, 4]`.

**The worst part.** v1 computed `r_needed = 2.7778` itself and printed the range
`[1.5, 4]` in the same artifact — but **never compared them**. Two numbers sat side
by side in one JSON file and the verdict contradicted both. The same class as the
harness defects: a claim that no check ever executed.

**The fix (v2).** A scan of the whole interval (20001 points) finds the root; a
self-check `scan_sees_1.25_on_stated_range`; a negative control — the same scan on
`[1.5, 2.5]` (which excludes the root) must report `contains_1.25 = false`. The
verdict was rewritten: 1.25 is **reachable**, but **only** at `r/Rₛ = 2.78`, which
the paper neither names nor singles out. The claim weakens from “impossible to
produce” to “produced only by an unjustified choice of radius”.

**What did NOT change.** The main H4 result (0.6 and 1.2 ms are not independent
estimates; they differ by the geometric factor 2.008) is untouched: it does not use
the redshift factor.

**Verification.** `h4_paper_arithmetic.py` → exit 0, all four self-checks true.
`verify_h4_independent.py` → `H4_INDEPENDENT_CONFIRMED` (does **not** import the
audited script; solves the equation **symbolically** with sympy, root `25/9`).
Negative control: corrupting the root in the artifact (2.7778 → 3.0) turns the
verifier red. Negative control 2: the archived v1 reproduces the false verdict on
the same input. The original is preserved at
`work/archive/h4_paper_arithmetic_v1_original.py` (sha256 `15f151c1…`).

### 2. Variant B — the arrow of time from boundary conditions

**File:** `work/lean/PastHypothesis.lean`. The logical skeleton of the past
hypothesis is formalized: the arrow does not come from T-symmetric dynamics but
requires a boundary condition that **breaks** T-symmetry.

| theorem | content |
|---|---|
| `rev_isTraj` | the reversal of a trajectory is again a trajectory |
| `arrow_vanishes_on_symmetric` | a self-reversing trajectory carries a zero arrow |
| `rev_allowed_of_invariant` | a T-invariant boundary preserves closure |
| `arrow_sum_zero` | **conclusion:** on a closed family the arrow sums to 0 |

**Non-vacuity (mandatory).** `demo_arrow_vanishes`; `paramArrow_isArrow` +
`paramArrow_ne_zero` (a non-zero arrow really exists); `demo_closure_is_the_line` —
the **decisive witness**: the same arrow on two families gives the sum `0` (closed)
versus `2` (not closed), so the closure hypothesis is exactly the line that
separates the conclusion; `Bfuture_not_invariant` +
`demo_boundary_breaks_reversal` — the boundary `{n | 0 ≤ n}` is not T-invariant.

**Verification.** `sh verify_lean.sh` → rc=0, `NO_SORRY_NO_AXIOM`; the axiom audit
of all 12 theorems yields only `[propext, Classical.choice, Quot.sound]`; a scan for
`native_decide`/`unsafe`/`implemented_by` is empty. **Negative controls:** removing
`TRSymmetric` from `rev_isTraj` → compile error; removing closure from
`arrow_sum_zero` → 4 errors; inserting a `sorry` → the scan catches it.

---

## Round 219 — H3 on the real GW150914 data, H4 against the paper’s text

Owner input: h5py 3.16.0 is installed, both GW150914 files are readable (524 288
samples, 32 s), check H3 on the data and H4 against the source materials, bring it
to a reproducible result, apply independent and negative controls, and do not claim
there are no gaps if something is untested — list what is confirmed, what is
refuted, and what limitations remain.

### 1. H3 on real data — CONFIRMED on a measured PSD

`work/h3_real_data.py` → exit 0, **15 of 15 thresholds**. Artifact
`work/artifacts/h3_real_data.json`.

| quantity | H1 | L1 |
|---|---:|---:|
| measured ASD(100 Hz) | 1.03·10⁻²³ | 9.85·10⁻²⁴ /√Hz |
| match-filter peak | 7.27 | 5.58 |
| peak at GPS | 1126259462.42 | 1126259462.42 |
| offset from the published merger | +19.6 ms | +19.2 ms |
| H1 − L1 | **0.37 ms** | |
| fraction of `τ₀` absorbed by the source | 1 − 8.9·10⁻¹⁶ | 1 − 1.1·10⁻¹⁵ |
| quadratic-term SNR at 1.2 ms | **0.174** | **0.127** |
| detectable amplitude at SNR 5 | 34.5 ms | 47.3 ms |

**Filter calibration — by measurement, not algebra.** Pure synthetic noise (64
realizations) gives a peak ρ of 4.0 on average (max 4.8) — the Rayleigh maximum for
this band. An injection with target SNR 20 is returned as **20.1**, with timing to
0.12 ms. The same filter on an off-source stretch gives 3.4/3.6.

**Why the peak is 7.27, not 20.** The published single-detector SNR of GW150914 ≈ 20
refers to the full IMR template. Here the template is the leading-order inspiral
only, truncated at 300 Hz, without merger and ringdown. Named, not hidden.

**Independent check.** `verify_h3_eigh.py`: own PSD estimator, delay directions from
the analytic expansion in `x = (πM_c f)^{2/3}`, projection through the eigenvalues
of the Gram matrix instead of QR. Agreement: absorbed fractions 0.88487345 and
0.87164874, residual SNRs to 1e-6.

### 2. H4 against the paper’s text — one of my claims is REFUTED

`work/h4_paper_text_audit.py` → exit 0.

**H4a stands.** The paper’s two “independent” numbers (0.6 and 1.2 ms) differ by
exactly the geometric factor `(ΔV_total)^{1/3}/Rₛ(62) = 2.0083`; the paper’s own
ratio is 2.0. The “agreement” between them is arithmetic.

**H4b revised, the claim weakened.** Direct reading of the text: `0.0006 s` → “of
order 1 millisecond”, `0.0012 s` → `0.0015 s`. Only **one** exact multiplier is
claimed: 1.5/1.2 = 1.25. The apparent second (1.0/0.6 = 1.667) is an artefact of the
rounded phrase, not a derived number.

**H4c refuted — my own error of round 217.** I asserted that 1.25 is reachable only
by the unjustified radius `r = 2.78 Rₛ`. **Wrong.** The volume-weighted mean of
`(1+z)` over the paper’s own stated shell `[1.5, 4] Rₛ` is **1.25034** — the paper’s
1.25 to 0.03 %. The natural weighting (created time ∝ created volume) reproduces the
paper’s underived number from its own geometry. The root `r = 25/9 Rₛ` is a red
herring: there the *pointwise* `1+z` equals 1.25, not the volume mean.

Independent check: `verify_independent_v3.py` recomputes this **symbolically**
(sympy, substitution `r = Rₛ₆₂·u`), agreement 1e-6.

**What remains.** H4a; and that 1.25 is introduced in the text with no derivation —
the reader is given a number and told only “with the additional gravitational
redshift”.

### 3. Six errors of this round — all caught by controls

| # | Error | How it showed up | Caught by |
|---|---|---|---|
| 1 | normalization `C = 4·df` (continuous FT) with an inverse FFT | SNR **4.6·10⁴** on pure noise | the “pure noise ≈ 1” control |
| 2 | template amplitude without the factor `fs` | template optimal SNR **1.7·10⁻⁷** | comparison with the analytic σ |
| 3 | no window | peak at the **segment edge** (ρ = 306) | the peak-position check |
| 4 | synthetic noise not filtered as the real data | noise peak 8–13 instead of 4.8 | the noise-peak distribution |
| 5 | absorbed fraction as `1 − e_res/e_tot` instead of `sqrt(...)` | verifier gave 0.783 against 0.885 | the verifier’s cross-check |
| 6 | the verifier gave each mass its **own** shell | factor **3.506** against 2.008 | comparison with the artifact |

The sixth is the most useful: the independent verifier disagreed loudly on a
quantity I had asserted since round 216, and it turned out the **verifier** was
wrong — the paper (4.3) evaluates all three masses on **one** shell
`[1.5 Rₛ(62), 4 Rₛ(62)]`. Fixed; 2.0083 then confirmed by a second, symbolic route.

**The PSD estimator changed its mind twice, both times by measurement.** Median
Welch gives ASD(100 Hz) = 1.03·10⁻²³ — the published aLIGO O1 sensitivity. The mean
gives 2.3·10⁻²², 22× higher: the release contains non-Gaussian glitches to which the
mean is not robust. The median was kept; both numbers are in the code.

### 4. Negative controls on the verifier itself

Eight single-field corruptions of the artifact, each in isolation — **all eight go
red**: the H3 peak, the residual SNR, the noise maximum, ASD(100 Hz), the geometric
factor, the volume mean, the 1.25 radius, the NC3 flag.

### 5. What is confirmed, refuted, and still open

**Confirmed (measured on real data).** H1 and H2: constant and linear delays are
absorbed by the source parameters to 10⁻¹⁴. H3: the quadratic term at 1.2 ms gives
an SNR of 0.13–0.17 on the measured PSD. H4a: the paper’s two numbers are one
number and the geometric factor 2.0083.

**Refuted.** My round-217 claim that 1.25 is reachable only by an unjustified
radius: the volume mean over the stated shell gives 1.25034. My round-217 framing
of “two different multipliers”: there is only one exact multiplier.

**Limitations that remain.**
1. The template is a leading-order inspiral, not IMR. H1/H2 do not depend on it
   (algebraically exact); the numerical value of H3 does.
2. The mirror copy of the data is not GWOSC-authenticated (the host cannot reach
   GWOSC). Only the owner’s host can lift this. — *Later resolved: the owner
   compared hashes against GWOSC directly and they matched; the caveat is
   withdrawn (round 220, and the amendment in `DATA_PROVENANCE.md`).*
3. The H4c result (1.25034) rests on the weight “time ∝ volume”. Other weights give
   1.27–1.30, so 1.25 is agreement under a natural weight, not a uniquely possible
   number.
4. `run_command` refuses commands containing the name of the protected home; all
   operations on `data/` went through `read_file`/`h5py`.

---

## Round 220 — higher PN, ET/CE, direct fit, homogeneous universe, four remarks

**Owner input:** (1) withdraw the GWOSC caveat (mirror and official hashes matched
one for one); (2) the next step (a): higher PN + an ET/CE estimate; (3) variant (c):
homogeneous universe as a theorem; (4) four remarks on the previous report.

### 0. The GWOSC caveat — withdrawn by measurement

The owner downloaded both files directly from GWOSC and compared hashes — they
matched. I confirmed the local bytes independently (`file_fingerprint`):
`81040e1e…dd1b97` (4 068 797 b) and `23a20702…e41b1f` (3 919 118 b). The caveat
“the mirror is not GWOSC-authenticated” is **withdrawn** in `REPORT.md` §7 and in
`data/DATA_PROVENANCE.md` (a dated amendment). The owner is right: it was a legal,
not a scientific, question; the correctness of the data already rested on three
independent physical checks.

### 1. Step (a), part 1: the degeneracy is not a leading-order artefact

Two **exact** identities hold at **any** PN order:

```
(I)  ∂Ψ/∂t_c = 2π f
(II) f ∂Ψ_chirp/∂f = M_c ∂Ψ_chirp/∂M_c      (Euler in v)
```

`(II)` holds because `Ψ_chirp` depends on `M_c` and `f` **only through**
`v = (πM_c η^{-3/5} f)^{1/3}`; no PN coefficient enters. Hence
`ΔΨ_lin = τ̇[M_c ∂Ψ/∂M_c + (t_c − t_ref) ∂Ψ/∂t_c]` — exactly a combination of the
source directions.

| basis | const | linear | quadratic SNR |
|---|---:|---:|---:|
| leading order | 1 − 6.7·10⁻¹⁶ | 1 − 2.2·10⁻¹⁶ | 0.0872 |
| 3.5PN `(M_c,t_c,φ_c)` | 1 − 1.0·10⁻¹⁵ | 1 − 5.6·10⁻¹⁶ | 0.0551 |
| 3.5PN `(M_c,η,t_c,φ_c)` | 1 − 1.0·10⁻¹⁵ | 1 − 5.6·10⁻¹⁶ | 0.0551 |

Identity `(II)` was verified **symbolically** with a negative control: adding
`ln(M_c/f)` breaks it, while `ln(M_c·f)` does **not** (because `M_c f` is a function
of `v`) — checked so the control is not vacuous. `pn_orders.py` exit 0.

**A property of the metric that had to be named.** The absorbed fraction is **not**
invariant under the choice of `t_ref` (0.9999967 versus 0.9530), while the residual
SNR is invariant to 2.3·10⁻¹⁰. So the absorbed fraction must **never** be quoted as
evidence — only the residual SNR. Entered as the threshold
`H_absorbed_fraction_NOT_reference_invariant`.

### 2. Step (a), part 2: ET/CE — and the honest ceiling

No real ET-D/CE curve exists in the container. Quoted: ET aims at “~factor of 10” in
a band from a few Hz (`src_e8a31a47672c`). Hence this is a study of **statistics**,
not a forecast:

| convention | 20 Hz aLIGO | 20 Hz ×10 | 5 Hz aLIGO | 5 Hz ×20 |
|---|---:|---:|---:|---:|
| amplitude tied to 1.2 ms at the band edge | 0.087 | 0.872 | 0.0047 | 0.094 |
| `τ₂` fixed | 0.087 | 0.872 | 7.625 | 152.5 |

The scaling in sensitivity is **exactly linear** (10.000000000 and 20.000000000).
**The sign of the wide-band benefit depends on the convention and flips**, because
the paper **does not specify `τ(t)`**. That is the ceiling of step (a).

### 3. Step (c): homogeneous universe — a theorem

`lean/HomogeneousUniverse.lean`, `obs H δ a = H a + δ a`:

* `creation_not_identifiable_unpinned` — with a **free** base expansion (dark-energy
  fit) `δ` is not identifiable;
* `creation_identifiable_pinned` — with **fixed** `H` it is identifiable;
* `baseline_unique_when_pinned` — the fixing is the line that separates the branches;
* `demo_unpinned_pair` / `demo_pinned_identifiable` — **non-vacuity**.

`sh verify_lean_cosmo.sh` → rc=0, `NO_SORRY_NO_AXIOM`, only the three standard
axioms. **Three negative controls go red:** removing the final step, asserting the
opposite, inserting a `sorry`.

### 4. The owner’s four remarks

**4.1 What the theory predicts — accepted, NOT closed.** The final phrasing was
narrowed: what is proven is degeneracy **within the adopted linear model**, not the
impossibility of testing the whole “now”-theory. New: the degeneracy is not an
artefact of truncation (§1).

**4.2 H3: the 37.3↔7.27 gap decomposed, and the delay fitted directly.**
`h3_direct_fit.py` exit 0.

| component | H1 | L1 | verdict |
|---|---:|---:|---|
| B1 calibration (injection 20) | 20.109 | 20.062 | **excluded** |
| B2 amplitude scale | 0.195 | 0.166 | remains: the data are 5.1× weaker |
| B3 fraction of σ² above ISCO (67.6 Hz) | 0.692 | 0.709 | measured |

**Direct fit of `τ₂`:** the profile at `τ₂ = 0` is 7.2732 / 5.5779; at the paper’s
`τ₂`, 7.2614 / 5.5795; the drop over the ±9.6 ms grid is **0.213 / 0.046** against a
1σ threshold. **The data do not constrain `τ₂` at all.** B4 (template phase
imperfection) cannot be separated from B3 on these data — named.

**4.3 H4: the weight — you are right, confirmed by measurement.**

| weight | mean (1+z) | difference from 1.25 |
|---|---:|---:|
| full proper volume | **1.250339** | **0.03 %** |
| created volume | **1.285956** | 2.9 % |
| flat volume | 1.241729 | 0.66 % |

My docstring **named** the created volume while the code **computed** the full one —
the same “word versus execution” gap as the harness defects. The claim is weakened:
“1.25 is reproduced by natural weighting”, not “derived”.

**4.4 Pre-registration: T8 brought into line.** The frozen T8 is the false-alarm
fraction over 100 realizations. `prereg_t8.py`:

| | value |
|---|---:|
| `z` of the signal (1.2 ms) | 0.0872 |
| σ(`z`) of the noise | 0.933 |
| max\|`z`\| of the noise | 2.13 |
| false alarms | **0 / 100** → **T8 passed** |

**But T8 alone is insufficient:** it is about the detector, not the effect — the
delay gives `z = 0.087`, 57× below the threshold of 5. The Fisher number was renamed
**T8b** and labelled post-hoc. Additionally: the delay column `h·2πf` in
`analysis.py` is **bit-identical** to the `t_c` column, so that Fisher matrix is
exactly rank-deficient, and 1.8·10²⁰ is a float64 roundoff floor (my recomputation
gives 3.1·10²⁰: the value does **not** reproduce, and is therefore not quoted as a
measurement).

**4.5 The stale “no h5py” — fixed.**

### 5. Independent verification of this round

`verify_new_results.py`: the PN identity **symbolically**; the projection via
**normal equations** instead of QR; the H4 weights by **symbolic integration**
instead of quad; the scaling and the profile flatness from the recorded numbers.
Then `--selftest`: corrupting each of the four artifacts — **all four go red**
(`SELFTEST_PASS`).

### 6. Regression

**19 of 20 commands exit 0.** The only non-zero is `mutation_control.py` — exit 1
of a **pre-existing** state (15/20; the file was not touched). `verify_lean.sh`
rc=0 (5 files), `verify_lean_cosmo.sh` rc=0, `verify_new_results.py --selftest` →
SELFTEST_PASS.

### 7. Nine errors of this round — all caught by controls

| # | Error | How it showed up | Caught by |
|---|---|---|---|
| 1 | 3PN log term outside the `v⁶` factor | the 3.5PN phase wrong by ~27×, `t_span` = **32.5 s** | the `t_span_plausible` threshold (0.5–2 s) |
| 2 | the identity via `np.gradient` | a spurious 0.46 % failure of an exact identity | sympy + an analytic derivative |
| 3 | the negative control `ln(M_c·f)` | the control did **not** go red (it is a function of `v`) | the control itself |
| 4 | a leading-order test against a 3.5PN direction | a spurious leak of 4.8·10⁻⁹ | re-deriving per order |
| 5 | `1/√(1−1/r)` with `r` in metres | all three weights = 1.000001 | the weights agreed instead of spreading |
| 6 | `f_isco` mixed units | 2.7·10³⁷ Hz | the fraction of σ² above ISCO came out exactly 0 |
| 7 | the 1σ interval returned `None` | looked like a bug, was the **result** | renamed `unconstrained_on_grid` |
| 8 | symbolic W2 returned a complex number | `TypeError` | asserting `Im < 1e-12` |
| 9 | the noise scale in T8 | 100/100 false alarms | the false-alarm rate itself |

Recorded in `NOTES.md`, section “Round 220”.

---

## Round 221 — four open questions closed by measurement

**Owner input:** “What remains open, stated honestly”, plus a table of four
questions. The owner said plainly what bothered him, and it was a hypothesis about
my error: 37.3 for leading order is higher than the full IMR gives (20).

### 0. The outcome in one table

| question | owner’s verdict | what the measurement showed |
|---|---|---|
| **37.3 vs 20** | closable, priority | **The owner is right in substance.** 37.3 reproduces exactly (not a normalization error), but it is an integral over 20–300 Hz, and **69–71 % of σ² lies above ISCO** (67.63 Hz). Cutting at ISCO gives **20.71 / 18.17** — the published value. Cutting at 250 Hz gives 36.58 / 32.87, so this is a **validity-region** effect, not a bandwidth effect. **But B2 ≠ B3:** if B2 were B3, the scale after the cut would go to 1; instead it goes 0.195 → 0.339 and stops. |
| **ET/CE sign** | no, without τ(t) | untouched: the answer stands, a dichotomy instead of a number. |
| **B4 vs B3** | partially | **separated by measurement.** The band cut explains ~1.7 of the 5.1 and leaves ~3 unexplained. The residue is a real open item. |
| **Grid ±9.6 ms** | closable | **closed two ways.** Expanded fourfold (±38.4 ms, 129 points) — still unconstrained; and a **physical bound of 1.667 ms** from the phase criterion, so the paper’s 1.2 ms takes **72 %** of the budget. |

### 1. Where 37.3 comes from (`h3_rho_opt_bands.py`)

Measure: `rho_opt(band) = sqrt(C · Σ_{k∈band} |H_k|²/S_k)`, the **same discrete
convention** as `h3_real_data.py`, so the number is directly comparable with the
recorded 37.31.

| band | H1 | L1 | fraction of full |
|---|---:|---:|---:|
| 20–300 Hz (as in the pipeline) | **37.3118** | **33.7097** | 1.000 |
| 20–250 Hz | 36.5850 | 32.8747 | 0.981 / 0.975 |
| **20–67.63 Hz (up to ISCO)** | **20.7098** | **18.1707** | **0.555 / 0.539** |

**Controls.** NC1 — the full band reproduces the recorded 37.31184231798413
**exactly** (rel < 1e-6): not a normalization error, the same quantity. NC2 — the
truncated bands are strictly smaller. NC3 — linearity in amplitude (2.000000000).

**Fraction of σ² above ISCO:** 0.6919 (H1), 0.7094 (L1) — matching the earlier
0.692 / 0.709, reproduced independently.

**The decisive check: B2 = B3 or not.** If “the data are 5.1× weaker” were the same
as “the template claims power above ISCO”, then cutting the band at the validity
edge would drive the preferred scale to 1. Measured:

| band | scale H1 | scale L1 |
|---|---:|---:|
| 20–300 Hz | **0.1948** | **0.1657** |
| 20–250 Hz | 0.2183 | 0.1578 |
| 20–67.63 Hz | **0.3385** | **0.2909** |

The deficit shrinks from 5.1× to **3.0× (H1) / 3.4× (L1)** — and **stops**. So B3
explains roughly a factor 1.7 of the 5.1 and leaves a **real factor ~3** that the
band does not explain. That is the B3/B4 separation obtained by measurement, not by
argument.

### 2. The τ₂ grid — expanded and bounded by physics (`h3_grid_bounds.py`)

**Expansion.** The grid ±9.6 ms → **±38.4 ms** (fourfold, 129 points), same filter
construction.

| | H1 | L1 |
|---|---:|---:|
| profile peak | 7.8106 | 5.5795 |
| at τ₂ = 0 | 7.2732 | 5.5779 |
| drop over the grid | **0.791** | **0.461** |
| constrained? | **no** | **no** |

**Controls.** NC1 — the coarse peak at τ₂ = 0 matches the recorded
`profile_at_zero` 7.273249726912924 / 5.577882654444806 to 1e-6. (The first version
compared against the recorded `peak_snr` 7.2739 / 5.5799 and gave a **false red** —
a different quantity, from a two-stage fine search over `t_c`.) NC2 — the profile at
τ₂ = 0 equals the coarse peak.

**Something important that cannot be left in the artifact.** On H1 the profile
maximum sits **at the grid edge** (−38.4 ms), i.e. it is still **rising** when the
scan stops. That is not a best value and not a bound — it means there is no interior
maximum at all. The strongest form of “unconstrained”.

**The physical bound.** The delay enters as the phase `2π f τ(t)`. A delay whose
phase sweeps more than π across the band is no longer a small perturbation of the
template — it is a different waveform. This gives a ceiling **independent of the
scan**:

```
max |2π f · τ₂ · (t(f) − t_ref)²| ≤ π    over the band
```

| | value |
|---|---:|
| τ₂ by the phase criterion | 0.0023370481 /s² |
| **amplitude** | **1.6667 ms** |
| ceiling from segment duration (4 s) | 4000 ms |
| claimed by the paper | 1.2 ms |
| **fraction of the budget used by the paper** | **72.0 %** |

Independent check: the maximum of the kernel is reached at the **upper edge** of the
band (300 Hz), and my recomputation with a different scan (2·10⁶ points) gives the
same τ₂ to 1e-6.

**Meaning.** The paper’s number is not merely unconstrained by the data — it is
**close to the edge of what can be called a small perturbation** of this template.

### 3. Independent verification of this round

`verify_rho_opt_bands.py` does not import the audited script. It builds its **own**
PSD (own Welch code, own log-log interpolant), computes `f_ISCO` **another way**
(total mass in SI rather than the chirp-mass formula), and takes a **different
quadrature** (continuous integral instead of a discrete sum).

| check | result |
|---|---|
| independent discrete sum == recorded | **rel = 0.00e+00 on all six cells** |
| f_ISCO another way | 67.630422 == 67.630422 (rel < 1e-9) |
| continuous vs discrete | within 3e-3 on all bands |
| convergence of the discrete sum as df→0 | 8.4e-4 → 2.8e-5 → 2.9e-6 (monotone) |
| scale after the cut moves by 1.5–2.5× | 1.737 / 1.755 |
| scale == peak/rho_opt | exact |

**Self-test: 4/4 corruptions go red** (`SELFTEST_PASS`).

**Four errors of this round, all caught by controls.**

| # | Error | How it showed up | Caught by |
|---|---|---|---|
| 1 | the verifier interpolated the PSD **linearly**, the audit **log-log** | a discrepancy of 1.7e-3 on **every** band | convergence did not converge; measuring the two interpolants gave 2.8e-3 versus 8.4e-4 |
| 2 | “convergence” on a **moving** band (f_ISCO not a multiple of 0.25 Hz) | a spurious non-monotone 6.3e-4 → 2.0e-4 → 1.4e-4 | the test itself: the band edge was moving |
| 3 | two first-draft thresholds were **wrong**: “ISCO band < 20” (it is 20.71) and “IMR band within 8 of 20” (it is 36.58) | two false red cells | the measured values themselves; the criteria were rewritten and the wrong ones **recorded**, not deleted |
| 4 | `h3_grid_bounds.py` NC1 compared against `peak_snr` instead of `profile_at_zero` | a false red | the comparison itself: 7.2739 versus 7.2732 |

None is in the code — all four are in the **statement of what is being compared**.

### 4. What remains open, stated

1. **The factor-~3 residue** in the 37.3↔7.27 gap after removing B1 and B3.
   Candidates — the event’s response versus a sky average, and the definition of
   “20” as an IMR number. Not separable without an IMR waveform. **The main open
   item.** *(Later closed: rounds 224/226.)*
2. **B4** (template phase imperfection) is separated from B3 **quantitatively**
   (B3 gives 1.7 of the 5.1) but **not named**.
3. **ET/CE** — the sign is still undetermined without τ(t). Untouched.
4. **Deriving τ(t) from the theory** — still not done; the paper does not state it.

---

## Round 224 — the factor-~3 residue closed with a real IMR template

**Owner input:** “Continue the saved unfinished work from the last point. Do not
repeat what is already done.” The last saved point was `msg221-note.md` §4.1: the
**factor-~3 residue in the 37.3↔7.27 gap**, declared a ceiling (“not separable
without an IMR waveform”).

**What was done:** the ceiling is lifted. phenomxpy 2.0.3 is available, so the
inspiral template was replaced by a full `IMRPhenomT`.

### 0. The outcome in one table

| question from §15.2 | before | after |
|---|---|---|
| does the IMR template lower ρ_opt? | unknown | **yes, by 15 %**: 37.31→**31.66** (H1), 33.71→**28.74** (L1) |
| is that enough for ~20? | — | **no.** The threshold “within 25 % of 20” **FAILED** (58 %). Recorded as a failed threshold |
| was the residue a template defect? | hypothesis | **more than half — yes.** Amplitude scale 0.1949→**0.5288** (H1), 0.1655→**0.4470** (L1) |
| what produces the residue? | “not separable” | **orientation.** ρ_opt = 20 at **ι ≈ 59.2°**, inside the physically allowed range |

### 1. Setup (`h3_imr_check.py`)

`IMRPhenomT` (2,2), `M_c = 28.096 M_⊙`, `D_L = 410 Mpc`, spins 0, `f_min = 20 Hz`,
`df = 0.25 Hz` — the **same grid** as the pipeline (`seg_n = 4·FS = 65536`), so the
numbers are directly comparable with the recorded ones. The measure is the **same
discrete convention** as `h3_rho_opt_bands.py`:

```
rho_opt(band)² = 4·df · Σ_{k∈band} |h~(f_k)|² / Sn(f_k)
```

### 2. The numbers

| band | H1 LO | H1 IMR | L1 LO | L1 IMR |
|---|---:|---:|---:|---:|
| 20–300 Hz | 37.3118 | **31.6595** | 33.7097 | **28.7350** |
| 20–250 Hz | 36.5850 | 30.9531 | 32.8747 | 27.9442 |
| 20–67.63 Hz | 20.7098 | 16.7296 | 18.1707 | 14.6924 |

**Preferred amplitude scale** (filter peak / ρ_opt):

| | LO | IMR |
|---|---:|---:|
| H1 | 0.1949 | **0.5288** |
| L1 | 0.1655 | **0.4470** |

The deficit shrinks from 5.1× to **1.9× (H1) / 2.2× (L1)**. So **more than half** of
the earlier residue was a leading-order defect (B4), not the data and not the PSD.
This is the quantitative “naming” of B4 that §15.5 item 2 lacked.

### 3. The residue is orientation

`inclination = 0` (face-on) gives the **maximum** amplitude, while the published
~20 is the number for the event’s actual orientation. Scan (H1):

| ι, ° | 0 | 15 | 30 | 45 | 60 | 90 |
|---|---:|---:|---:|---:|---:|---:|
| ρ_opt | 31.66 | 30.60 | 27.70 | 23.74 | **19.79** | 15.83 |

**ρ_opt = 20 is reached at ι ≈ 59.19°.** Monotone decreasing; the crossing point is
reproduced by the independent verifier (19.9817 at ι = 59.1934°).

**Meaning:** the residue stops being an “unexplained factor ~3” and becomes
“orientation + the definition of the published number”. Stronger than “not
separable without IMR”, weaker than “the theory is refuted”: **the data do not
contradict the model, and the residue is explained by detector physics.**

### 4. Three errors of this round

| # | Error | How it showed up | Caught by |
|---|---|---|---|
| 1 | multiplied the inspiral amplitude by `fs` and did not divide by `fs` in the matched filter | NC1 gave rho_opt = 297708 instead of 37.31 — off by 8000× | **NC1** |
| 2 | NC5 asserted “the model is strictly zero below its own `f_min`” | **false**: phenomxpy carries a tail to 2.6·10⁻²³; measured 0.52 % of band power | the control itself (it did not go red where it should have) |
| 3 | the verifier built the model at `df = 0.5` and interpolated linearly | **19 %** of the band power lost | the 2.2 % disagreement with the audit |

**The third is the most useful.** The FD waveform phase rotates by ~π between
adjacent samples at `df = 0.5 Hz`, so linear interpolation of a complex function is
not interpolation. Measured: at **shared** frequencies the two grids agree to
**1e-16**, while after interpolation the power ratio is **0.810**. So **the verifier
was wrong and the audit was right** — established by measurement, not by argument.
The verifier was rewritten on the native grid; the false NC5 claim and the wrong
1e-3 threshold are **recorded in the code, not deleted**.

### 5. Controls (all pass, exit 0)

| control | result |
|---|---|
| NC1 inspiral rho_opt == recorded | 37.31184231798413 / 33.709713860505396, rel < 1e-6 |
| NC2 IMR rho_opt finite and > 0 in all bands | yes |
| NC3 linearity in amplitude | 2.000000000 |
| NC4 zero waveform → ρ = 0 | 0.0 |
| NC5 power below the model’s `f_min` | **0.52 %** of band (< 1 %) |

**Controls and findings are separated on purpose.** A failed **control** means a
broken instrument and a void run; a failed **finding** is a result. The exit code
gates on the controls only; the threshold `T_IMR_full_band_within_25pct_of_20` is
**expectedly failed**, and that is the answer.

### 6. Independent verification (`verify_h3_imr.py`)

Imports nothing from the audit: **own** Welch, **own** log-log PSD interpolant,
**own** `IMRPhenomT` call, **own** quadrature.

| check | result |
|---|---|
| independent IMR rho_opt == recorded | **rel = 0.00e+00** (H1 and L1) |
| inspiral rho_opt == recorded 37.31 | yes |
| inclination scan monotone | yes |
| crossing point ρ = 20 | 59.19340049167065 == recorded |
| ρ at that point (independently) | **19.9817** |
| internal identity `scale == peak/ρ_opt` | yes |
| **self-test: 5 corruptions of the artifact** | **5/5 go red** (`SELFTEST_PASS`) |

The artifact reproduces **byte-for-byte** on a re-run: sha256 `36ed9094…`.

### 7. What remains open, stated

1. **Spins = 0** and the base number at `ι = 0` are simplifications named in the
   artifact; the IMR ρ_opt here is an **upper bound**, not the event’s recovered
   SNR. *(Round 226 closed this with spins over the published envelope.)*
2. **The delay is still not fitted** to the strain: only the template’s optimal SNR
   is measured here. *(Round 226 closed this.)*
3. **ET/CE** — unchanged: the sign is undetermined without `τ(t)`.
4. **Deriving `τ(t)` from the theory** — not done; the paper does not state it.

---

## Round 226 — both items of §16.5 closed

**Owner input:** (1) run IMRPhenomT **with spins**, so that “upper bound” becomes
“recovered SNR”; (2) fit the delay `τ₂` **to the strain** with the full IMR template
instead of the inspiral — “this is what the paper needs”.

### 0. The outcome in one table

| question | before (§16.5) | after |
|---|---|---|
| spins | “upper bound” at ι=0, spin 0 | **range over the published envelope**: 30.02–32.49 (H1), 27.16–29.54 (L1) |
| is the published SNR reached? | untested | **not by spins**; by inclination yes for H1 (ι ≈ 59.8°), no for L1 (min 13.58 at ι = 90°) |
| fitting τ₂ to the strain | inspiral, “unconstrained” | **full IMR: a one-sided bound appears** (drop 2.93 H1 / 2.14 L1) |
| power control | absent | **present**: 20 ms is recovered, 1.2 ms is not |

### 1. The comparison target was corrected

Since round 216 the project compared against “the published ~20”. No such number
exists: **19.5 (H1)** and **13.3 (L1)** are the re-weighted single-detector SNRs
(arXiv:1602.03839); 24 (arXiv:1602.03837) and 25.1 (arXiv:1602.03840) are
**network** values. Comparing a single-detector ρ_opt to a network value is a
category error.

### 2. Item 1 — `h3_imr_spins.py`

The spin envelope is **not a rectangle**: a1 and a2 are linked through
χ_eff = −0.07 ± 0.17. The corner (0.69, 0.89) carries χ_eff = +0.78 and is excluded
by measurement. Scan: `|a1| ≤ 0.69, |a2| ≤ 0.89, χ_eff ∈ [−0.24, +0.09]`, 58 points.

| set | χ_eff | ρ_opt (H1) | /ρ_opt(0) |
|---|---:|---:|---:|
| zero | 0.000 | 31.6595 | 1.0000 |
| χ_eff median | −0.070 | 31.1073 | 0.9826 |
| medians | +0.373 | 34.8362 | 1.1003 |
| bounds (outside χ_eff) | +0.778 | 38.5973 | 1.2191 |
| bounds_low (outside χ_eff) | −0.778 | 26.5692 | 0.8392 |

**Spins do not reach the published SNR**, and the largest spin effect moves the
number **up**. Inclination explains H1: 19.5 at **ι ≈ 59.8°**. The inclination
multiplier is **(1+cos²ι)/2 exactly** and independent of spins (3.3·10⁻¹⁶) — so the
joint range is exact, not a sampled approximation.

The published SNR ratio 1.466 versus 1.102 (optimal) and 1.303 (filter-extracted):
the residue is the **antenna response**, which the model lacks.

### 3. Item 2 — `h3_direct_fit_imr.py`

**A decision, not cosmetics:** `τ(t)` needs the time at each frequency. The inspiral
formula runs to zero (`t_PN(250 Hz) = −1.0 ms`), while the IMR **plateaus**
(`t_IMR(250 Hz) = −155 ms`). The ratio goes 1.03 → **155**. Both conventions were
computed.

| template / kernel / grid | drop H1 | drop L1 | bound |
|---|---:|---:|---|
| LO / pn / ±9.6 ms | 0.213 | 0.046 | no |
| LO / pn / ±38.4 ms | 0.791 | 0.461 | no |
| IMR / pn / ±9.6 ms | 0.668 | 0.506 | no |
| **IMR / pn / ±38.4 ms** | **2.925** | **2.144** | **yes** |
| **IMR / imr / ±38.4 ms** | **1.371** | **1.047** | **yes** |

**The answer to the referee splits by template — that is the finding.** With the
leading inspiral there is no constraint on any grid. With the full IMR on the
±38.4 ms grid the profile drops by more than 1σ in one direction: a **one-sided
bound** of ≈ −10…−13 ms (H1) and +19…+37 ms (L1). The first time the project
extracted any constraint on `τ₂` from the data rather than from a plausibility
argument.

**Power control:** the same fit recovers a 20 ms injection (21.75 / 20.17) and does
**not** recover 1.2 ms (0.15 / −0.56 at a spread of 4–5 ms).

### 4. My five errors, all caught by a control

1. `irfft` instead of `ifft` → the peak underestimated by 16384× (**NC6**).
2. `B = zeros(seg_n)` with an rfft-length mask → a silent length mismatch (**NC6**).
3. `template_time_series` divided by `fs` → the round trip gave σ/fs (**NC7**).
4. A ±9.6 ms grid for a 20 ms injection → the grid edge as the answer (**NC4**).
5. Verifier: a two-stage peak search **silently missed the maximum** (7.8106 instead
   of 8.2011) → deleted.

**A false criterion recorded, not deleted:** NC6 of the first draft (“more power →
not a lower peak”) is **false** — the spins raise ρ_opt 31.66 → 34.84 while the peak
falls 16.74 → 13.51, because the phase changes. Replaced by Cauchy–Schwarz.

### 5. Independent verification — `verify_h3_round226.py`

Own Welch, own PSD, own IMRPhenomT, own quadrature (explicit matrix), own phase
scan. **All eight checks rel = 0.00e+00** (V1, V2, V4, V5, V6) or 1.15·10⁻¹⁴ (V5b,
V7); V3 = 4.44·10⁻¹⁶. **Self-test: 7/7 corruptions go red.**

**NC1b** — LO+pn+expanded reproduces the recorded `h3_grid_bounds.json` to rel =
0.0. Without it the new numbers would not be comparable with the old ones.

Both artifacts reproduce **byte-for-byte** (`d967c59e…`, `53b40ebc…`). Regression:
**20/20 scripts exit 0**, Lean rc=0 in both files, `mutation_control.py` exit 1 —
pre-existing (15/20, file untouched).

### 6. What remains open

1. **L1 does not reach 13.3** even at ι = 90° (min 13.58) — antenna response.
2. **`τ(t)` is not derived from the theory** — the paper does not state it. The main
   open item.
3. Aligned spins: `χ_p < 0.71` is not covered.
4. Phase-only delay; the amplitude is not delayed.
5. A joint fit over (M_c, η, spins, ι) would widen the profile.

---

## Round 227 — joint H1+L1 fit with antenna responses

See `NOTES.md`, “Round 227”, for the full account. Summary: the four owner items of
msg_00227 were done (joint fit first and bound second; antenna responses from the
published sky position instead of an inclination scan; §1/§7 rewritten and the
chronology moved; a separate claim ledger `CLAIMS.md`). Two structural facts
emerged — the detector response collapses to one complex coefficient
`C_d = F₊A − iF×B`, and `φ_c` is therefore not identifiable in a coherent fit with a
fixed sky position. Six errors, all caught by controls; the most useful was the
hand-written HEALPix implementation, which was **deleted** rather than debugged.

---

## Round 228 — ledger revisions per the owner’s review

The owner found six discrepancies between `CLAIMS.md`/`REPORT.md` and the facts —
none of them physics, all of them **formulation**. The first was found from the
ledger’s **legend**, not from the numbers. All six are fixed; the fourth category
**[DERIVED symbolically]** was introduced; B1 became an upper estimate; “+3.6 ms”
was defined; “0.6 %” was softened; the fit resolution ≈ 5 ms was stated; and the
§9 ↔ E7 contradiction was removed. Two verifier errors were caught by the criterion
going red on correct data — the same “comparing something other than what I think”
class. Full account in `NOTES.md`, “Round 228”, and `CLAIMS.md` §G.

---

## Round 229 — the publication package, and the convention finding

**Owner input (msg_00229):** produce a PDF paper, translate the Russian artifacts
into English (Russian kept, English primary), state the authorship plainly, and say
explicitly that the agent did the work and walked the whole path.

**Produced:** a self-contained `publication/` package — the paper in English
(primary) and Russian, both Markdown and PDF; English translations of the report,
claim ledger, pre-registration, error log and round summaries; README, artifact
index with hashes, manifest, `CITATION.cff`, `.zenodo.json`, licences; and an
independent verifier (`verify_publication.py`, 156 checks, 0 failures) with a
self-test (`selftest_publication.py`, 10/10 corruptions caught). The Russian
originals are untouched.

**The substantive finding — the π-bound depends on the convention by 15×.** An
independent re-derivation of the phase bound returned **24.9635 ms** where the
artifact gives **1.6667 ms**. The cause is not arithmetic but the time-reference
convention: the artifact uses the paper’s picture (delay zero at the band start,
growing toward merger), which gives 1.667 ms and makes the paper’s 1.2 ms 72 % of
the budget; placing the vertex at coalescence (delay zero at merger, backwards
relative to the paper) gives 24.96 ms under the same criterion. The artifact is
right **for the paper’s convention**, and the paper now states the bound as
convention-dependent, with both numbers given. This is the same question that round
226 answered substantively for the group delay (`t_IMR/t_PN` from 1.03 to 155);
again, without a named convention the bound is undetermined by a factor of 15.

**A hole the self-test found in its own verifier:** the first corruption case (C1)
initially *survived*, because the corrupted token occurred twice in the paper. The
check was testing presence, not placement — fixed by using a token that occurs
exactly once, with the limitation documented rather than hidden.
