# NOTES — what went wrong on the way (v1 → round 228)

*English translation of `work/NOTES.md`. The Russian original is kept alongside;
the English text is the primary one for publication.*

This file records my own errors, in order, because the report is only worth
reading if the path to it is visible. Every entry was found by a check, not by
re-reading my own code.

## v1 of `analysis.py` — four errors, all caught by looking at the numbers

| # | Error | How it showed up | Fix |
|---|---|---|---|
| 1 | `Rs62_m = 2*G*MF/C**2` with `MF` already in **seconds** (geometric units) | `Rs62_km = 4.5e-34` — a black hole smaller than a proton | `Rs = 2M` in geometric units → 183 km |
| 2 | Residual computed as a **phase-space** norm, not a strain norm | `unmodelled_snr = 9.7e32` — an SNR larger than the number of atoms in the Sun | multiply by the strain amplitude `h(f)` |
| 3 | `np.linalg.lstsq` on raw columns spanning ~40 orders of magnitude | `max_abs_residual = 2e8` for a direction that is *exactly* in the span | orthonormalise via SVD first |
| 4 | “False-alarm rate” test compared a **phase-space** SNR to a threshold of 1 | 100 % false alarms on pure noise | deleted; replaced by a deterministic flat-likelihood test |

## v2 — the SNR was still wrong

The projection was correct, but the statistic was the norm of a *phase* residual.
A phase residual is not an SNR. Rebuilt around `u = h · ΔΨ` so that
`‖(I−P)u‖_weighted` is a genuine matched-filter SNR, directly comparable to the
detection threshold ≈ 5.

## v3 — the PSD was not calibrated

The analytic PSD gave an in-band signal SNR of **129** for GW150914, whereas the
published single-detector value is **20**. Any SNR computed against it was
inflated by 6.5×. Fixed by calibrating the PSD with one constant so the injected
signal’s SNR equals 20. All headline SNRs are now calibrated.

## v4 — three broken controls

| Control | What was wrong |
|---|---|
| K9 | the “free” parameters absorbed the injected signal, so the residual was flat for the *wrong* reason — it looked like a pass and was a failure |
| K3 | threshold `> 5` was arbitrary; at `τ̇ = 0.01` the non-absorbable part is only SNR 1.8 |
| K8 | the noise was scaled so that its own SNR was ~63, not 20 |

Fixed: K9 now fits only `(t_c, φ_c)` with `M_c` fixed; K3 is reported at two
values of `τ̇`; the noise is normalised to the published SNR.

## v5 — the mutation control found 8 holes in my own test

The automated campaign (`mut_dcd2ff7529de`) returned a **null result**: its
coverage tracer reported “0 executed lines” and it classified all 20 mutants as
unexecuted without running one. That is an instrument failure, not evidence, so I
wrote `mutation_control.py` by hand.

First pass: **11/20 killed, 8 survived.** The survivors were real holes — the
ratio-based thresholds were invariant under exactly the faults that mattered:

| Survivor | Why it survived | Closed by |
|---|---|---|
| M4 `t(f)` coefficient 5/256 → 5/128 | rescales `t(f)`; the projection is span-invariant | S6 (exact algebraic identity) |
| M5 strain exponent −7/6 → −7/5 | rescales `h(f)`; span-invariant | S3 (independent amplitude route) |
| M8 PSD calibration dropped | rescales `w(f)`; span-invariant | S2 (SNR must equal 20) |
| M14 `Rs = 2MC` → `MC` | `Rs62_km` is reported but never asserted | S1 (pins 183.148 km) |
| M16 `D_L` 410 → 41 | calibration absorbs it | S2 + S3 |
| M3 `dΨ/dM_c` factor 5/3 → 3/5 | rescales a basis column | *genuinely equivalent* |
| M6 chirp coefficient 3/128 → 3/64 | rescales `Ψ_chirp` and the delay together | *genuinely equivalent* |
| M19 sign flip in `ΔΨ` | `absorbed_fraction` uses `‖P u‖/‖u‖`, sign-blind | *genuinely equivalent* |

Second pass after adding S1–S6: **15/20 killed, 4 survived.**

The four survivors were then **classified by measurement, not by assertion**
(`equivalence_check.py` runs each mutant and diffs the full artifact against the
baseline, with a tolerance-aware comparison and a two-tier split between headline
claims and scale-dependent diagnostics):

* **M3, M6, M16, M19 — equivalent on every headline claim.** `absorbed_fraction`,
  `argmin_candidate`, `K8_residual_spread` and the geometry numbers are identical
  to the baseline to within 1e-12 absolute / 1e-9 relative. M19 (sign flip) is
  identical in *every* leaf.

**One honest limitation this exposed.** The Fisher condition number and its
singular values are **not** invariant under these faults: M16 moves the condition
number from 1.8e20 to 4.5e26 while every headline claim stands. That is expected —
the condition number is a property of the parameter basis, and rescaling a basis
column changes it. The consequence for the report is concrete: **the condition
number is reported as an order of magnitude (> 1e12, i.e. numerically singular),
not as a measurement.** The same caveat applies to
`K9_over_K8_sharpness_ratio`, which is a ratio of two residuals and therefore
inherits the scaling of the K9 basis direction.

## v2 of `h4_paper_arithmetic.py` — a sampling artifact, found by the owner

The owner read the report and caught a concrete error in §5/H4: I wrote that the
factor 1.25 is unreachable “at any radius in the stated range 1.5–4 Rs”. **That
was false.** The code evaluated the redshift factor `1+z = 1/sqrt(1 - Rs/r)` at
only four radii `{1.5, 2, 3, 4}` and tested `|1+z - 1.25| < 0.02`, which those
four samples miss. The analytic root of `1+z(r) = 1.25` is `r/Rs = 2.7778`, and it
lies **inside** `[1.5, 4]`.

The sharpest part: **v1 computed `r_needed = 2.7778` itself and printed the range
`[1.5, 4]` in the same artifact, but never compared them.** Two numbers sat next
to each other in one JSON file and the verdict contradicted both. This is the same
failure mode as the harness turns: a claim that no check ever executed.

Fixed in v2: the scan covers the whole interval (20001 points), locates the root,
and carries a negative control on `[1.5, 2.5]` (which excludes the root) that must
report `contains_1.25 = false`. `verify_h4_independent.py` recomputes the root
**symbolically** with sympy (the audited script does it by hand) and disagrees if
the artifact’s root is wrong — verified by corrupting the artifact, which turns
the verifier red. The original v1 is archived at
`work/archive/h4_paper_arithmetic_v1_original.py` (sha256 `15f151c1…`) and
reproduces the false verdict on the same input.

**What did NOT change:** the main H4 result — 0.6 ms and 1.2 ms are not
independent estimates — stands, because it never used the redshift factor.

## What the independent verifier caught

`verify_independent.py` found a **sign error in my own assertion** about the `τ̇`
direction: I wrote that the `τ̇` term is `−M_c · dΨ/dM_c`, but sympy gives
`+M_c · dΨ/dM_c`. The algebra in the analysis was right; the check text was wrong.
That is precisely the failure mode an independent path exists to catch.

---

# Round 219 — H3 on the real strain, and H4 against the paper’s own text

The owner said h5py was installed, both GW150914 files were readable, and asked
for H3 to be checked on H1/L1 and H4 against the source materials. Six errors were
made and caught on the way; five of them were caught by controls I had already
written, and one by an independent verifier.

## H3 on real data — five errors, all caught by controls

| # | Error | How it showed up | Caught by |
|---|---|---|---|
| 1 | Matched-filter normalization `C = 4*df` (the continuous-FT form) with an inverse-FFT reconstruction | **noise-only peak SNR ≈ 4.6·10⁴** — an SNR no detector can produce | the noise-only control (must be ~1) |
| 2 | Template amplitude convention: used `h~(f)` directly as the rfft, missing the factor `fs` | template optimal SNR **1.7·10⁻⁷** instead of ~45 | the analytic σ comparison |
| 3 | No window on the segment | peak ran to the **segment edge** (ρ = 306) because the rectangular edges leak into the band | the check that the peak lands near the merger |
| 4 | Synthetic validation noise generated from the full PSD but **not** high-passed, while the real data was | noise-only peak 8–13 instead of ~4.8 | the calibration’s own noise-only distribution |
| 5 | Verifier computed the absorbed fraction as `1 − e_res/e_tot` (absorbed **energy**) instead of `sqrt(1 − e_res/e_tot)` (= ‖Pu‖/‖u‖) | verifier said 0.783 where the artifact said 0.885 | the verifier’s own cross-check |

**What settled each one was a measurement, not algebra.** `mf_norm_probe.py`
settled the normalization by generating noise from a known PSD and reading off
which candidate gives ρ ≈ 1 and recovers an SNR-20 injection as 20; the answer was
`H_k = fs·h~(f_k)` with `C = 4·df/fs²`. The continuous-FT form `4·df` is correct
only for `X_cont = dt·X_rfft`, so it is off by `fs²`.

**The PSD estimator was also settled by measurement, and I changed my mind
twice.** `psd_whiten_check.py` compares the FFT-binned in-band variance with the
PSD integral. With a median estimator the real off-source data carries ~6× more
in-band power than the PSD predicts, while a synthetic series generated from that
PSD reproduces the integral to 0.98 — so the excess is in the data (non-Gaussian
transients in the released strain), not estimator bias. I first switched to the
mean estimator on that basis; then the mean gave ASD(100 Hz) = 2.3·10⁻²² /rtHz,
**22× the published aLIGO O1 sensitivity**, because the mean is not robust to those
same transients. Reverted to the median, which gives 1.03·10⁻²³ — the published
value. The final choice is stated in the code with both numbers.

## H4 against the paper’s text — one of my own earlier claims is refuted

Reading the two source paragraphs directly (not from memory) showed that the
paper’s chain is: bare 0.6 ms → “of order 1 millisecond”, and bare 1.2 ms →
1.5 ms. Two things follow:

* **H4b is REVISED, and it weakens the criticism.** Only **one** precise factor is
  claimed (1.5/1.2 = 1.25). The apparent second factor, 1.0/0.6 = 1.667, is an
  artefact of the paper writing “of order 1 millisecond” — a rounded phrase.
  Treating it as a second derived factor was over-reading the text.
* **H4c is REFUTED, and it is my own error from turn 217.** I claimed the factor
  1.25 was obtainable “only by an unjustified choice of radius” r = 2.78 Rs. That
  is wrong. The **volume-weighted mean of (1+z) over the paper’s own stated shell
  [1.5, 4] Rs is 1.25034** — the paper’s 1.25 to 0.03 %. The natural weighting
  (time created ∝ volume created) reproduces the paper’s unexplained number from
  the paper’s own stated geometry. The `r = 25/9 Rs` root is a red herring: it is
  where the *pointwise* redshift equals 1.25, which is not what a volume-weighted
  average gives.

  This was found by writing the new audit, not by re-reading my old one — and it
  is recorded as a correction of my own earlier claim, with the old claim quoted
  in the artifact next to its refutation.

**What survives H4:** H4a. The paper’s two BARE numbers are not independent
estimates — their ratio is the geometric factor (ΔV_total)^{1/3}/Rs(62) = 2.0083,
and the paper’s own bare ratio is 2.0. So “the two estimates agree” is arithmetic.
And the 1.25, while now explicable, is still introduced with no derivation in the
text: the reader is given a number and told only “with the additional
gravitational redshift”.

## One more error, this time in the new independent verifier

`verify_independent_v3.py` initially reported the geometric factor as **3.506**
against the artifact’s 2.008. The verifier was wrong: it gave each of the three
masses its **own** shell `[1.5 Rs(M), 4 Rs(M)]`, whereas the paper’s eq. (4.3)
evaluates all three on the **same** physical shell `[1.5 Rs(62), 4 Rs(62)]` — that
is the whole point of subtracting the smaller holes’ volume excess from the region
the 62 M☉ hole now occupies. Fixed, and the 2.0083 was then confirmed by a second,
symbolic route (`r = Rs62·u` substitution, `k = 29/62` and `36/62`).

**This is the useful kind of disagreement.** The verifier failing loudly on a
quantity I had asserted since turn 216 is exactly what it is for; had it agreed, I
would have learned nothing about which of the two was right.

## Negative controls on the verifier itself

Eight single-field corruptions of the artifacts, each applied in isolation: every
one turned `verify_independent_v3.py` red (H3 peak SNR, H3 residual SNR,
noise-only maximum, ASD(100 Hz), the geometric factor, the volume-weighted mean,
the 1.25 radius, and the NC3 shell-dependence flag). A verifier that cannot fail
proves nothing, so this was checked rather than asserted.

---

# Round 220 — higher PN orders, ET/CE, direct fit, and the owner’s four review points

The owner accepted the GWOSC provenance point (hashes matched, so the caveat is
withdrawn) and asked for the next step: (a) higher PN orders plus an ET/CE
estimate, and optionally (v) the homogeneous-universe statement as a theorem. The
same message carried four review points on my earlier report.

## Errors of this round, all caught by controls or by re-reading the numbers

| # | Error | How it showed up | Caught by |
|---|---|---|---|
| 1 | The 3PN log term `−(6848/21)·ln(4v)` was subtracted **outside** the `v⁶` factor instead of multiplying it | the 3.5PN phase at 20 Hz was wrong by ~27×, and the reported band time span came out **32.5 s** instead of 0.84 s | the `t_span_plausible` self-check (0.5–2 s) |
| 2 | The identity check `f·dΨ/df = M_c·dΨ/dM_c` used `np.gradient` for the derivative | a spurious 0.46 % “failure” of an identity that is exact | a sympy route + an analytic derivative through `v` |
| 3 | The symbolic **negative control** used `ln(M_c·f)` | the control did **not** break the identity, because `M_c·f = v³/(πη^{-3/5})` IS a function of `v` alone — so the control was vacuous | the control itself failed to go red, which is why it was noticed |
| 4 | A leading-order test compared a 3.5PN delay direction against a leading-order basis | a spurious `4.8·10⁻⁹` leak where the identity requires machine zero | re-deriving per-order (each order tested in its own basis) |
| 5 | `h4_weight_audit.py` applied `1/sqrt(1 − 1/r)` with `r` in **metres** | all three weights came out `1.000001` — the `Rs/r` term underflowed to zero | the three weights were then *identical*, which contradicted the expected spread |
| 6 | `f_isco` mixed geometric and SI units (`c³/(6√6 π G M)`) | ISCO frequency `2.7·10³⁷ Hz` | the B3 fraction-of-σ² above ISCO came out exactly 0 |
| 7 | The direct-fit 1σ interval returned `None` when the profile never dropped by 1 | looked like a bug; it is the **result** — the data do not constrain `τ₂` at all | reported as `unconstrained_on_grid` instead of `None` |
| 8 | `verify_new_results.py`’s symbolic W2 mean returned a complex number | `TypeError: Cannot convert complex to float` | asserted the imaginary part < 1e-12 instead of silently dropping it |
| 9 | `prereg_t8.py` scaled the noise so the residual norm came out ~7·10²³ | 100/100 false alarms — an obviously broken statistic | the false-alarm rate itself (should be ~0) |

## What the owner’s four review points changed

* **Point 1 (what the theory predicts) — accepted, not closed.** The honest ceiling
  is now written into §14.1: what is proved is degeneracy *in the adopted linear
  model*, not that the whole “now” theory is untestable. What this round *adds* is
  that the degeneracy is not a leading-order artefact (§12.1).
* **Point 2 (H3 does not fit the delay to the signal) — done.** `h3_direct_fit.py`
  fits `τ₂` directly and decomposes the 37.3 ↔ 7.27 gap. Result: B1 (calibration)
  is excluded, B2 (amplitude) is measured at 0.195, B3 (band truncation) accounts
  for 69 % of σ². The direct fit **cannot** constrain `τ₂` on a ±9.6 ms grid.
* **Point 3 (H4’s 1.25 depends on an unjustified weight) — confirmed by
  measurement.** 1.2503 (full proper volume) vs **1.2860** (created volume) vs
  1.2417 (flat). My own docstring named the created volume and my code computed the
  full one — the same “word versus execution” gap as the harness defects.
* **Point 4 (report vs preregistration) — fixed.** The frozen T8 (false-alarm rate
  on 100 noise realisations) is now implemented as written and passes with 0/100;
  the Fisher condition number is renamed T8b and labelled post-hoc. Additionally
  found: `analysis.py`’s delay column is **bit-identical** to its `t_c` column, so
  its 4×4 Fisher matrix is exactly rank-deficient and its reported 1.8·10²⁰ is a
  float64 roundoff floor — my recomputation gives 3.1·10²⁰ on the same bytes, i.e.
  the value does not reproduce, which is why it is reported as a restatement of the
  exact degeneracy rather than a measurement.

## A property of the metric that had to be stated, not just noticed

The **absorbed fraction is not invariant** under the choice of delay reference
`t_ref` (0.9999967 at `t_ref = t(f_min)` versus 0.9530 at `t_ref = t(f_max)`),
while the **residual SNR is invariant** to 2.3·10⁻¹⁰. So the absorbed fraction must
never be quoted as the evidence for H3; only the residual may. This is now a
threshold in `pn_orders.py` (`H_absorbed_fraction_NOT_reference_invariant`),
because a metric that moves when the physics does not is a trap.

## Independent verification of this round

`verify_new_results.py` recomputes every new headline number by a different route:
the PN identity **symbolically** with sympy; the projection via **normal
equations** instead of QR; the H4 weights by **symbolic integration** instead of
quad; the ET/CE scaling exactness and the direct-fit flatness from the recorded
numbers. Then `--selftest` corrupts each of the four artifacts in turn and
requires the corresponding verifier to go red: **all four did** (`SELFTEST_PASS`).
The new Lean file `HomogeneousUniverse.lean` compiles with `NO_SORRY_NO_AXIOM`,
audits to only the three standard axioms, and carries three negative controls that
all fail as required.

---

# Round 221 — the owner’s objection about 37.3, and the τ₂ grid

The owner read §12.2 and §14.2 and made a substantive objection rather than a
style note: **the template’s optimal SNR of 37.3 (leading-order inspiral,
M_c = 28 M☉, D_L = 410 Mpc) is HIGHER than the published full-IMR single-detector
SNR (~20).** Either the normalisation is wrong, or 37.3 is an integral over a band
that includes the region where the inspiral formula is already invalid — in which
case “optimal SNR 37.3” is not physically meaningful and B2 is B3 in other
clothes.

He was right about the number and wrong about the identity, and both halves came
out of the measurement.

## Four errors, and every one of them is in what was being COMPARED, not in code

| # | Error | How it showed up | Caught by |
|---|---|---|---|
| 1 | `verify_rho_opt_bands.py` interpolated the PSD **linearly**; the audited pipeline interpolates it **log-log** (`psd_on_grid`) | a residual of **1.7e-3 on every band**, including the smooth ones — a systematic, not a band-edge effect | the convergence test refused to converge; measuring the two interpolants directly gave 2.8e-3 (linear) against 8.4e-4 (log-log) |
| 2 | The “convergence” test refined `df` on a band whose top edge **moved**, because f_ISCO = 67.6304 Hz is not a multiple of 0.25 Hz | a spurious non-monotone 6.3e-4 → 2.0e-4 → 1.4e-4 | the test’s own numbers: a convergence test on a moving interval is not a convergence test |
| 3 | Two of my first-draft criteria were simply **wrong**: “ISCO band < 20” (it is 20.71) and “IMR band within 8 of 20” (it is 36.58) | two false red cells | the measured values themselves. The wrong criteria are **recorded in the code**, not deleted, because the second one is the finding: cutting at 250 Hz does *not* fix the number, which is what makes the reduction a validity effect rather than a bandwidth effect |
| 4 | `h3_grid_bounds.py` NC1 compared the coarse-grid peak against the recorded `peak_snr` 7.2739 instead of the recorded `profile_at_zero` 7.2732 | a false red on a control that should pass | the comparison itself: 7.2739 comes from a two-stage FINE search over `t_c` and is a different quantity |

**The pattern across all four:** the arithmetic was right and the *identity of the
thing being compared* was wrong. That is the same failure mode the owner has been
finding in my work for several rounds — and this time two of the four were caught
by controls I had written for exactly that purpose, which is the first time the
controls fired on me rather than the reverse.

## What the measurement settled

* **37.3 is not a normalisation error.** `h3_rho_opt_bands.py` reproduces
  37.31184231798413 exactly (rel < 1e-6) in the pipeline’s own convention.
* **It is a band integral, and the owner’s reading is correct.** 69–71 % of the
  template’s σ² lies above the ISCO frequency. Cut at ISCO: **20.71 / 18.17** — the
  published value. Cut at 250 Hz: 36.58 / 32.87 — so the reduction is a validity
  effect, not a bandwidth effect.
* **But B2 is not B3 in other clothes.** The decisive test: if it were, cutting the
  band would drive the data’s preferred amplitude scale to 1. Measured 0.195 →
  0.339 (H1) and 0.166 → 0.291 (L1): the deficit shrinks from 5.1× to 3.0× and
  **stops**. B3 accounts for ~1.7 of the 5.1 and leaves a real factor ~3
  unexplained. That residue is the honest open item.
* **The grid item closed both ways.** Expanded ±9.6 → ±38.4 ms: still
  unconstrained. And a physics bound that no scan would have found: a quadratic
  delay must sweep less than π of phase across the band, capping the amplitude at
  **1.667 ms** — the paper’s 1.2 ms uses **72 %** of that budget.
* **And one thing that had to be said rather than left in the artifact:** on H1 the
  profile maximum sits AT the grid edge (−38.4 ms), i.e. it is still rising when
  the scan stops. That is not a best fit and not a bound — there is no interior
  maximum at all, which is the strongest form of “unconstrained”.

## Independent verification of this round

`verify_rho_opt_bands.py` builds its **own** PSD (own Welch code, own log-log
interpolant), computes f_ISCO by a **different algebra** (total mass in SI rather
than the chirp-mass formula), and uses a **different quadrature** (continuous
integral rather than discrete sum). Its primary check — an independently built
discrete sum against the recorded number — agrees at **rel = 0.00e+00 on all six
cells**. Four single-field corruptions of the artifact all turn it red
(`SELFTEST_PASS`). Both new artifacts reproduce **byte-for-byte** on a re-run.

---

# Round 224 — the factor-~3 residue closed with a real IMR template

The ceiling declared in §15.2 (“without an IMR waveform the residue cannot be
separated”) is lifted: phenomxpy 2.0.3 is installed. `h3_imr_check.py` replaces the
leading-order inspiral with a full `IMRPhenomT` on the pipeline’s OWN grid
(df = 0.25 Hz) and in the SAME discrete convention.

## What the measurement settled

* **The IMR template is 15 % lower** than the inspiral: 37.3118 → **31.6595** (H1),
  33.7097 → **28.7350** (L1) in 20–300 Hz.
* **But that is not enough to reach ~20.** The threshold “within 25 % of 20”
  **FAILED** (31.66 vs 20 is 58 %). Recorded as a failed threshold, not hidden.
* **The data’s preferred amplitude scale moved to 0.5288 (H1) / 0.4470 (L1)**,
  from 0.1949 / 0.1655. So **more than half** of the old 5.1× deficit was
  leading-order template imperfection (B4) — this is the quantitative naming of B4
  that §15.5 item 2 was missing.
* **The residue is orientation, and it is measured.** rho_opt falls monotonically
  with inclination: 31.66 (0°) → 19.79 (60°) → 15.83 (90°), crossing rho_opt = 20
  at **ι ≈ 59.19°**, inside the physically allowed range.

## Three errors of this round

| # | Error | How it showed up | Caught by |
|---|---|---|---|
| 1 | multiplied the inspiral amplitude by `fs` and did not divide by `fs` in the matched filter | NC1 gave rho_opt = 297708 instead of 37.31 — off by 8000× | **NC1** |
| 2 | NC5 asserted “the model is strictly zero below its own f_min” | **false**: phenomxpy carries a tail to 2.6e-23; measured 0.52 % of band power | the control itself (it did not go red where it should have) |
| 3 | the verifier built the model at df = 0.5 and interpolated linearly | **19 %** of the band power destroyed | the 2.2 % disagreement with the audit |

**The third is the useful one.** The FD phase rotates by ~π between adjacent
0.5 Hz samples, so linear interpolation of a complex oscillating function is not
interpolation. Measured: at SHARED frequencies the two grids agree to **1e-16**,
while after interpolation the power ratio is **0.810**. So the **verifier was wrong
and the audit was right**, established by measurement rather than argument. The
false NC5 claim and the wrong 1e-3 threshold are **recorded in the code, not
deleted**.

## Controls vs findings, kept apart on purpose

A failed CONTROL means a broken instrument and a void run; a failed FINDING is a
result. The exit code gates on the controls only. All five controls pass: NC1
reproduces the recorded 37.31184231798413 / 33.709713860505396 to <1e-6, NC2 IMR
rho_opt finite and positive in every band, NC3 amplitude linearity (2.000000000),
NC4 zero waveform → 0, NC5 below-f_min power 0.52 % (< 1 %).

## Independent verification

`verify_h3_imr.py` imports nothing from the audit: own Welch, own log-log PSD
interpolant, own `IMRPhenomT` call, own quadrature. Independent IMR rho_opt matches
the recorded value at **rel = 0.00e+00** on both detectors; the crossing
inclination reproduces to 1e-6 and independently gives rho_opt = 19.9817 there.
Five single-field corruptions of the artifact all turn it red (`SELFTEST_PASS`),
and the artifact is restored byte-identically. The artifact reproduces
**byte-for-byte** on a re-run: sha256 `36ed9094…`.

## What remains open

1. **Spins = 0** and the base number at ι = 0 are stated simplifications; the IMR
   rho_opt is an upper bound, not the event’s recovered SNR.
2. **The delay is still not fitted** to the strain; only the template’s optimal SNR
   is measured here.
3. **ET/CE** — unchanged: the sign is undetermined without tau(t).
4. **Deriving tau(t) from the theory** — not done; the paper does not state it.

---

# Round 226 — both items of §16.5 closed

## Five errors, all caught by a control

1. **`np.fft.irfft` instead of `np.fft.ifft`** in the matched filter. `irfft`
   assumes a conjugate-symmetric spectrum and returns `2·Re(Σ)`, a different
   quantity. The peak was underestimated by exactly 16384 (= fs), and a 20 ms
   injection was “recovered” at the grid edge. Caught by **NC6** (comparison with
   an explicit matrix).
2. **`B = np.zeros(seg_n)` indexed by an rfft-length mask** (`seg_n//2+1`). NumPy
   does not object to a mask shorter than the array — it silently takes the prefix,
   so the error was quiet: all profiles ~10⁵ too small, the injection always at the
   edge. Caught by **NC6**.
3. **`template_time_series` divided by `fs`.** The noiseless round trip gave σ/fs =
   0.0019 instead of σ = 31.66. Caught by **NC7** (a round trip with target SNR 25
   returned 0.0015).
4. **A ±9.6 ms grid for a 20 ms injection.** “Recovered 9.6 ms” is the grid edge,
   not a result. Caught by **NC4** (a `recovered_at_scan_edge` check was added).
5. **Verifier: two-stage peak search.** A coarse pass with stride·8, then a
   refinement in a ±16 window. On H1 it **silently missed the true maximum**:
   returned 7.8106 where the correct value is 8.2011 (5 %). A speed optimisation
   that changes the answer is not an optimisation. Deleted; the fast path is now
   exact (inverse DFT), and the matrix path is kept as a spot check.

## False criteria, recorded rather than deleted

**NC6 of the first draft** (`h3_imr_spins.py`): “a template with more band power
cannot have a lower matched-filter peak on the same data”. **Measured: false.** The
published spins raise ρ_opt 31.66 → 34.84, while the filter peak **falls** 16.74 →
13.51: the spin also changes the template’s **phase**, so the overlap with the data
is not a function of power alone. The criterion was replaced by the Cauchy–Schwarz
inequality (peak ≤ ρ_opt), which is true and checkable.

**NC1b** — not a false criterion but the strongest control: LO+pn+expanded fit must
reproduce the recorded numbers of `h3_grid_bounds.json`. Without it the new IMR
numbers would not be comparable with the old ones.

## Verifier errors found by measurement

* **V2** first compared against **unconstrained** angles (|a1|, |a2| at their own
  bounds), which the artifact deliberately excludes via χ_eff. Discrepancy 12–19 %.
  Fixed: the verifier builds its own envelope with the χ_eff constraint.
* **V5/V6** first used a sqrt-hann window instead of tukey 0.25 and disagreed by
  1.6 % on the profile drop. Measured: `Σw²` = 32767 (sqrt-hann) versus 55295
  (tukey 0.25). The window is part of the **definition** of the pipeline;
  independence comes from the quadrature, not the window.

## Numbers that did not exist in the project before this round

* Published single-detector SNRs: **19.5 (H1)** and **13.3 (L1)**, not “~20”.
* 24 and 25.1 are **network** quantities; comparing a single-detector ρ_opt to them
  is a category error.
* The published SNR ratio 1.466 versus the ratio of optimal template SNRs 1.102 and
  the filter-extracted 1.303; the residue is the antenna response.
* The inclination multiplier is **(1+cos²ι)/2 exactly** and independent of spins.

---

# Round 227 — joint H1+L1 fit with antenna responses

Owner input (msg_00227): (1) make a **joint** H1+L1 fit with **one** `τ₂` and free
`(M_c, t_c, φ_c)`, and **only after that** speak of a bound; (2) use the antenna
responses `F₊, F×` for the GW150914 position **instead of** scanning over the
inclination `ι`; (3) rewrite §1 and §7 to the current state, update the file table,
move the chronology to `NOTES.md`; (4) draw up a separate list of claims tagged
“proven / measured under such assumptions / hypothesis”.

Item (4) was done as a separate file, `work/CLAIMS.md` — and it is what determined
the shape of §1: the headline may rest only on **[PROVEN]**, while every numerical
bound is tagged **[MEASURED]** with its assumptions listed.

## Two structural facts without which the fit cannot be assembled

**The first — the phenomxpy convention.** Measured: `IMRPhenomT` returns the
polarizations as

```
hp(ι, φ) = A(ι)·e^{+2iφ}·hp₀ ,   hc(ι, φ) = B(ι)·e^{+2iφ}·hc₀ ,   hc₀ = −i·hp₀
```

(`A = (1+cos²ι)/2`, `B = cos ι`; the identity `hc₀ = −i·hp₀` was verified, maximum
deviation 3.3·10⁻⁴). Hence the detector strain

```
F₊·hp + F×·hc = e^{2iφ}·hp₀·( F₊·A − i·F×·B )
```

i.e. **the whole response transformation collapses to one complex coefficient**
`C_d = F₊,d·A(ι) − i·F×,d·B(ι)`.

**The second — from the first.** Since `h_d = C_d·h` with a common `h`, the
coherent network statistic has the form

```
ρ_coh² = | Σ_d conj(C_d)·ζ_d |² / Σ_d |C_d|²·s_d²
```

where `ζ_d` is the complex matched-filter output for the waveform `(1,0)` and `s_d²`
its norm. `ζ_d` does **not** depend on `C_d`, so it is computed once per `(M_c, τ₂)`
per detector, after which the whole `(ι, ψ)` grid is free. The fit becomes genuinely
coherent: the amplitude and phase ratio of the two detectors is **predicted** by the
sky position, not a free handle.

**A consequence named as a finding rather than hidden.** `φ_c` multiplies `C_d` of
each detector by the **same** phase, so it is identically degenerate with an overall
phase and **not identifiable** in a coherent fit with a fixed sky position. It is
profiled, not quoted.

## Sky position: the reading chosen by measurement, not by trusting a library

The position is taken from the **published** map `LALInference_skymap.fits.gz`
(LOSC, P1500227). Two independent HEALPix readers — `healpy` and `astropy_healpix`
— agree on every pixel position to **1.8·10⁻¹⁵ rad**.

But their arrays **differ**, and that turned out to be a trap. The FITS table stores
`PROB` as a 2-D array; `data['PROB'].reshape(-1)` (row-major) and the array returned
by `healpy.read_map` place the probability on **different** parts of the sky. No
library was taken on trust — the reading was chosen **by physics**: the published
constraint (arXiv:1602.03840) requires the source to lie on the annulus of constant
H1−L1 arrival-time difference = 6.9 (+0.5/−0.4) ms.

| reading | concentration on the annulus | at the delay |
|---|---:|---:|
| row-major + NESTED | **0.1635** | **−6.90 ms** ← accepted |
| `healpy.read_map` + NESTED | 0.0233 | +3.70 ms |
| shuffled map (20 runs) | 0.0061 avg, 0.0064 max | — |

The accepted reading is **25×** the null of a shuffled map. It was also measured
that the projection requires a GMST rotation (the map is in celestial coordinates,
the detector vertices in the Earth frame): without the rotation the peak sits at
−3.4 ms with concentration 0.026.

**Checks on the position.** ML position RA = 134.80°, Dec = −69.79°; 90 % area
**616.4 deg²** (published 610), 50 % — 149.1 deg² (published 150);
`τ_H1 − τ_L1 = −6.898 ms` (published 6.9). And the decisive one: the optimal SNRs
at `ι = 0` are 25.31 (H1) and 17.16 (L1), their ratio **1.4749** against the
published single-detector **19.5/13.3 = 1.4662** — agreement 0.6 %. That is the
check that the sky position is the right one. (But see `CLAIMS.md` B7: this is an
agreement, **not** a 0.6 %-precision test — the published numbers carry ~7 % noise,
so the real checks are the delay, the area, and the shuffle control.)

## What the joint fit gave

| delay kernel | grid | profile drop | one-sided bound |
|---|---|---:|---|
| pn | ±9.6 ms | 1.375σ | none |
| pn | **±38.4 ms** | **3.496σ** | **yes**: +3.6 ms |
| imr | ±9.6 ms | 0.852σ | none |
| imr | **±38.4 ms** | **2.352σ** | **yes**: +2.4 ms |

Peak network SNR: **22.93** (pn kernel) and **23.07** (imr kernel); at `τ₂ = 0` —
**22.17**. For comparison, the published network SNR of GW150914 is **24**
(arXiv:1602.03837) — the model reproduces the network quantity to ~4 %, whereas the
earlier single-detector numbers were compared with single-detector 19.5/13.3.

**Power control — without which the bound is worth nothing.** The same fit
**recovers** a 20 ms injected delay (20.3 ± 5.1 ms) and **does not resolve** the
paper’s 1.2 ms: spread **4.6 ms** at an amplitude of 1.2 ms.

## My errors of this round — six, all caught by a control

| # | error | how it showed up | caught by |
|---|---|---|---|
| 1 | light-travel time between the detectors not accounted for | coherent SNR 12.0 versus incoherent 21.1 — coherent **lower**, which cannot be if the phase relation is right | NC9 + a direct measurement of the sign |
| 2 | the sign of the light-travel term was guessed, not measured | the first variant gave 14.3; the correct `exp(+2πifτ_d)` gives **21.08** against an incoherent 21.10 | measurement over the `(ι, ψ)` grid |
| 3 | NC1 compared the IMR template against the **recorded number of the LO template** | a control failing for an irrelevant reason | cross-check with `h3_direct_fit_imr.json` per template |
| 4 | NC8 “not recovered if the mean differs from 1.2 ms by more than 1 ms” — a **false criterion** | the mean can land on 1.2 ms by chance with a spread of 4.6 ms, which is the **opposite** of resolution | replaced by an honest one: spread greater than the signal |
| 5 | hand-written HEALPix (two attempts) | pixel 0 came out at `θ = 0.09°` instead of `89.93°` | instead of debugging, two libraries + a physical test |
| 6 | verifier conflated the conventions `rfft(x)` and `rfft(x)/fs` | every check went red with a factor `fs` | V6 (the recorded numbers) |

**The fifth is the most useful.** Rewriting HEALPix by hand is exactly the class of
error this project has caught round after round; so the hand-written implementation
was **deleted**, not debugged, and replaced by two independent libraries plus a
physical test of the reading.

**The fourth is recorded, not deleted:** the false criterion stays in the code with
a note, as do the previous four.

## Independent verification (`verify_h3_joint.py`)

Imports nothing from the audit. Own PSD (own Welch), own GMST (a different algebraic
form), own beam algebra (the tensor `D` written component by component), an
**explicit matrix** `exp(2πift)` instead of an inverse DFT, own quadrature.

| check | result |
|---|---|
| V1 antenna responses | 5.4·10⁻⁸ |
| V2 ratio of optimal SNRs | **rel = 3.9·10⁻¹¹** |
| V3 network SNR at `τ₂ = 0` | **rel = 7.3·10⁻¹¹** |
| V4 profile in `τ₂`, both kernels | peak 1.4·10⁻⁵, drop 6.5·10⁻⁴ |
| V5 round trip | 7.5·10⁻⁶ |
| V6 recorded numbers | **rel = 0.0** |
| **self-test: corrupt the artifact** | **goes red** |

**Measured separately, not hidden:** the PSD estimator is part of the *definition*
of the quantity, not an implementation detail. With a **mean** Welch instead of a
median one the single-detector SNRs move by tens of percent, because the O1 release
contains non-Gaussian glitches (already measured in round 219). Independence comes
from the **quadrature** and the **code**, not from changing the estimator; so the
verifier uses the same estimator, and the sensitivity to it is named as a number.

## What remains open

1. **`τ(t)` is still not derived from the theory** — the paper does not state it.
   The main open item of the whole cycle; until it exists, what is proven is
   degeneracy *within the adopted linear model*, not the impossibility of testing
   the whole “now”-theory.
2. **Aligned spins**; precession `χ_p < 0.71` is not covered.
3. **Phase-only delay** — the amplitude is not delayed.
4. **`η` fixed** at the published 29/36; a joint fit over `(M_c, η, spins, ι)` would
   widen the profile. (Round 228 correction: for the claim “a bound exists” this is
   the **unsafe** direction, which is why the bound is labelled an upper estimate.)
5. **The factor-~3 residue** in the 37.3↔7.27 gap — closed by the full IMR
   (`h3_imr_check.py`, round 224); the obsolete “not separable” phrasing is
   withdrawn.
6. **L1 does not reach 13.3** even at `ι = 90°` (minimum 13.58).

Artifacts: `artifacts/h3_joint_fit.json` (`art_32983e9ae39a`, sha256 `1daa8407…`),
`artifacts/sky_samples.npz` (`art_d8f822f55257`, sha256 `0afea4a2…`). Both
reproduce **byte-for-byte** on a re-run. Claim ledger — `work/CLAIMS.md`.

---

# Round 228 — ledger revisions per the owner’s review

The owner (msg_00228) found six discrepancies between `CLAIMS.md`/`REPORT.md` and
the facts. None was an error of physics — all six were errors of **formulation**.
The first was found from the **ledger’s legend**, not from the numbers.

## What was fixed

1. **A4–A6 stood under the tag [PROVEN] without a Lean file.** A fourth category,
   **[DERIVED symbolically]** (sympy + a negative control, no formal proof), was
   introduced. §F was rewritten: “in Lean, A1–A3 are proven”.
2. **B1 is an upper estimate, not a measured bound.** Free `η` and spins would
   absorb more of the quadratic term (with `η` fixed and zero spins, 88 % is
   already absorbed). Added to the assumptions: the sky position is fixed at the
   maximum-likelihood point (area 616 deg²).
3. **E4 was rewritten.** “The safe side” is true for “no bound” but not for “a
   bound exists”. For B1, fixing `η`/spins works in the unsafe direction.
4. **“+3.6 ms” was defined.** It is the right crossing of `peak − 1σ` at an interior
   peak. On ±9.6 ms the peak is at the left edge → the right crossing is +5.4 ms;
   on ±38.4 ms the peak is interior (−18.0 ms) → +3.6 ms.
5. **B7: 0.6 % is agreement, not precision.** The ratio 19.5/13.3 carries a ~7 %
   tolerance. The real checks on the position are the delay, the area and the
   shuffle control.
6. **§1: the fit resolution is ≈ 5 ms** (spread 4.6 ms at 1.2 ms, 5.1 ms at 20 ms);
   69 ms is tagged as a **leading-order** result.
7. **The §9 ↔ E7 contradiction was removed.** One phrasing kept: the factor-~3
   residue is **closed** by the full IMR.

## Two errors of this round — both in the verifier, not the report

**Error 1 (V3).** The first criterion searched `lean/*.lean` for the **word**
“delay”/“tau” and went red — because `Identifiability.lean` is *devoted* to an
additive delay in the abstract, so the word is there, while the **content** of
A4–A6 (the coefficients `τ₀/τ₁/τ₂`, the chirp phase, Euler’s theorem) is not. The
criterion was rewritten to test content; a positive control V3e was added: the word
“delay” **must** be present — otherwise the test would not distinguish content from
wording.

**Error 2 (V6).** The regex `0\.(?:8849|8716)` did not find H1, because the
recorded value is `0.88487345…`, i.e. the digits `88487`, not `8849`. Numbers are
now read from the JSON **by key**, not by pattern.

Both are the same class as before: “I am comparing something other than what I
think I am”. They were caught because the criterion went red on **correct** data,
which forced a check of the criterion rather than the data.

## Verification

`verify_claims_228.py` (imports nothing from the audit): V1–V6 recompute all six
revisions from the frozen artifacts; three negative controls. `--selftest` corrupts
a **copy** of `h3_joint_fit.json` (driving the peak to the edge) and requires the
verifier to go red — **SELFTEST_PASS** (V1a and V1c go red).

---

# Round 229 — the publication package, and a finding about conventions

**Owner input (msg_00229):** produce a PDF paper and translate the Russian
artifacts into English, keeping the Russian alongside with English primary; it will
then be published on GitHub and Zenodo; the authors are Aiodam (autonomous research
agent) and Oxunjon Ubaydullayev; and it must be stated plainly that the work was
done by the agent — that the agent walked the whole path itself.

## What was produced

`publication/` — a self-contained package: the paper in English (primary) and
Russian, both as Markdown and PDF; English translations of the report, the claim
ledger, the pre-registration, this error log and the round summaries; a README, an
artifact index with hashes, a manifest, `CITATION.cff`, `.zenodo.json`, licences;
and an independent verifier with its own self-test. The Russian originals in
`work/` are untouched.

Authorship is stated plainly in the paper (§0), the README and `CITATION.cff`: the
work was carried out by the autonomous research agent Aiodam as primary
investigator, with Oxunjon Ubaydullayev as co-author who set the question and
supplied the data.

## The finding: the π-bound depends on the convention by a factor of 15

This is the substantive discovery of the round, and it came out of a **check**, not
out of writing prose.

While re-deriving the π-phase bound in a fresh process (without reading the
artifact), I got **24.9635 ms** where the artifact gives **1.6667 ms** — a factor of
**15**. I did not smooth this over; I worked out which of the two is right.

**The cause is the time-reference convention, not arithmetic.** The paper’s delay
`τ(t)` is a function of the **emission time** `t`, with `t = 0` at coalescence. The
paper’s picture is that the effect **accumulates** toward merger — small early in
the inspiral, largest at merger. In those variables:

* **the paper’s convention:** `τ = τ₂ (t − t_min)²`, vertex at the **start** of the
  band (delay 0 at 20 Hz, growing toward merger) → the π-criterion gives
  `τ₂ = 0.0023370 s⁻²`, amplitude **1.667 ms**, and the paper’s 1.2 ms uses **72 %**
  of the budget;
* **the opposite convention:** `τ = τ₂ t²`, vertex **at coalescence** (delay 0 at
  merger, largest at 20 Hz — backwards relative to the paper) → the same criterion
  gives **24.96 ms**.

Checked independently: in the artifact’s convention the added phase is monotone
across the band, so its maximum and its range coincide (1.667 ms); in the opposite
one the vertex at zero gives 24.96 ms; ratio 14.98.

**Conclusion.** The artifact is **right for the paper’s convention**; my
re-derivation was right for its own. So “1.667 ms / 72 %” is a statement **about the
paper’s convention**, not a convention-free fact, and the paper now says so (§5.5)
with 24.96 ms given alongside. It became limitation 8 in §9.

**Why this matters rather than being pedantry.** The same question — “which
convention?” — was already answered substantively in this project: in round 226 the
group-delay ratio `t_IMR/t_PN` moved from 1.03 to 155 depending on whether the
inspiral or the IMR form was used (`h3_direct_fit_imr.py`), and both conventions
were computed rather than one chosen. The same rule now had to be applied to the
phase bound — and it was substantive again: **without a named convention the
π-bound is undetermined by a factor of 15.**

## What the verifier checks, and a hole it found in itself

`verify_publication.py` imports nothing from the audit. It reads the frozen
artifacts **by key** and the published documents, and requires every headline number
in the paper to be the artifact’s number. It also checks the author names, the
scope caveat, the four ledger categories, that the metadata parses, that every hash
in `ARTIFACTS.md` matches the bytes on disk, and **Lean traceability** — every
theorem name the paper cites must be **defined** in `lean/*.lean`.

Result: **156 checks, 0 failures** — `PUBLICATION_CONSISTENT`.

`selftest_publication.py` corrupts **copies** in an isolated tree and requires the
verifier to go red, while requiring it to stay green on the untouched tree.
**10 of 10 cases behaved as required** (`SELFTEST_PASS`), including C8 (a theorem
renamed in Lean, so the paper would cite a proof that no longer exists) and C9 (the
convention caveat deleted).

**The self-test found a hole in the verifier itself:** C1 initially **survived** —
the corrupted token `7.27` occurs twice in the paper, so replacing one occurrence
left the check satisfied. That is exactly the class this project keeps catching:
the check was testing **presence**, not what the author thought it was testing.
Fixed: C1 now uses a token occurring exactly **once** (`0.174`), and the limitation
(presence, not placement) is written into the self-test’s docstring rather than
hidden.

## What remains, and what was not done

1. **The DOI and repository URL are not filled in** — the owner does that at deposit
   time; explicit placeholders are marked in `.zenodo.json` and `CITATION.cff`.
2. **The licence** is set to CC BY 4.0 (text) + MIT (code) as a customary default
   for Zenodo deposits; this is the owner’s choice, named in the README and in the
   `.zenodo.json` notes rather than assumed silently.
3. **Publication itself was not performed** — neither GitHub nor Zenodo. I have no
   authority for that, and the owner said he would publish it.
4. **`τ(t)` is still not derived from the theory** — the main open item, unchanged.
