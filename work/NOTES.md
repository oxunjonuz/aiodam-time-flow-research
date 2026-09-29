# NOTES — what went wrong on the way (v1 → v5)

This file records my own errors, in order, because the report is only worth
reading if the path to it is visible. Every entry was found by a check, not by
re-reading my own code.

## v1 of `analysis.py` — four errors, all caught by looking at the numbers

| # | Error | How it showed up | Fix |
|---|---|---|---|
| 1 | `Rs62_m = 2*G*MF/C**2` with `MF` already in **seconds** (geometric units) | `Rs62_km = 4.5e-34` — a black hole smaller than a proton | `Rs = 2M` in geometric units → 183 km |
| 2 | Residual computed as a **phase-space** norm, not a strain norm | `unmodelled_snr = 9.7e32` — an SNR larger than the number of atoms in the Sun | multiply by the strain amplitude `h(f)` |
| 3 | `np.linalg.lstsq` on raw columns spanning ~40 orders of magnitude | `max_abs_residual = 2e8` for a direction that is *exactly* in the span | orthonormalise via SVD first |
| 4 | "False-alarm rate" test compared a **phase-space** SNR to a threshold of 1 | 100 % false alarms on pure noise | deleted; replaced by a deterministic flat-likelihood test |

## v2 — the SNR was still wrong

The projection was correct, but the statistic was the norm of a *phase*
residual. A phase residual is not an SNR. Rebuilt around `u = h · ΔΨ` so that
`‖(I−P)u‖_weighted` is a genuine matched-filter SNR, directly comparable to the
detection threshold ≈ 5.

## v3 — the PSD was not calibrated

The analytic PSD gave an in-band signal SNR of **129** for GW150914, whereas the
published single-detector value is **20**. Any SNR computed against it was
inflated by 6.5×. Fixed by calibrating the PSD with one constant so the injected
signal's SNR equals 20. All headline SNRs are now calibrated.

## v4 — three broken controls

| Control | What was wrong |
|---|---|
| K9 | the "free" parameters absorbed the injected signal, so the residual was flat for the *wrong* reason — it looked like a pass and was a failure |
| K3 | threshold `> 5` was arbitrary; at `τ̇ = 0.01` the non-absorbable part is only SNR 1.8 |
| K8 | the noise was scaled so that its own SNR was ~63, not 20 |

Fixed: K9 now fits only `(t_c, φ_c)` with `M_c` fixed; K3 is reported at two
values of `τ̇`; the noise is normalised to the published SNR.

## v5 — the mutation control found 8 holes in my own test

The automated campaign (`mut_dcd2ff7529de`) returned a **null result**: its
coverage tracer reported "0 executed lines" and it classified all 20 mutants as
unexecuted without running one. That is an instrument failure, not evidence, so
I wrote `mutation_control.py` by hand.

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
baseline, with a tolerance-aware comparison and a two-tier split between
headline claims and scale-dependent diagnostics):

* **M3, M6, M16, M19 — equivalent on every headline claim.** `absorbed_fraction`,
  `unmodelled_snr`, `argmin_candidate`, `K8_residual_spread` and the geometry
  numbers are identical to the baseline to within 1e-12 absolute / 1e-9 relative.
  M19 (sign flip) is identical in *every* leaf.

**One honest limitation this exposed.** The Fisher condition number and its
singular values are **not** invariant under these faults: M16 moves the condition
number from 1.8e20 to 4.5e26 while every headline claim stands. That is expected
— the condition number is a property of the parameter basis, and rescaling a
basis column changes it. The consequence for the report is concrete: **the
condition number is reported as an order of magnitude (> 1e12, i.e. numerically
singular), not as a measurement.** The same caveat applies to
`K9_over_K8_sharpness_ratio`, which is a ratio of two residuals and therefore
inherits the scaling of the K9 basis direction.

## v2 of `h4_paper_arithmetic.py` — a sampling artifact, found by the owner

The owner read the report and caught a concrete error in §5/H4: I wrote that the
factor 1.25 is unreachable "at any radius in the stated range 1.5–4 Rs". **That
was false.** The code evaluated the redshift factor `1+z = 1/sqrt(1 - Rs/r)` at
only four radii `{1.5, 2, 3, 4}` and tested `|1+z - 1.25| < 0.02`, which those
four samples miss. The analytic root of `1+z(r) = 1.25` is `r/Rs = 2.7778`, and
it lies **inside** `[1.5, 4]`.

The sharpest part: **v1 computed `r_needed = 2.7778` itself and printed the range
`[1.5, 4]` in the same artifact, but never compared them.** Two numbers sat next
to each other in one JSON file and the verdict contradicted both. This is the
same failure mode as the harness turns: a claim that no check ever executed.

Fixed in v2: the scan covers the whole interval (20001 points), locates the root,
and carries a negative control on `[1.5, 2.5]` (which excludes the root) that
must report `contains_1.25 = false`. `verify_h4_independent.py` recomputes the
root **symbolically** with sympy (the audited script does it by hand) and
disagrees if the artifact's root is wrong — verified by corrupting the artifact,
which turns the verifier red. The original v1 is archived at
`work/archive/h4_paper_arithmetic_v1_original.py` (sha256 `15f151c1…`) and
reproduces the false verdict on the same input.

**What did NOT change:** the main H4 result — 0.6 ms and 1.2 ms are not
independent estimates — stands, because it never used the redshift factor.

## What the independent verifier caught

`verify_independent.py` found a **sign error in my own assertion** about the
`τ̇` direction: I wrote that the `τ̇` term is `−M_c · dΨ/dM_c`, but sympy gives
`+M_c · dΨ/dM_c`. The algebra in the analysis was right; the check text was
wrong. That is precisely the failure mode an independent path exists to catch.

---

# Round 219 — H3 on the real strain, and H4 against the paper's own text

The owner said h5py was installed, both GW150914 files were readable, and asked
for H3 to be checked on H1/L1 and H4 against the source materials. Six errors
were made and caught on the way; five of them were caught by controls I had
already written, and one by an independent verifier.

## H3 on real data — five errors, all caught by controls

| # | Error | How it showed up | Caught by |
|---|---|---|---|
| 1 | Matched-filter normalization `C = 4*df` (the continuous-FT form) with an inverse-FFT reconstruction | **noise-only peak SNR ≈ 4.6·10⁴** — an SNR no detector can produce | the noise-only control (must be ~1) |
| 2 | Template amplitude convention: used `h~(f)` directly as the rfft, missing the factor `fs` | template optimal SNR **1.7·10⁻⁷** instead of ~45 | the analytic σ comparison |
| 3 | No window on the segment | peak ran to the **segment edge** (ρ = 306) because the rectangular edges leak into the band | the check that the peak lands near the merger |
| 4 | Synthetic validation noise generated from the full PSD but **not** high-passed, while the real data was | noise-only peak 8–13 instead of ~4.8 | the calibration's own noise-only distribution |
| 5 | Verifier computed the absorbed fraction as `1 − e_res/e_tot` (absorbed **energy**) instead of `sqrt(1 − e_res/e_tot)` (= ‖Pu‖/‖u‖) | verifier said 0.783 where the artifact said 0.885 | the verifier's own cross-check |

**What settled each one was a measurement, not algebra.** `mf_norm_probe.py`
settled the normalization by generating noise from a known PSD and reading off
which candidate gives ρ ≈ 1 and recovers an SNR-20 injection as 20; the answer
was `H_k = fs·h~(f_k)` with `C = 4·df/fs²`. The continuous-FT form `4·df` is
correct only for `X_cont = dt·X_rfft`, so it is off by `fs²`.

**The PSD estimator was also settled by measurement, and I changed my mind
twice.** `psd_whiten_check.py` compares the FFT-binned in-band variance with the
PSD integral. With a median estimator the real off-source data carries ~6× more
in-band power than the PSD predicts, while a synthetic series generated from that
PSD reproduces the integral to 0.98 — so the excess is in the data (non-Gaussian
transients in the released strain), not estimator bias. I first switched to the
mean estimator on that basis; then the mean gave ASD(100 Hz) = 2.3·10⁻²² /rtHz,
**22× the published aLIGO O1 sensitivity**, because the mean is not robust to
those same transients. Reverted to the median, which gives 1.03·10⁻²³ — the
published value. The final choice is stated in the code with both numbers.

## H4 against the paper's text — one of my own earlier claims is refuted

Reading the two source paragraphs directly (not from memory) showed that the
paper's chain is: bare 0.6 ms → "of order 1 millisecond", and bare 1.2 ms →
1.5 ms. Two things follow:

* **H4b is REVISED, and it weakens the criticism.** Only **one** precise factor
  is claimed (1.5/1.2 = 1.25). The apparent second factor, 1.0/0.6 = 1.667, is
  an artefact of the paper writing "of order 1 millisecond" — a rounded phrase.
  Treating it as a second derived factor was over-reading the text.
* **H4c is REFUTED, and it is my own error from turn 217.** I claimed the factor
  1.25 was obtainable "only by an unjustified choice of radius" r = 2.78 Rs. That
  is wrong. The **volume-weighted mean of (1+z) over the paper's own stated shell
  [1.5, 4] Rs is 1.25034** — the paper's 1.25 to 0.03%. The natural weighting
  (time created ∝ volume created) reproduces the paper's unexplained number from
  the paper's own stated geometry. The `r = 25/9 Rs` root is a red herring: it is
  where the *pointwise* redshift equals 1.25, which is not what a volume-weighted
  average gives.

  This was found by writing the new audit, not by re-reading my old one — and it
  is recorded as a correction of my own earlier claim, with the old claim quoted
  in the artifact next to its refutation.

**What survives H4:** H4a. The paper's two BARE numbers are not independent
estimates — their ratio is the geometric factor (ΔV_total)^{1/3}/Rs(62) = 2.0083,
and the paper's own bare ratio is 2.0. So "the two estimates agree" is arithmetic.
And the 1.25, while now explicable, is still introduced with no derivation in the
text: the reader is given a number and told only "with the additional
gravitational redshift".

## One more error, this time in the new independent verifier

`verify_independent_v3.py` initially reported the geometric factor as **3.506**
against the artifact's 2.008. The verifier was wrong: it gave each of the three
masses its **own** shell `[1.5 Rs(M), 4 Rs(M)]`, whereas the paper's eq. (4.3)
evaluates all three on the **same** physical shell `[1.5 Rs(62), 4 Rs(62)]` —
that is the whole point of subtracting the smaller holes' volume excess from the
region the 62 M☉ hole now occupies. Fixed, and the 2.0083 was then confirmed by a
second, symbolic route (`r = Rs62·u` substitution, `k = 29/62` and `36/62`).

**This is the useful kind of disagreement.** The verifier failing loudly on a
quantity I had asserted since turn 216 is exactly what it is for; had it agreed,
I would have learned nothing about which of the two was right.

## Negative controls on the verifier itself

Eight single-field corruptions of the artifacts, each applied in isolation: every
one turned `verify_independent_v3.py` red (H3 peak SNR, H3 residual SNR,
noise-only maximum, ASD(100 Hz), the geometric factor, the volume-weighted mean,
the 1.25 radius, and the NC3 shell-dependence flag). A verifier that cannot fail
proves nothing, so this was checked rather than asserted.

---

# Round 220 — higher PN orders, ET/CE, direct fit, and the owner's four review points

The owner accepted the GWOSC provenance point (hashes matched, so the caveat is
withdrawn) and asked for the next step: (a) higher PN orders plus an ET/CE
estimate, and optionally (v) the homogeneous-universe statement as a theorem. The
same message carried four review points on my earlier report.

## Errors of this round, all caught by controls or by re-reading the numbers

| # | Error | How it showed up | Caught by |
|---|---|---|---|
| 1 | The 3PN log term `−(6848/21)·ln(4v)` was subtracted **outside** the `v⁶` factor instead of multiplying it | the 3.5PN phase at 20 Hz was wrong by ~27×, and the reported band time span came out **32.5 s** instead of 0.84 s | the `t_span_plausible` self-check (0.5–2 s) |
| 2 | The identity check `f·dΨ/df = M_c·dΨ/dM_c` used `np.gradient` for the derivative | a spurious 0.46 % "failure" of an identity that is exact | a sympy route + an analytic derivative through `v` |
| 3 | The symbolic **negative control** used `ln(M_c·f)` | the control did **not** break the identity, because `M_c·f = v³/(πη^{-3/5})` IS a function of `v` alone — so the control was vacuous | the control itself failed to go red, which is why it was noticed |
| 4 | A leading-order test compared a 3.5PN delay direction against a leading-order basis | a spurious `4.8·10⁻⁹` leak where the identity requires machine zero | re-deriving per-order (each order tested in its own basis) |
| 5 | `h4_weight_audit.py` applied `1/sqrt(1 − 1/r)` with `r` in **metres** | all three weights came out `1.000001` — the `Rs/r` term underflowed to zero | the three weights were then *identical*, which contradicted the expected spread |
| 6 | `f_isco` mixed geometric and SI units (`c³/(6√6 π G M)`) | ISCO frequency `2.7·10³⁷ Hz` | the B3 fraction-of-σ² above ISCO came out exactly 0 |
| 7 | The direct-fit 1σ interval returned `None` when the profile never dropped by 1 | looked like a bug; it is the **result** — the data do not constrain `τ₂` at all | reported as `unconstrained_on_grid` instead of `None` |
| 8 | `verify_new_results.py`'s symbolic W2 mean returned a complex number | `TypeError: Cannot convert complex to float` | asserted the imaginary part < 1e-12 instead of silently dropping it |
| 9 | `prereg_t8.py` scaled the noise so the residual norm came out ~7·10²³ | 100/100 false alarms — an obviously broken statistic | the false-alarm rate itself (should be ~0) |

## What the owner's four review points changed

* **Point 1 (what the theory predicts) — accepted, not closed.** The honest
  ceiling is now written into §14.1: what is proved is degeneracy *in the adopted
  linear model*, not that the whole "now" theory is untestable. What this round
  *adds* is that the degeneracy is not a leading-order artefact (§12.1).
* **Point 2 (H3 does not fit the delay to the signal) — done.** `h3_direct_fit.py`
  fits `τ₂` directly and decomposes the 37.3 ↔ 7.27 gap. Result: B1 (calibration)
  is excluded, B2 (amplitude) is measured at 0.195, B3 (band truncation) accounts
  for 69 % of σ². The direct fit **cannot** constrain `τ₂` on a ±9.6 ms grid.
* **Point 3 (H4's 1.25 depends on an unjustified weight) — confirmed by
  measurement.** 1.2503 (full proper volume) vs **1.2860** (created volume) vs
  1.2417 (flat). My own docstring named the created volume and my code computed
  the full one — the same "word versus execution" gap as the harness defects.
* **Point 4 (report vs preregistration) — fixed.** The frozen T8 (false-alarm
  rate on 100 noise realisations) is now implemented as written and passes with
  0/100; the Fisher condition number is renamed T8b and labelled post-hoc.
  Additionally found: `analysis.py`'s delay column is **bit-identical** to its
  `t_c` column, so its 4×4 Fisher matrix is exactly rank-deficient and its
  reported 1.8·10²⁰ is a float64 roundoff floor — my recomputation gives 3.1·10²⁰
  on the same bytes, i.e. the value does not reproduce, which is why it is
  reported as a restatement of the exact degeneracy rather than a measurement.

## A property of the metric that had to be stated, not just noticed

The **absorbed fraction is not invariant** under the choice of delay reference
`t_ref` (0.9999967 at `t_ref = t(f_min)` versus 0.9530 at `t_ref = t(f_max)`),
while the **residual SNR is invariant** to 2.3·10⁻¹⁰. So the absorbed fraction
must never be quoted as the evidence for H3; only the residual may. This is now a
threshold in `pn_orders.py` (`H_absorbed_fraction_NOT_reference_invariant`),
because a metric that moves when the physics does not is a trap.

## Independent verification of this round

`verify_new_results.py` recomputes every new headline number by a different
route: the PN identity **symbolically** with sympy; the projection via **normal
equations** instead of QR; the H4 weights by **symbolic integration** instead of
quad; the ET/CE scaling exactness and the direct-fit flatness from the recorded
numbers. Then `--selftest` corrupts each of the four artifacts in turn and
requires the corresponding verifier to go red: **all four did** (`SELFTEST_PASS`).
The new Lean file `HomogeneousUniverse.lean` compiles with `NO_SORRY_NO_AXIOM`,
audits to only the three standard axioms, and carries three negative controls
that all fail as required.

---

# Round 224 — the factor-~3 residue closed with a real IMR template

The ceiling declared in §15.2 ("without an IMR waveform the residue cannot be
separated") is lifted: phenomxpy 2.0.3 is installed. `h3_imr_check.py` replaces
the leading-order inspiral with a full `IMRPhenomT` on the pipeline's OWN grid
(df = 0.25 Hz) and in the SAME discrete convention.

## What the measurement settled

* **The IMR template is 15 % lower** than the inspiral: 37.3118 -> **31.6595**
  (H1), 33.7097 -> **28.7350** (L1) in 20-300 Hz.
* **But that is not enough to reach ~20.** The threshold "within 25 % of 20"
  **FAILED** (31.66 vs 20 is 58 %). Recorded as a failed threshold, not hidden.
* **The data's preferred amplitude scale moved to 0.5288 (H1) / 0.4470 (L1)**,
  from 0.1949 / 0.1655. So **more than half** of the old 5.1x deficit was
  leading-order template imperfection (B4) — this is the quantitative naming of
  B4 that §15.5 item 2 was missing.
* **The residue is orientation, and it is measured.** rho_opt falls monotonically
  with inclination: 31.66 (0 deg) -> 19.79 (60 deg) -> 15.83 (90 deg), crossing
  rho_opt = 20 at **iota ~ 59.19 deg**, inside the physically allowed range.

## Three errors of this round

| # | Error | How it showed up | Caught by |
|---|---|---|---|
| 1 | multiplied the inspiral amplitude by `fs` and did not divide by `fs` in the matched filter | NC1 gave rho_opt = 297708 instead of 37.31 — off by 8000x | **NC1** |
| 2 | NC5 asserted "the model is strictly zero below its own f_min" | **false**: phenomxpy carries a tail to 2.6e-23; measured 0.52 % of band power | the control itself (it did not go red where it should have) |
| 3 | the verifier built the model at df = 0.5 and interpolated linearly | **19 %** of the band power destroyed | the 2.2 % disagreement with the audit |

**The third is the useful one.** The FD phase rotates by ~pi between adjacent
0.5 Hz samples, so linear interpolation of a complex oscillating function is not
interpolation. Measured: at SHARED frequencies the two grids agree to **1e-16**,
while after interpolation the power ratio is **0.810**. So the **verifier was
wrong and the audit was right**, established by measurement rather than argument.
The false NC5 claim and the wrong 1e-3 threshold are **recorded in the code, not
deleted**.

## Controls vs findings, kept apart on purpose

A failed CONTROL means a broken instrument and a void run; a failed FINDING is a
result. The exit code gates on the controls only. All five controls pass:
NC1 reproduces the recorded 37.31184231798413 / 33.709713860505396 to <1e-6,
NC2 IMR rho_opt finite and positive in every band, NC3 amplitude linearity
(2.000000000), NC4 zero waveform -> 0, NC5 below-f_min power 0.52 % (< 1 %).

## Independent verification

`verify_h3_imr.py` imports nothing from the audit: own Welch, own log-log PSD
interpolant, own `IMRPhenomT` call, own quadrature. Independent IMR rho_opt
matches the recorded value at **rel = 0.00e+00** on both detectors; the crossing
inclination reproduces to 1e-6 and independently gives rho_opt = 19.9817 there.
Five single-field corruptions of the artifact all turn it red (`SELFTEST_PASS`),
and the artifact is restored byte-identically. The artifact reproduces
**byte-for-byte** on a re-run: sha256 `36ed9094…`.

## What remains open

1. **Spins = 0** and the base number at iota = 0 are stated simplifications; the
   IMR rho_opt is an upper bound, not the event's recovered SNR.
2. **The delay is still not fitted** to the strain; only the template's optimal
   SNR is measured here.
3. **ET/CE** — unchanged: the sign is undetermined without tau(t).
4. **Deriving tau(t) from the theory** — not done; the paper does not state it.

---

# Round 221 — the owner's objection about 37.3, and the τ₂ grid

The owner read §12.2 and §14.2 and made a substantive objection rather than a
style note: **the template's optimal SNR of 37.3 (leading-order inspiral,
M_c = 28 M☉, D_L = 410 Mpc) is HIGHER than the published full-IMR
single-detector SNR (~20).** Either the normalisation is wrong, or 37.3 is an
integral over a band that includes the region where the inspiral formula is
already invalid — in which case "optimal SNR 37.3" is not physically meaningful
and B2 is B3 in other clothes.

He was right about the number and wrong about the identity, and both halves came
out of the measurement.

## Four errors, and every one of them is in what was being COMPARED, not in code

| # | Error | How it showed up | Caught by |
|---|---|---|---|
| 1 | `verify_rho_opt_bands.py` interpolated the PSD **linearly**; the audited pipeline interpolates it **log-log** (`psd_on_grid`) | a residual of **1.7e-3 on every band**, including the smooth ones — a systematic, not a band-edge effect | the convergence test refused to converge; measuring the two interpolants directly gave 2.8e-3 (linear) against 8.4e-4 (log-log) |
| 2 | The "convergence" test refined `df` on a band whose top edge **moved**, because f_ISCO = 67.6304 Hz is not a multiple of 0.25 Hz | a spurious non-monotone 6.3e-4 → 2.0e-4 → 1.4e-4 | the test's own numbers: a convergence test on a moving interval is not a convergence test |
| 3 | Two of my first-draft criteria were simply **wrong**: "ISCO band < 20" (it is 20.71) and "IMR band within 8 of 20" (it is 36.58) | two false red cells | the measured values themselves. The wrong criteria are **recorded in the code**, not deleted, because the second one is the finding: cutting at 250 Hz does *not* fix the number, which is what makes the reduction a validity effect rather than a bandwidth effect |
| 4 | `h3_grid_bounds.py` NC1 compared the coarse-grid peak against the recorded `peak_snr` 7.2739 instead of the recorded `profile_at_zero` 7.2732 | a false red on a control that should pass | the comparison itself: 7.2739 comes from a two-stage FINE search over `t_c` and is a different quantity |

**The pattern across all four:** the arithmetic was right and the *identity of the
thing being compared* was wrong. That is the same failure mode the owner has been
finding in my work for several rounds — and this time two of the four were caught
by controls I had written for exactly that purpose, which is the first time the
controls fired on me rather than the reverse.

## What the measurement settled

* **37.3 is not a normalisation error.** `h3_rho_opt_bands.py` reproduces
  37.31184231798413 exactly (rel < 1e-6) in the pipeline's own convention.
* **It is a band integral, and the owner's reading is correct.** 69–71 % of the
  template's σ² lies above the ISCO frequency. Cut at ISCO: **20.71 / 18.17** —
  the published value. Cut at 250 Hz: 36.58 / 32.87 — so the reduction is a
  validity effect, not a bandwidth effect.
* **But B2 is not B3 in other clothes.** The decisive test: if it were, cutting
  the band would drive the data's preferred amplitude scale to 1. Measured
  0.195 → 0.339 (H1) and 0.166 → 0.291 (L1): the deficit shrinks from 5.1× to
  3.0× and **stops**. B3 accounts for ~1.7 of the 5.1 and leaves a real factor ~3
  unexplained. That residue is the honest open item.
* **The grid item closed both ways.** Expanded ±9.6 → ±38.4 ms: still
  unconstrained. And a physics bound that no scan would have found: a quadratic
  delay must sweep less than π of phase across the band, capping the amplitude at
  **1.667 ms** — the paper's 1.2 ms uses **72 %** of that budget.
* **And one thing that had to be said rather than left in the artifact:** on H1
  the profile maximum sits AT the grid edge (−38.4 ms), i.e. it is still rising
  when the scan stops. That is not a best fit and not a bound — there is no
  interior maximum at all, which is the strongest form of "unconstrained".

## Independent verification of this round

`verify_rho_opt_bands.py` builds its **own** PSD (own Welch code, own log-log
interpolant), computes f_ISCO by a **different algebra** (total mass in SI rather
than the chirp-mass formula), and uses a **different quadrature** (continuous
integral rather than discrete sum). Its primary check — an independently built
discrete sum against the recorded number — agrees at **rel = 0.00e+00 on all six
cells**. Four single-field corruptions of the artifact all turn it red
(`SELFTEST_PASS`). Both new artifacts reproduce **byte-for-byte** on a re-run.


---

## Ход 226 — ошибки и ложные критерии

### Пять ошибок, все пойманы контролем

1. **`np.fft.irfft` вместо `np.fft.ifft`** в матч-фильтре. `irfft` предполагает
   сопряжённо-симметричный спектр и возвращает `2·Re(Σ)`, то есть другую
   величину. Пик занижался ровно в 16384 раза (= fs), и инжекция 20 мс
   «восстанавливалась» на краю сетки. Поймано **NC6** (сверка с явной матрицей).
2. **`B = np.zeros(seg_n)` проиндексирован маской длины rfft** (`seg_n//2+1`).
   NumPy не возражает против маски короче массива — он молча берёт префикс,
   поэтому ошибка была тихой: все профили ~10⁵ меньше, инжекция всегда на краю.
   Поймано **NC6**.
3. **`template_time_series` делил на `fs`.** Безшумный круговой прогон давал
   σ/fs = 0.0019 вместо σ = 31.66. Поймано **NC7** (круговой прогон с целевым
   С/Ш 25 возвращал 0.0015).
4. **Сетка ±9.6 мс для инжекции 20 мс.** «Восстановлено 9.6 мс» — это край
   сетки, а не результат. Поймано **NC4** (добавлена проверка
   `recovered_at_scan_edge`).
5. **Верификатор: двухступенчатый поиск пика.** Коарс-проход по stride·8, затем
   уточнение в окне ±16. На H1 **тихо пропустил истинный максимум**: вернул
   7.8106 там, где правильное значение 8.2011 (5 %). Оптимизация скорости,
   которая меняет ответ, — не оптимизация. Удалена; быстрый путь теперь точный
   (обратный ДПФ), матричный оставлен точечной сверкой.

### Ложные критерии, записанные, а не удалённые

**NC6 первого черновика** (`h3_imr_spins.py`): «шаблон с большей полосовой
мощностью не может иметь меньший пик матч-фильтра на тех же данных».
**Измерено: ложно.** Опубликованные спины поднимают ρ_opt 31.66 → 34.84, а пик
фильтра **падает** 16.74 → 13.51: спин меняет и **фазу** шаблона, поэтому
перекрытие с данными не есть функция одной мощности. Критерий заменён на
неравенство Коши — Шварца (пик ≤ ρ_opt), которое истинно и проверяемо.

**NC1b** — не ложный критерий, а самый сильный контроль: LO+pn+расширенная
подгонка обязана воспроизвести записанные числа `h3_grid_bounds.json`. Без него
новые IMR-числа несравнимы со старыми.

### Ошибки верификатора, найденные измерением

- **V2** сначала сравнивал с **неограниченными** углами (|a1|, |a2| на своих
  границах), которые артефакт намеренно исключает по χ_eff. Расхождение 12–19 %.
  Исправлено: верификатор строит свой конверт с ограничением χ_eff.
- **V5/V6** сначала использовал окно sqrt-hann вместо tukey 0.25 и расходился
  на 1.6 % по падению профиля. Измерено: `Σw²` = 32767 (sqrt-hann) против 55295
  (tukey 0.25). Окно — часть **определения** пайплайна; независимость даёт
  квадратура, а не окно.

### Числа, которых не было в проекте до этого хода

- Опубликованные однометекторные С/Ш: **19.5 (H1)** и **13.3 (L1)**, не «~20».
- 24 и 25.1 — **сетевые** величины; сравнение с ними однометекторного ρ_opt —
  ошибка категории.
- Отношение опубликованных С/Ш 1.466 против отношения оптимальных С/Ш шаблона
  1.102 и извлечённых фильтром 1.303; остаток — отклик антенн.
- Множитель наклонения равен **(1+cos²ι)/2 точно** и не зависит от спинов.

---

# Хронология ходов 216–227 (перенесена из REPORT.md, ход 227)

Это хроника: что делалось в каком ходе и что при этом сломалось. Текущее
состояние результата — в `REPORT.md`; разметка утверждений — в `CLAIMS.md`.

### 12. Ход 220: шаг (а) — высшие PN-порядки, ET/CE, и инстанцирование критерия

### 12.1 Вырождение НЕ является артефактом ведущего порядка

Ваш пункт 1 требовал показать, что Lean-критерий действительно **инстанцирован**
конкретными физическими векторами, а не связан со статьёй словесно. Это
делается двумя **точными тождествами**, справедливыми при **любом** PN-порядке
для любой гладкой фазы `Ψ_chirp = F(v, η)`, `v = (π M f)^{1/3}`:

```
(I)   ∂Ψ/∂t_c = 2π f                                  (t_c входит линейно)
(II)  f ∂Ψ_chirp/∂f = M_c ∂Ψ_chirp/∂M_c               (теорема Эйлера по v)
```

`(II)` верно потому, что `Ψ_chirp` зависит от `M_c` и `f` **только через `v`**,
а `v` однородна степени 1 по `(M_c^{1/3}, f^{1/3})`. **Ни один PN-коэффициент в
это не входит**, поэтому никакой порядок усечения его не нарушает. Следствие:

```
ΔΨ_лин = 2π f τ̇ (t(f) − t_ref)
       = τ̇ [ M_c ∂Ψ/∂M_c + (t_c − t_ref) ∂Ψ/∂t_c ]
```

— **ровно** комбинация направлений параметров источника, при любом PN-порядке.

**Проверено.** `pn_orders.py` → exit 0. Тождество `(II)` проверено **символьно**
(sympy, усечённый ряд) с **негативным контролем**: добавка `ln(M_c/f)`, не
являющаяся функцией `v`, **ломает** тождество (а `ln(M_c f)` — не ломает, потому
что `M_c f` есть функция `v`; это тоже проверено, чтобы контроль не был
вакуумным). Численно на полном 3.5PN: относительное расхождение `3.0·10⁻⁹`.

| базис | const поглощено | linear поглощено | квадратичный остаток С/Ш |
|---|---:|---:|---:|
| ведущий порядок, `(M_c,t_c,φ_c)` | 1 − 6.7·10⁻¹⁶ | 1 − 2.2·10⁻¹⁶ | **0.0872** |
| **3.5PN**, `(M_c,t_c,φ_c)` | 1 − 1.0·10⁻¹⁵ | 1 − 5.6·10⁻¹⁶ | 0.0551 |
| **3.5PN**, `(M_c,η,t_c,φ_c)` | 1 − 1.0·10⁻¹⁵ | 1 − 5.6·10⁻¹⁶ | **0.0551** |

**Вывод: H1/H2 — теорема, а не приближение ведущего порядка.** Численное
значение квадратичного остатка меняется (0.087 → 0.055 при 3.5PN), но вывод
«постоянная и линейная задержки поглощены точно» не меняется вообще.

**Найденное здесь свойство метрики.** Доля поглощения **не** инвариантна
относительно выбора точки отсчёта `t_ref` (0.9999967 при `t_ref = t(f_min)`
против 0.9530 при `t_ref = t(f_max)`), тогда как **остаточный С/Ш инвариантен**
до 2.3·10⁻¹⁰. Значит **долю поглощения нельзя приводить как доказательство** —
только остаточный С/Ш. Это записано в пороги как `H_absorbed_fraction_NOT_reference_invariant`.

### 12.2 ET/CE: что можно сказать, и где предел

Реальной кривой ET-D / CE в контейнере нет (сайт ET и файлы LALSuite недоступны).
Цитируемое: ET нацелен на «~factor of 10 compared to the existing detectors» в
полосе от нескольких Гц (источник `src_e8a31a47672c`). Поэтому `et_ce_scaling.py`
делает **не** прогноз для ET, а исследование чувствительности **статистики** по
двум ручкам: масштаб шума и нижняя граница полосы.

| конвенция | 20 Гц, aLIGO | 20 Гц, ×10 | 5 Гц, aLIGO | 5 Гц, ×20 |
|---|---:|---:|---:|---:|
| амплитуда привязана к 1.2 мс на краю полосы | 0.087 | 0.872 | 0.0047 | 0.094 |
| `τ₂` зафиксирован | 0.087 | 0.872 | 7.625 | 152.5 |

* Постоянная и линейная задержки поглощены **точно** во **всех** клетках — это
  тождество, а не число.
* Остаток квадратичного члена масштабируется **ровно линейно** по
  чувствительности (отношение 10.000000000 при ×10, 20.000000000 при ×20).
* **Знак выгоды от широкой полосы зависит от конвенции** и переворачивается:
  при фиксированной амплитуде — падает (0.087→0.0047), при фиксированном `τ₂` —
  растёт (0.087→7.625). Статья **не задаёт** функцию `τ(t)`, поэтому знак
  ET/CE-выгоды ею не определён. Это честный потолок шага (а).

---

### 13. Ход 220: шаг (в) — однородная вселенная как теорема

Статья сама признаёт (§6), что в однородной вселенной эффект неотличим от
тёмной энергии. `work/lean/HomogeneousUniverse.lean` переводит это признание в
теорему в той же аддитивной форме, что и `Identifiability.lean`:

наблюдаемое `obs H δ a = H a + δ a`, где `H` — базовое расширение, `δ` — вклад
«созданного времени»; измеряется только `obs`.

| теорема | содержание |
|---|---|
| `creation_not_identifiable_unpinned` | если базовое расширение `H` **свободно** (как при подгонке тёмной энергии), то `δ` **не** идентифицируем: любое `δ` воспроизводится сдвинутым базисом |
| `creation_identifiable_pinned` | если `H` зафиксировано независимо, `δ` **идентифицируем** |
| `baseline_unique_when_pinned` | при зафиксированном `H` согласие наблюдаемых влечёт согласие `δ` — фиксация и есть та линия, что разделяет ветви |
| `demo_unpinned_pair` | **анти-вакуумность**: конкретная пара `δ ≠ δ'` с разными базисами даёт **то же** наблюдаемое |
| `demo_pinned_identifiable` | вторая ветвь тоже непуста |

**Проверка.** `sh verify_lean_cosmo.sh` → rc=0, `NO_SORRY_NO_AXIOM`, аудит
аксиом — только `[propext, Classical.choice, Quot.sound]`.
**Три негативных контроля, все краснеют:** снять финальный шаг у
pinned-теоремы → `unsolved goals`; заявить **обратное** первой теореме → ошибка
компиляции; подсунуть `sorry` → скан ловит.

**Смысл, и граница.** Это теорема о **идентифицируемости**, а не о физике
тёмной энергии: она говорит, что из истории расширения нельзя выделить вклад
«созданного времени», **пока базовое расширение подгоняется** — то есть ровно
то, что признаёт статья. Она ничего не говорит о том, верна ли теория.

---

### 14. Ход 220: разбор ваших четырёх замечаний

### 14.1 Что именно предсказывает теория — согласен, это остаётся открытым

Ваше замечание принято и **не** закрыто. Lean доказывает критерий для
абстрактной линейной модели, а расчёт применяет его к фазе ведущего порядка;
связь с **конкретной** задержкой статьи не выведена, потому что статья называет
свой расчёт эвристическим и **не задаёт `τ(t)`**. Что сделано в этом ходе:
§12.1 показывает, что вырождение **не** артефакт усечения — оно следует из
симметрии `v`, а не из числа членов. Что остаётся: **вывести `τ(t)` из теории**
либо сузить итоговое утверждение. Я формулирую его теперь так: *доказано
вырождение постоянного и линейного членов в принятой линейной модели, а не
невозможность проверить всю «now»-теорию*.

### 14.2 H3 — зазор 37.3 ↔ 7.27 объяснён и разделён на части

`h3_direct_fit.py` → exit 0. Зазор между оптимальным С/Ш шаблона (37.3 H1,
33.7 L1) и пиком фильтра (7.27 / 5.58) разложен **измерением**:

| компонента | H1 | L1 | вердикт |
|---|---:|---:|---|
| B1 калибровка PSD (инжекция с целевым 20) | 20.109 | 20.062 | **исключена** — фильтр не смещён |
| B2 предпочитаемый данными масштаб амплитуды | 0.195 | 0.166 | **остаётся**: данные в 5.1× слабее модельной амплитуды |
| B3 доля σ² выше ISCO (67.6 Гц) | 0.692 | 0.709 | **измерена**: 69 % мощности шаблона лежит там, где инспиральный шаблон недействителен |

То есть зазор — это **не** ошибка калибровки (B1 исключена) и не загадка: он
распадается на (B2) несовпадение амплитуды и (B3) усечение полосы. **B4**
(несовершенство фазы шаблона) отделить от B3 на этих данных нельзя — это
названо, а не спрятано.

**Прямая подгонка задержки к сигналу — сделана.** Сканирование `τ₂` по сетке
±9.6 мс амплитуды, на каждом `τ₂` — максимум матч-фильтра по `(t_c, φ_c)`:

| | H1 | L1 |
|---|---:|---:|
| профиль при `τ₂ = 0` | 7.2732 | 5.5779 |
| профиль при `τ₂` статьи (1.2 мс) | 7.2614 | 5.5795 |
| падение профиля на всей сетке | **0.213** | **0.046** |
| ограничен ли `τ₂` на сетке | **нет** | **нет** |

**Результат: данные не ограничивают `τ₂` вообще** — профиль не падает на 1σ
(падение 0.21 и 0.05 при пороге 1.0) на сетке шириной ±9.6 мс, и значения при
`τ₂ = 0` и при 1.2 мс неразличимы (разница 0.012 и 0.002). Это **прямое**
измерение того же вывода, что даёт проекция, но полученное подгонкой к записи, а
не проекцией направления. Отсюда честная формулировка: не «верхний предел
34.5 мс», а «на сетке ±9.6 мс данные не предпочитают ни одного значения».

### 14.3 H4 — вес: ваше замечание подтверждено измерением

`h4_weight_audit.py` → exit 0. Три веса на **одной** заявленной оболочке
`[1.5, 4] Rₛ(62)`:

| вес | среднее `(1+z)` | отличие от 1.25 |
|---|---:|---:|
| **полный собственный объём** `4πr²/√(1−Rₛ/r)` | **1.250339** | **0.03 %** |
| созданный объём `4πr²(1/√(1−Rₛ/r) − 1)` | **1.285956** | 2.9 % |
| плоский (координатный) объём `4πr²` | 1.241729 | 0.66 % |

**Вы правы.** Совпадение с 1.25 — свойство **одного** (естественного, но
незаявленного) веса. Вес «созданного объёма», который называет докстринг самого
`h4_paper_text_audit.py`, даёт **1.286**, а не 1.25. Претензия ослаблена:
*«1.25 воспроизводится естественным взвешиванием»*, а не *«1.25 выведено»*.
Обратите внимание: мой `h4_paper_text_audit.py` **объяснял** вес как созданный
объём, а **считал** по полному — то же расхождение «слово против исполнения»,
что и в дефектах harness. Три веса расходятся на 0.044 (3.5 % от 1.25), и это
записано как негативный контроль.

### 14.4 Пререгистрация — расхождение T8 устранено

`prereg_t8.py` → exit 0. Вы правы: замороженный T8 — «доля ложных срабатываний
на чистом шуме (100 реализаций) ≤ 5 %», а `analysis.py:409` реализовал под этим
именем число обусловленности Фишера. Это разные величины, и первая не
исполнялась.

**Теперь исполнена как записана.** Статистика `z = ⟨белый шум, единичное
непоглощённое направление задержки⟩`; при чистом шуме `z ~ N(0,1)` **точно**.
100 реализаций при пороге обнаружения 5:

| | значение |
|---|---:|
| `z` сигнала (задержка 1.2 мс) | **0.0872** |
| среднее `z` шума | 0.049 |
| σ(`z`) шума | 0.933 |
| max\|`z`\| шума | 2.13 |
| ложных срабатываний | **0 / 100** |
| **T8 по пререгистрации** | **пройден** |

**Но T8 сам по себе недостаточен, и это важно сказать.** Он говорит о
**детекторе**, а не об эффекте: чистый шум никогда не достигает `z = 5`, а
задержка статьи даёт `z = 0.087` — в **57 раз ниже того же порога**. T8 был бы
выполнен и прибором, который вообще не видит сигнал. Поэтому `z` сигнала
приведён рядом.

**T8b переименован.** Число обусловленности Фишера — **post-hoc** проверка,
добавленная после заморозки; она приводится отдельно и **никогда** не как
доказательство по замороженному порогу. Дополнительно найдено: колонка задержки
`h·2πf` в `analysis.py` **побитово совпадает** с колонкой `t_c` — что и есть
содержание H1, — поэтому та матрица Фишера **точно вырождена**, а число 1.8·10²⁰
есть дно округления float64, а не физическая величина. Мой пересчёт даёт
3.1·10²⁰ на тех же байтах: **roundoff-значение не воспроизводится**, и именно
поэтому оно приводится как пересказ точного вырождения, а не как измерение.

### 14.5 Устаревшее «`h5py` нет» — исправлено

§7 п.3 и §10 уже описывали использование `h5py`; строка о его отсутствии в §10
заменена (см. выше). `h5py` 3.16.0 доступен и использован.

---

### 15. Ход 221: четыре открытых вопроса, закрытых измерением

Владелец прочитал §12.2 и §14.2 и выдвинул возражение по существу: **оптимальный
С/Ш 37.3 для ведущего порядка при `M_c = 28 M_⊙`, `D_L = 410 Мпк` выше, чем
полный IMR-шаблон даёт (20).** Либо нормировка, либо 37.3 — интеграл по полосе,
включающей область, где инспиральная формула недействительна. И если так, то
«оптимальный С/Ш 37.3» физически не осмыслен, а B2 — это B3 в другой одежде.
Плюс отдельно: сетка ±9.6 мс.

### 15.1 Откуда 37.3 (`h3_rho_opt_bands.py`, `art_3f466cf5f33f`)

`rho_opt(band) = sqrt(C·Σ_{k∈band}|H_k|²/S_k)` в **той же** дискретной конвенции,
что и пайплайн, поэтому число прямо сравнимо с записанным.

| полоса | H1 | L1 | доля |
|---|---:|---:|---:|
| 20–300 Гц | **37.3118** | **33.7097** | 1.000 |
| 20–250 Гц | 36.5850 | 32.8747 | 0.981 / 0.975 |
| **20–67.63 Гц (ISCO)** | **20.7098** | **18.1707** | **0.555 / 0.539** |

**NC1 — полная полоса воспроизводит 37.31184231798413 точно** (rel < 1e-6).
Значит это **не** ошибка нормировки. Это интеграл по 20–300 Гц, и **69–71 % σ²
лежит выше ISCO** (0.6919 / 0.7094 — совпадает с прежними 0.692 / 0.709).
Обрезка на ISCO даёт **20.71 / 18.17** — опубликованное значение. Обрезка на
250 Гц даёт 36.58 / 32.87, то есть снижение — эффект **области недействительности**,
а не ширины полосы. **Владелец прав: 37.3 — не физически осмысленное число.**

### 15.2 Но B2 — не B3 в другой одежде

Если бы «данные в 5.1× слабее» было тем же, что «шаблон заявляет мощность выше
ISCO», обрезка увела бы предпочитаемый масштаб к 1:

| полоса | масштаб H1 | масштаб L1 |
|---|---:|---:|
| 20–300 Гц | **0.1948** | **0.1657** |
| 20–250 Гц | 0.2183 | 0.1578 |
| 20–67.63 Гц | **0.3385** | **0.2909** |

Дефицит сжимается с 5.1× до **3.0× / 3.4×** и **останавливается**. B3 объясняет
≈фактор 1.7 из 5.1 и оставляет **подлинный фактор ~3**. Это и есть разделение
B3/B4 **измерением**: не «отделить нельзя», а «отделено количественно». Остаток
кандидатно объясняется откликом события vs усреднением по небу и тем, что «20» —
число полного IMR; **без IMR-формы не разделяется.** Это главный открытый пункт.

### 15.3 Сетка τ₂ — расширена и ограничена физикой (`art_2343189a11d1`)

**Расширение:** ±9.6 → **±38.4 мс** (вчетверо, 129 точек). Падение профиля
**0.791 / 0.461** при пороге 1σ — **всё ещё не ограничено**. Контроли: грубый пик
при τ₂ = 0 совпадает с записанным `profile_at_zero` (7.273249726912924 /
5.577882654444806) до 1e-6.

**И то, что нельзя оставить в артефакте:** у H1 максимум профиля стоит **на краю
сетки** (−38.4 мс) — профиль ещё **растёт**, когда скан останавливается. Это не
лучшее значение и не граница: внутреннего максимума нет вовсе. Самая сильная форма
«не ограничено».

**Физическая граница** (не зависит от скана): фаза `2πf·τ₂·(t(f)−t_ref)²` не
должна заметать больше π по полосе.

| | значение |
|---|---:|
| τ₂ по критерию фазы | 0.0023370481 /с² |
| **амплитуда** | **1.6667 мс** |
| потолок по длительности сегмента | 4000 мс |
| статья | 1.2 мс |
| **доля бюджета** | **72.0 %** |

Максимум ядра — на верхнем краю полосы (300 Гц); пересчёт другим сканом
(2·10⁶ точек) даёт то же τ₂ до 1e-6. **Смысл:** число статьи не просто не
ограничено данными — оно близко к краю того, что вообще можно назвать малым
возмущением этого шаблона.

### 15.4 Независимая проверка и мои ошибки

`verify_rho_opt_bands.py` (`art_45a21834782d`) не импортирует аудит: свой PSD
(свой Уэлч, свой log-log интерполянт), свой путь для `f_ISCO` (через полную массу
в СИ), другая квадратура (непрерывный интеграл).

| проверка | результат |
|---|---|
| независимая дискретная сумма == записанной | **rel = 0.00e+00** во всех шести ячейках |
| `f_ISCO` другим путём | 67.630422 == 67.630422 |
| непрерывный vs дискретный | в пределах 3e-3 на всех полосах |
| сходимость при df→0 | 8.4e-4 → 2.8e-5 → 2.9e-6, монотонно |
| **селфтест** | **4/4 порчи краснеют** |

**Четыре моих ошибки, все пойманы контролями, ни одна не в коде — все в
утверждении о том, что сравнивается:**

1. верификатор интерполировал PSD **линейно**, аудит — **log-log**; расхождение
   1.7e-3 на каждой полосе (измерено: два интерполянта дают 2.8e-3 против 8.4e-4);
2. «сходимость» на **движущейся** полосе (f_ISCO не кратно 0.25 Гц) — ложная
   немонотонность;
3. два первых порога **неверны**: «ISCO-полоса < 20» (она 20.71) и «IMR-полоса в
   пределах 8 от 20» (она 36.58) — критерии переписаны, а неверные **записаны**,
   а не удалены;
4. NC1 в `h3_grid_bounds.py` сравнивал с `peak_snr` вместо `profile_at_zero` —
   ложное красное (7.2739 против 7.2732).

### 15.5 Что осталось открытым

1. **Остаток фактора ~3** в зазоре после снятия B1 и B3 — главный открытый пункт.
2. **B4** отделён от B3 количественно, но **не поименован**.
3. **ET/CE** — знак по-прежнему не определён без τ(t); не трогал.
4. **Вывести τ(t) из теории** — по-прежнему не сделано.

---

### 16. Ход 224: остаток фактора ~3 закрыт настоящим IMR-шаблоном

Потолок, объявленный в §15.2 («без IMR-волновой формы остаток не разделяется»),
снят: `phenomxpy` 2.0.3 установлен. `h3_imr_check.py` (`art_...`) заменяет
инспиральный шаблон полным `IMRPhenomT` на **той же сетке df = 0.25 Гц**, что и
пайплайн, в **той же дискретной конвенции** `rho_opt² = 4·df·Σ|h~|²/Sn`.

### 16.1 Что показало измерение

| полоса | H1 LO | H1 IMR | L1 LO | L1 IMR |
|---|---:|---:|---:|---:|
| 20–300 Гц | 37.3118 | **31.6595** | 33.7097 | **28.7350** |
| 20–250 Гц | 36.5850 | 30.9531 | 32.8747 | 27.9442 |
| 20–67.63 Гц | 20.7098 | 16.7296 | 18.1707 | 14.6924 |

**IMR-шаблон ниже инспирального на 15 %** — но этого мало. Замена шаблона
**не** приводит число к опубликованным ~20: порог «в пределах 25 % от 20»
**не выполнен** (31.66 против 20 — это 58 %). Это записано как **невыполненный
порог**, а не спрятано.

**Зато предпочитаемый данными масштаб амплитуды сдвинулся к 1:**
0.1949 → **0.5288** (H1), 0.1655 → **0.4470** (L1). То есть **больше половины**
прежнего дефицита 5.1× был дефектом шаблона (B4), а не данных и не PSD.

### 16.2 Остаток — это ориентация, и он измерен

`inclination = 0` (face-on) — максимальная амплитуда. Развёртка по наклонению:

| ι, ° | 0 | 15 | 30 | 45 | 60 | 90 |
|---|---:|---:|---:|---:|---:|---:|
| ρ_opt (H1) | 31.66 | 30.60 | 27.70 | 23.74 | **19.79** | 15.83 |

**ρ_opt = 20 достигается при ι ≈ 59.2°** — внутри физически допустимого
диапазона. Это переводит остаток из «необъяснённого фактора ~3» в
«ориентация плюс определение опубликованного числа», причём количественно.
Независимая проверка на этой же точке даёт 19.98.

### 16.3 Контроли и три моих ошибки

Контроли (все проходят, exit 0): NC1 инспиральный rho_opt воспроизводит
записанное 37.31184231798413 / 33.709713860505396 **точно** (rel < 1e-6);
NC2 IMR rho_opt конечен и > 0 во всех полосах; NC3 линейность по амплитуде
(2.000000000); NC4 нулевая волновая форма → ρ = 0; NC5 мощность ниже f_min
модели — 0.52 % от полосной, ниже 1 %.

| # | Ошибка | Как проявилась | Кем поймана |
|---|---|---|---|
| 1 | умножал инспиральную амплитуду на `fs`, а в матч-фильтре не делил на `fs` | NC1 дал 297708 вместо 37.31 — расхождение в 8000× | NC1 |
| 2 | NC5 утверждал «модель строго нулевая ниже своего f_min» | **ложно**: phenomxpy даёт хвост 2.6·10⁻²³; замер — 0.52 % мощности | сам контроль |
| 3 | верификатор строил модель на df = 0.5 и интерполировал линейно | терялось **19 %** полосной мощности: фаза вращается на ~π между соседними отсчётами, линейная интерполяция — не интерполяция | расхождение 2.2 % с аудитом |

Третья — самая полезная: **верификатор ошибся, а аудит был прав**, и это
установлено измерением (на общих частотах сетки совпадают до 1e-16, а после
интерполяции отношение мощностей 0.810). Пороги, оказавшиеся неверными,
**записаны в коде, а не удалены**.

### 16.4 Независимая проверка

`verify_h3_imr.py` не импортирует аудит: свой Уэлч, свой log-log интерполянт
PSD, свой вызов `IMRPhenomT`, своя квадратура.

| проверка | результат |
|---|---|
| независимый IMR rho_opt == записанному | **rel = 0.00e+00** (H1 и L1) |
| развёртка по наклонению монотонна | да |
| точка пересечения ρ = 20 | 59.19340049167065 == записанной |
| ρ при этой точке | 19.9817 (ожидалось ~20) |
| селфтест: 5 порч артефакта | **5/5 краснеют**, артефакт восстановлен байт-в-байт |

Артефакт воспроизводится **байт-в-байт** при повторном прогоне
(sha256 `36ed9094…`).

### 16.5 Что осталось открытым

1. **Спины = 0** и `ι = 0` в базовых числах — упрощение, названное в артефакте.
2. **Задержка по-прежнему не подогнана** к записи; здесь измерен только
   оптимальный С/Ш шаблона.
3. **ET/CE** — без изменений: знак не определён без τ(t).
4. **Вывести τ(t) из теории** — не сделано.


---

### 17. Ход 226: два пункта, названные в §16.5, закрыты

Вход владельца (msg_00226): (1) прогнать IMRPhenomT **со спинами**, чтобы
формулировка «верхняя граница» сменилась на «восстановленный С/Ш»; (2) подогнать
задержку τ₂ **к записи** полным IMR-шаблоном вместо инспирала — «это то, что
нужно для paper».

### 17.1 Пункт 1 — спины и наклонение (`h3_imr_spins.py`, `art_d8f0b19dced3`)

**Первое, что пришлось исправить в постановке.** Проект с хода 216 сравнивал
свои числа с «опубликованным ~20». Такого числа нет ни в одной статье:

| величина | значение | источник |
|---|---:|---|
| re-weighted С/Ш H1 | **19.5** | arXiv:1602.03839 |
| re-weighted С/Ш L1 | **13.3** | arXiv:1602.03839 |
| сетевой С/Ш | 24 | arXiv:1602.03837 |
| сетевой оптимальный С/Ш | 25.1 ± 1.7 | arXiv:1602.03840 |

24 и 25.1 — **сетевые** величины; сравнивать с ними однометекторный ρ_opt —
ошибка категории. Цель — 19.5 и 13.3. Круглое «20» выведено из обихода.

**Второе — конверт по спинам не прямоугольник.** a1 и a2 ограничены по
отдельности, но **совместно** связаны через χ_eff = −0.07 ± 0.17. Угол
(0.69, 0.89) несёт χ_eff = +0.78, что измерение исключает. Скан идёт по
`|a1| ≤ 0.69, |a2| ≤ 0.89, χ_eff ∈ [−0.24, +0.09]`; исключённый угол оставлен
в именованных наборах, чтобы разница была видна, а не спрятана.

| набор | a1 | a2 | χ_eff | ρ_opt (H1) | ρ_opt/ρ_opt(0) |
|---|---:|---:|---:|---:|---:|
| zero | 0.00 | 0.00 | 0.000 | 31.6595 | 1.0000 |
| χ_eff median | −0.07 | −0.07 | −0.070 | 31.1073 | 0.9826 |
| medians | +0.32 | +0.44 | +0.373 | 34.8362 | 1.1003 |
| medians_low | +0.03 | +0.04 | +0.034 | 31.9401 | 1.0089 |
| **bounds (вне χ_eff)** | +0.69 | +0.89 | +0.778 | 38.5973 | 1.2191 |
| **bounds_low (вне χ_eff)** | −0.69 | −0.89 | −0.778 | 26.5692 | 0.8392 |

**Ответ.** По опубликованному конверту спинов ρ_opt лежит в **30.02 … 32.49**
(H1) и **27.16 … 29.54** (L1). Опубликованное однометекторное значение
**не достигается**: спины не объясняют остаток, и наибольший спиновой эффект
двигает число **вверх**, прочь от цели. Наклонение объясняет остаток для H1:
19.5 достигается при **ι ≈ 59.8°** (при χ_eff median). Для L1 13.3 не
достигается даже при ι = 90° (минимум 13.58) — это названо открытым пунктом,
а не сглажено.

**Побочно измеренное отношение.** Опубликованные однометекторные С/Ш относятся
как 19.5/13.3 = **1.466**, тогда как собственные оптимальные С/Ш шаблона —
как 31.66/28.74 = **1.102**, а извлечённые фильтром — как **1.303**. Разница
между 1.303 и 1.466 — это **отклик антенн двух площадок**, которого в модели
нет вовсе. Названо числом, а не отговоркой.

**Точное тождество, найденное по ходу.** Множитель наклонения равен
**(1+cos²ι)/2 ровно** и **не зависит от спинов** (проверено при трёх наборах
спинов, совпадение 3.3·10⁻¹⁶). Поэтому совместный диапазон —
это точно `spin_range × [0.5, 1]`, а не выборочная аппроксимация. Это
контроль NC8, и он же делает формулировку честной.

### 17.2 Пункт 2 — подгонка τ₂ к записи полным IMR (`h3_direct_fit_imr.py`, `art_fd571ae9919b`)

**Одно содержательное решение, и оно не косметическое.** Для τ(t) нужен момент
времени на каждой частоте. Прежняя подгонка брала инспиральную формулу
`t_PN(f) = −(5/256) M_c^{−5/3} (πf)^{−8/3}`, которая уходит в 0 с ростом f:
`t_PN(250 Гц) = −1.0 мс`. **IMR-волновая форма так не делает** — измерено:
её групповая задержка `t_IMR(f) = −(1/2π) dΨ/df` **выходит на плато**:
`t_IMR(250 Гц) = −155 мс`. Отношение `t_IMR/t_PN` идёт от 1.03 при 20 Гц до
**155** при 250 Гц. Обе конвенции посчитаны и обе приведены: выбрать одну
значило бы выбрать ответ.

| | H1, drop | L1, drop | односторонняя граница |
|---|---:|---:|---|
| LO-шаблон, ядро pn, сетка ±9.6 мс | 0.213 | 0.046 | нет |
| LO-шаблон, ядро pn, сетка ±38.4 мс | 0.791 | 0.461 | нет |
| **IMR-шаблон, ядро pn, сетка ±9.6 мс** | 0.668 | 0.506 | нет |
| **IMR-шаблон, ядро pn, сетка ±38.4 мс** | **2.925** | **2.144** | **да** |
| IMR-шаблон, ядро imr, сетка ±38.4 мс | **1.371** | **1.047** | **да** |

**Ответ рецензенту, и он расщепляется по шаблону — это и есть находка.**
(1) С **ведущим инспиралом** профиль не ограничен ни на одной сетке: drop
0.79 / 0.46 против порога 1.0. (2) С **полным IMRPhenomT** сетка ±9.6 мс
по-прежнему не ограничена (0.27 / 0.26), но на сетке вчетверо шире профиль
**падает больше чем на 1σ в одну сторону** — 2.93 (H1) и 2.14 (L1) для ядра
pn, 1.37 и 1.05 для ядра imr. Значит с правильной волновой формой утверждение
перестаёт быть «ограничений нет вовсе» и становится **односторонней границей**
примерно на −10…−13 мс (H1) и +19…+37 мс (L1). Это **первый раз**, когда
проект извлёк из данных хоть какое-то ограничение на τ₂, а не из аргумента
физичности.

**Контроль мощности — то, без чего это ничего не стоит.** Тот же фит
**восстанавливает** инжектированную задержку 20 мс (21.75 мс H1, 20.17 мс L1)
и **не восстанавливает** заявленные статьёй 1.2 мс (0.15 / −0.56 мс при
разбросе 4–5 мс). Значит «не ограничено» здесь означает «ниже порога этого
прибора», а не «фит слеп».

### 17.3 Мои ошибки этого хода — пять, все пойманы контролем

| # | Ошибка | Как проявилась | Кем поймана |
|---|---|---|---|
| 1 | `np.fft.irfft` вместо `ifft` для комплексного матч-фильтра | пик занижен в 16384 раза; инжекция 20 мс «восстанавливалась» на краю сетки | **NC6** |
| 2 | `B = zeros(seg_n)` проиндексирован маской длины rfft | тихое несовпадение длин, все профили ~10⁵ меньше | **NC6** |
| 3 | `template_time_series` делил на `fs` | безшумный круговой прогон давал σ/fs вместо σ | **NC7** |
| 4 | сетка ±9.6 мс для инжекции 20 мс | «восстановлено 9.6 мс» = край сетки как ответ | **NC4** |
| 5 | верификатор: двухступенчатый поиск пика | **тихо пропустил истинный максимум** (7.8106 вместо 8.2011) | расхождение с аудитом |

**Пятая — самая полезная.** Оптимизация скорости, которая меняет ответ, — не
оптимизация. Она удалена; быстрый путь теперь точный (обратный ДПФ), а
матричный путь оставлен точечной сверкой (V5b, совпадение 1.15·10⁻¹⁴).

**Ложный критерий первого черновика записан, а не удалён.** NC6 в первой
редакции утверждал «шаблон с большей полосовой мощностью не может иметь
меньший пик матч-фильтра». **Измерено: ложно.** Добавление опубликованных
спинов поднимает ρ_opt с 31.66 до 34.84, а пик фильтра **падает** с 16.74 до
13.51, потому что спин меняет и **фазу** шаблона. Это четвёртый неверный
критерий в проекте; как и предыдущие, он остаётся в коде.

### 17.4 Независимая проверка (`verify_h3_round226.py`)

Не импортирует ни один из двух аудитов: свой Уэлч, свой PSD-интерполянт, свой
вызов IMRPhenomT, своя квадратура (явная матрица `exp(2πift)`), своя
развёртка фазы для ядра задержки.

| проверка | результат |
|---|---|
| V1 IMR ρ_opt (H1, L1) | **rel = 0.00e+00** |
| V2 конверт по спинам с ограничением χ_eff | **rel = 0.00e+00**, 58 точек |
| V3 множитель наклонения = (1+cos²ι)/2 | 4.44·10⁻¹⁶ |
| V4 ядро `t_IMR` при 20 и 250 Гц | **rel = 0.00e+00** |
| V5 LO+pn+расширенная: пик и падение | **rel = 0.00e+00** |
| V5b обратное ДПФ против явной матрицы | 1.15·10⁻¹⁴ |
| V6 IMR+pn+расширенная: падение > 1σ | **rel = 0.00e+00**, подтверждено |
| V7 сверка с записанным `h3_grid_bounds.json` | 1.15·10⁻¹⁴ |
| **селфтест: 7 порч артефактов** | **7/7 краснеют** |

### 18. Ход 227: совместный фит H1+L1 с откликами антенн

Вход владельца (msg_00227): (1) сделать **совместный** фит H1+L1 с **одним** `τ₂`
и свободными `(M_c, t_c, φ_c)`, и **только после этого** говорить о границе;
(2) подставить отклики антенн `F₊, F×` для положения GW150914 **вместо**
сканирования по наклонению `ι`; (3) переписать §1 и §7 под текущее состояние,
обновить таблицу файлов, хронологию вынести в `NOTES.md`; (4) нарисовать
отдельный список утверждений с пометкой «доказано / измерено при таких
допущениях / гипотеза».

Пункт (4) выполнен отдельным файлом `work/CLAIMS.md` — и именно он определил
форму §1: заголовок имеет право опираться только на **[ДОКАЗАНО]**, а все
численные границы помечены **[ИЗМЕРЕНО]** с перечисленными допущениями.

### 18.1 Два структурных факта, без которых фит не собрать

**Первый — соглашение phenomxpy.** Измерено: `IMRPhenomT` возвращает поляризации
как

```
hp(ι, φ) = A(ι)·e^{+2iφ}·hp₀ ,   hc(ι, φ) = B(ι)·e^{+2iφ}·hc₀ ,   hc₀ = −i·hp₀
```

(`A = (1+cos²ι)/2`, `B = cos ι`; тождество `hc₀ = −i·hp₀` проверено, максимальное
отклонение 3.3·10⁻⁴). Поэтому штамм на детекторе

```
F₊·hp + F×·hc = e^{2iφ}·hp₀·( F₊·A − i·F×·B )
```

то есть **всё откликовое преобразование сворачивается в один комплексный
коэффициент** `C_d = F₊,d·A(ι) − i·F×,d·B(ι)`.

**Второй — отсюда же.** Так как `h_d = C_d·h` при общем `h`, когерентная
сетевая статистика имеет вид

```
ρ_coh² = | Σ_d conj(C_d)·ζ_d |² / Σ_d |C_d|²·s_d²
```

где `ζ_d` — комплексный выход матч-фильтра для волновой формы `(1,0)`, а `s_d²`
её норма. `ζ_d` **не зависит** от `C_d`, поэтому считается один раз на
`(M_c, τ₂)` на детектор, а вся сетка `(ι, ψ)` после этого бесплатна. Фит
становится действительно когерентным: отношение амплитуд и фаз двух детекторов
**предсказано** положением на небе, а не является свободной ручкой.

**Следствие, названное как находка, а не спрятанное.** `φ_c` умножает `C_d`
каждого детектора на **одну и ту же** фазу, значит он тождественно вырожден с
общей фазой и **не идентифицируем** в когерентном фите с фиксированным
положением. Он профилируется, а не цитируется.

### 18.2 Положение на небе: чтение выбрано измерением, а не доверием библиотеке

Положение берётся из **опубликованной** карты `LALInference_skymap.fits.gz`
(LOSC, P1500227). Два независимых читателя HEALPix — `healpy` и
`astropy_healpix` — согласны на положении каждого пикселя до **1.8·10⁻¹⁵ рад**.

Но сами массивы у них **разные**, и это оказалось ловушкой. FITS-таблица хранит
`PROB` двумерным массивом; `data['PROB'].reshape(-1)` (по строкам) и массив,
который возвращает `healpy.read_map`, кладут вероятность на **разные** участки
неба. Ни одной библиотеке не поверили на слово — чтение выбрано **физикой**:
опубликованное ограничение (arXiv:1602.03840) требует, чтобы источник лежал на
кольце постоянной разности времён прихода H1−L1 = 6.9 (+0.5/−0.4) мс.

| чтение | концентрация на кольце | на задержке |
|---|---:|---:|
| построчное + NESTED | **0.1635** | **−6.90 мс** ← принято |
| `healpy.read_map` + NESTED | 0.0233 | +3.70 мс |
| перемешанная карта (20 прогонов) | 0.0061 средн., 0.0064 макс. | — |

Принятое чтение — **25×** нуля перемешанной карты. Попутно измерено, что
проекция требует вращения GMST (карта в небесных координатах, вершины
детекторов — в земной системе): без вращения пик стоит на −3.4 мс с
концентрацией 0.026.

**Проверки положения.** ML-положение RA = 134.80°, Dec = −69.79°; 90 % площадь
**616.4 deg²** (опубликовано 610), 50 % — 149.1 deg² (опубликовано 150);
`τ_H1 − τ_L1 = −6.898 мс` (опубликовано 6.9). И решающая: оптимальные С/Ш при
`ι = 0` равны 25.31 (H1) и 17.16 (L1), их отношение **1.4749** против
опубликованного однометекторного **19.5/13.3 = 1.4662** — совпадение 0.6 %.
Это и есть проверка, что положение на небе — то самое.

### 18.3 Что дал совместный фит

| ядро задержки | сетка | падение профиля | односторонняя граница |
|---|---|---:|---|
| pn | ±9.6 мс | 1.375σ | нет |
| pn | **±38.4 мс** | **3.496σ** | **да**: +3.6 мс |
| imr | ±9.6 мс | 0.852σ | нет |
| imr | **±38.4 мс** | **2.352σ** | **да**: +2.4 мс |

Пик сетевого С/Ш: **22.93** (ядро pn) и **23.07** (ядро imr); при `τ₂ = 0` —
**22.17**. Для сравнения: опубликованный сетевой С/Ш GW150914 равен **24**
(arXiv:1602.03837) — то есть модель воспроизводит сетевую величину на уровне
~4 %, тогда как прежние однометекторные числа сравнивались с однометекторными
19.5/13.3.

**Контроль мощности — то, без чего граница ничего не стоит.** Тот же фит
**восстанавливает** инжектированную задержку 20 мс (20.3 ± 5.1 мс) и **не
разрешает** заявленные статьёй 1.2 мс: разброс **4.6 мс** при самой амплитуде
1.2 мс. Значит «не ограничено на записанной сетке» означает «ниже порога этого
прибора», а не «фит слеп».

### 18.4 Мои ошибки этого хода — шесть, все пойманы контролем

| # | ошибка | как проявилась | кем поймана |
|---|---|---|---|
| 1 | не учтён световой ход между детекторами | когерентный С/Ш 12.0 против некогерентного 21.1 — когерентный **ниже**, чего быть не может при верном соотношении фаз | NC9 + прямое измерение знака |
| 2 | знак светового хода угадан, а не измерен | первый вариант дал 14.3; правильный `exp(+2πifτ_d)` даёт **21.08** при некогерентном 21.10 | измерение по сетке `(ι, ψ)` |
| 3 | NC1 сравнивал IMR-шаблон с **записанным числом LO-шаблона** | «падение» контроля по причине, не имеющей отношения к делу | сверка с `h3_direct_fit_imr.json` по шаблонам |
| 4 | NC8 «не восстановлено, если среднее отличается от 1.2 мс больше чем на 1 мс» — **ложный критерий** | среднее может случайно лечь на 1.2 мс при разбросе 4.6 мс, что есть **противоположность** разрешению | заменён на честный: разброс больше сигнала |
| 5 | ручная реализация HEALPix (две попытки) | пиксель 0 выходил на `θ = 0.09°` вместо `89.93°` — сверено с healpy | отказ от ручной реализации в пользу двух библиотек |
| 6 | верификатор: спутаны конвенции `rfft(x)` и `rfft(x)/fs` | все проверки краснели с фактором `fs` | V6 (записанные числа) |

**Пятая — самая полезная.** Переписывать HEALPix руками — ровно тот класс
ошибки, который этот проект ловит ход за ходом; поэтому ручная реализация
**удалена**, а не отлажена, и заменена двумя независимыми библиотеками плюс
физическим тестом на правильность чтения.

**Четвёртая записана, а не удалена:** ложный критерий остаётся в коде с
пометкой, как и предыдущие четыре неверных критерия проекта.

### 18.5 Независимая проверка (`verify_h3_joint.py`)

Не импортирует аудит. Свой PSD (свой Уэлч), своя GMST (другая алгебраическая
форма), своя алгебра пучка (тензор `D` выписан покомпонентно), **явная матрица**
`exp(2πift)` вместо обратного ДПФ, своя квадратура.

| проверка | результат |
|---|---|
| V1 отклики антенн | 5.4·10⁻⁸ |
| V2 отношение оптимальных С/Ш | **rel = 3.9·10⁻¹¹** |
| V3 сетевой С/Ш при `τ₂ = 0` | **rel = 7.3·10⁻¹¹** |
| V4 профиль по `τ₂`, оба ядра | peak 1.4·10⁻⁵, drop 6.5·10⁻⁴ |
| V5 круговой прогон | 7.5·10⁻⁶ |
| V6 записанные числа | **rel = 0.0** |
| **селфтест: порча артефакта** | **краснеет** |

**Отдельно измерено, а не спрятано:** оценщик PSD — часть **определения**
величины, не деталь реализации. При **среднем** Уэлче вместо медианного
однометекторные С/Ш уезжают на десятки процентов, потому что релиз O1 содержит
негауссовы глитчи (уже измерено в ходе 219). Независимость даёт **квадратура**
и **код**, а не смена оценщика; поэтому верификатор использует тот же оценщик, а
чувствительность к нему названа числом.

### 18.6 Что осталось открытым

1. **`τ(t)` по-прежнему не выведен из теории** — статья его не задаёт. Главный
   открытый пункт всего цикла; без него доказано вырождение *в принятой линейной
   модели*, а не невозможность проверить всю «now»-теорию.
2. **Спины выровненные**; прецессия `χ_p < 0.71` не покрыта.
3. **Задержка только фазовая** — амплитуда не задерживается.
4. **`η` зафиксирована** на опубликованном 29/36; совместный фит по
   `(M_c, η, спины, ι)` расширил бы профиль. Утверждение сделано в безопасную
   сторону (меньше свободных параметров).
5. **Остаток фактора ~3** в зазоре 37.3↔7.27 — назван числом, без IMR-формы не
   разделяется.
6. **L1 не достигает 13.3** даже при `ι = 90°` (минимум 13.58).

Артефакты: `artifacts/h3_joint_fit.json` (`art_32983e9ae39a`, sha256
`1daa8407…`), `artifacts/sky_samples.npz` (`art_d8f822f55257`, sha256
`0afea4a2…`). Оба воспроизводятся **байт-в-байт** при повторном прогоне.
Список утверждений — `work/CLAIMS.md`.


### 17.5 Что осталось открытым, и названо

1. **L1 не достигает 13.3** даже при ι = 90°: минимум 13.58. Отношение
   опубликованных С/Ш (1.466) больше отношения оптимальных С/Ш шаблона (1.102)
   — разница приходится на отклик антенн, которого в модели нет.
2. **τ(t) по-прежнему не выведен из теории** — статья его не задаёт. Это
   остаётся главным открытым пунктом всего цикла.
3. **Спины выровненные**: IMRPhenomT не покрывает опубликованную границу
   прецессии χ_p < 0.71. Названное упрощение.
4. **Задержка — только фазовая**: амплитуда не задерживается. Названный выбор.
5. **Совместный фит** по (M_c, η, спины, наклонение) расширил бы профиль;
   утверждение сделано в безопасную сторону (меньше свободных параметров).

---

# Ход 227 — ошибки, ложный критерий, и что было удалено, а не отлажено

## Шесть ошибок, все пойманы контролем

1. **Не учтён световой ход между детекторами.** Когерентный сетевой С/Ш выходил
   **12.0** против некогерентного **21.1** — когерентный оказался **ниже**
   некогерентного, чего при верном соотношении фаз быть не может (Коши — Шварц
   даёт равенство, когда фазы совпадают). Поймано **NC9**.
2. **Знак светового хода угадан, а не измерен.** Первый вариант
   (`exp(−2πifτ_d)`) дал 14.3; правильный `exp(+2πifτ_d)` даёт **21.08** при
   некогерентном 21.10. Знак определён **измерением по сетке `(ι, ψ)`**, а не
   аргументом. Ошибка была в том, что я сначала «вывел» знак рассуждением и
   получил его неверным.
3. **NC1 сравнивал IMR-шаблон с записанным числом LO-шаблона.** Записанное
   «7.273249726912924» относится к **ведущему инспиралу**, а «16.74195690868695»
   — к **полному IMR**; это разные шаблоны в одном артефакте. Контроль падал по
   причине, не имеющей отношения к делу. Исправлено: каждый шаблон сверяется с
   тем числом, которое он произвёл.
4. **Ложный критерий NC8 (первый черновик).** «Инжекция 1.2 мс не восстановлена,
   если восстановленное среднее отличается от 1.2 мс больше чем на 1 мс».
   **Ложно:** среднее может случайно лечь на 1.2 мс при разбросе 4.6 мс, а это
   **противоположность** разрешению. Заменён на честный: разброс больше сигнала
   (`std > 1.2 мс`). Это **пятый** неверный критерий в проекте; он записан, а не
   удалён.
5. **Ручная реализация HEALPix — две попытки, обе неверны.** Пиксель 0 выходил
   на `θ = 0.09°` вместо `89.93°` (сверено с healpy). Это ровно тот класс
   ошибки, который проект ловит ход за ходом. Решение: ручная реализация
   **удалена**, а не отлажена; взяты две независимые библиотеки
   (`healpy`, `astropy_healpix`), а правильность **чтения** FITS-таблицы выбрана
   физическим тестом (концентрация на кольце 6.9 мс против перемешанной карты).
6. **Верификатор: спутаны конвенции `rfft(x)` и `rfft(x)/fs`.** Все проверки
   краснели с фактором `fs` (V6 показывал 203956 против 16.74). Поймано V6 —
   сверкой с **записанными** числами.

## Ложный критерий, записанный, а не удалённый

**NC8 первого черновика** (см. пункт 4 выше): критерий «не восстановлено»
формулировался через **среднее**, тогда как правильная величина — **разброс**.
Критерий оставлен в истории этого раздела как пятый неверный критерий проекта
(после четырёх, записанных в ходах 219–226).

## Что было удалено, а не отлажено

**Ручная реализация NESTED HEALPix** (`healpix_nested.py` остался как
вспомогательный модуль, но в фите не используется). Две неверные попытки — это
сигнал: переписывать стандартную библиотеку руками дороже, чем взять две
независимые и проверить **результат** физикой. Правильность чтения карты теперь
держится не на доверии библиотеке, а на том, что принятое чтение кладёт **94 %**
вероятности на опубликованное кольцо 6.9 мс при **25×** превышении над
перемешанной картой.

## Числа, которых не было в проекте до этого хода

- Положение на небе GW150914: **RA = 134.80°, Dec = −69.79°**, 90 % площадь
  **616.4 deg²** (опубликовано 610), 50 % — 149.1 deg² (опубликовано 150).
- `τ_H1 − τ_L1 = −6.898 мс` из геометрии LAL и положения на небе.
- Отношение оптимальных С/Ш H1/L1 = **1.4749** против опубликованного
  однометекторного **1.4662** (совпадение 0.6 %).
- Когерентный сетевой С/Ш: пик **22.93** (ядро pn) / **23.07** (ядро imr), при
  `τ₂ = 0` — **22.17**; опубликованный **сетевой** 24.
- Односторонняя граница на `τ₂` на сетке ±38.4 мс: падение **3.496σ** (pn) и
  **2.352σ** (imr).
- `φ_c` **не идентифицируем** в когерентном фите с фиксированным положением —
  структурный факт, а не численный.

---

# Ход 228 — правки реестра по разбору владельца

Владелец (msg_00228) нашёл шесть расхождений между `CLAIMS.md`/`REPORT.md` и
фактами. Ни одно не было ошибкой физики — все были ошибками **формулировки**.
Первое он нашёл по **легенде реестра**, а не по числам.

## Что исправлено

1. **A4–A6 стояли под меткой [ДОКАЗАНО], не имея Lean-файла.** Введена
   четвёртая категория **[ВЫВЕДЕНО символьно]** (sympy + негативный контроль,
   формального доказательства нет). §F переписан: «в Lean доказаны A1–A3».
2. **B1 — оценка сверху, не измеренная граница.** Свободные `η` и спины
   поглотят ещё часть квадратичного члена (при фиксированных `η`, нулевых
   спинах уже 88 % поглощается). В допущения добавлено: положение на небе
   зафиксировано в максимуме правдоподобия (область 616 deg²).
3. **E4 переписано.** «Безопасная сторона» верна для «границы нет», но не для
   «граница есть». Для B1 фиксация `η`/спинов работает в опасную сторону.
4. **«+3.6 мс» определено.** Это правая точка пересечения `peak − 1σ` при
   внутреннем пике. На ±9.6 мс пик на левом краю → правая точка +5.4 мс; на
   ±38.4 мс пик внутренний (−18.0 мс) → +3.6 мс.
5. **B7: 0.6 % — согласие, не точность.** У отношения 19.5/13.3 допуск ~7 %.
   Настоящая проверка положения — задержка, площадь, контроль с перемешиванием.
6. **§1: разрешение фита ≈ 5 мс** (разброс 4.6 мс при 1.2 мс, 5.1 мс при
   20 мс); 69 мс помечено как результат **ведущего порядка**.
7. **Противоречие §9 ↔ E7 снято.** Оставлена одна формулировка: остаток
   фактора ~3 **закрыт** полным IMR.

## Две мои ошибки этого хода — обе в верификаторе, не в отчёте

**Ошибка 1 (V3).** Первый критерий искал в `lean/*.lean` **слово** «delay»/«tau»
и покраснел — потому что `Identifiability.lean` *посвящён* аддитивной задержке
в абстракции, и слово там есть, тогда как **содержания** A4–A6 (коэффициенты
`τ₀/τ₁/τ₂`, чирповая фаза, теорема Эйлера) нет. Критерий переписан на контент;
добавлен положительный контроль V3e: слово «delay» **обязано** присутствовать —
иначе тест не различал бы контент и слово.

**Ошибка 2 (V6).** Регулярка `0\.(?:8849|8716)` не нашла H1, потому что
записанное значение — `0.88487345…`, то есть цифры `88487`, а не `8849`.
Числа читаются из JSON **по ключу**, а не по образцу.

Обе — тот же класс, что и раньше: «сравниваю не то, что думаю». Пойманы тем,
что критерий покраснел на **верных** данных, и это заставило проверить сам
критерий, а не данные.

## Проверка

`verify_claims_228.py` (не импортирует ни один аудит): V1–V6 пересчитывают все
шесть правок из замороженных артефактов; три негативных контроля. `--selftest`
портит копию `h3_joint_fit.json` (пик уводится на край) и требует, чтобы
верификатор покраснел — **SELFTEST_PASS** (V1a и V1c краснеют).

---

# Ход 230 — красивый PDF, фигуры по реальной записи, и два дефекта в селфтесте

Вход владельца (msg_00230): сделать PDF красивым с хорошим шрифтом, сгенерировать
волновые формы, относящиеся к работе, и исправить формулировку авторства — соавтор
нашёл **шесть** ошибок, остальное агент нашёл сам.

## Что сделано

* **Типографика.** `preamble.tex` + `build_paper.sh`: pandoc + xelatex, основной
  шрифт **P052** (клон Palatino, есть кириллица), `microtype`, цветные линейки
  разделов, рамка абстракта, колонтитулы. Было — DejaVu Serif (экранный шрифт).
* **Пять фигур** (`make_figures.py` → `figures/*.pdf`) по **реальной** записи
  GW150914 и замороженным артефактам: (1) H1/L1 вокруг слияния в полосе 35–350 Гц,
  (2) измеренная ASD против значений артефакта, (3) отбелённая запись против
  полного `IMRPhenomT`, (4) ρ_opt против наклонения, (5) совместный профиль по τ₂.
  Скрипт печатает пересчитанные числа, чтобы их можно было сверить с артефактами
  (например `rho_opt_IMR = 31.659512701612417` — совпадает с записанным).
* **Авторство посчитано, а не вспомнено.** `count_errors.py` считает ошибки в
  `NOTES.md` и печатает использованное правило: **47** записанных ошибок, из них
  соавтор нашёл **11** (шесть в финальной вычитке реестра — ход 228, одна раньше —
  артефакт дискретизации в H4, ход 217, четыре — четыре самопротиворечия в
  собранном PDF, ход 230), остальные **36** — контроли, верификаторы и селфтесты
  агента. Дефект вёрстки (формула за полем) в число 47 **не входит** — см. пометку
  в разделе «Ход 230, третья часть» ниже. Раньше в бумаге стояло «catching a number
  of errors» — число без числа.
* **Новые проверки:** `verify_pdfs.py` (77 проверок: шрифты встроены, нужный шрифт
  использован, отсутствующие в P052 глифы **присутствуют в текстовом слое**,
  фигуры на месте и их подписи доходят до PDF) и `selftest_pdfs.py` (7/7 порч
  краснеют). В `verify_publication.py` добавлена проверка `MANIFEST.md` против
  байтов на диске.

## Два дефекта в селфтесте, найденные при этой работе

**Дефект 1 — C4 был вакуумным.** Проверка искала строку `autonomous research
agent` где угодно в документе, но эта строка встречается в `paper_en.md`
**дважды** (в YAML-списке авторов и в абзаце раскрытия). Удаление абзаца
оставляло проверку зелёной. Измерено прямым экспериментом: правка одного
вхождения → `PUBLICATION_CONSISTENT`. То есть контроль, который не мог упасть, —
**тот же класс, что и всё, что владелец ловил в этом проекте**. Исправлено:
проверка требует **всё предложение**, а не фразу.

**Дефект 2 — восстановление после C3 было сломано.** `edit(p, "stated", "not
derived", 1)` заменял **первое** вхождение слова «stated», а первое вхождение —
это `unstated` в абстракте. Файл оставался испорченным (`unnot derived`), и
следующий контроль C4 краснел **по чужой причине**. Измерено: после пары
порча/восстановление C3 файл не совпадал с исходным. Исправлено: восстановление
привязано к целому предложению, и после **каждой** пары добавлена проверка
байтового равенства (`restore after Cn is byte-exact`) — иначе следующий контроль
ничего не значит.

Оба дефекта — в **инструменте проверки**, а не в отчёте, и оба нашлись только
потому, что я стал править текст, который эти проверки читают.

## Что осталось открытым

Не изменилось: **τ(t) из теории не выведен** — статья его не задаёт. Это
по-прежнему главный открытый пункт цикла.

## Ход 230, третья часть: переполнение рамки, найденное владельцем

> **Пометка: это дефект ВЁРСТКИ, а не исследования.** Он записан здесь для
> полноты, но в число **47** ошибок (см. `count_errors.py`) **не входит** и ни
> одной из сторон не засчитывается: формула, вышедшая за поле страницы, — это
> сборка PDF, а не утверждение о физике. Поэтому разделение в §0 статьи —
> **11 / 36** (соавтор / агент), в **трёх** вычитках, а не 12 / 35 в четырёх.

Владелец (op_f3e47e6add64) открыл уже собранный PDF и увидел, что формула

    (ΔV_total)^{1/3}/R_s(62) = 2.0083

выходит за правую рамку страницы, и попросил проверить **все** рамки, а важные
формулы при необходимости вывести отдельно по центру.

**Что оказалось.** Переполнений было **два класса**, и оба подтвердились
измерением, а не на глаз:

1. **Формула внутри пункта списка.** Она была одним неразрывным токеном длиной
   ~113 pt в строке, где оставалось ~40 pt. Исправлено: формула вынесена на
   отдельную строку (отступ в 6 пробелов в Markdown), как и просил владелец.
2. **URL в списке литературы** (`https://gwosc.org/eventapi/html/event/GW150914/v3`)
   — один неразрывный токен длиннее строки. Исправлено двумя средствами:
   `xurl` (разрыв URL по любому символу) и угловые скобки вокруг адреса.

**Отдельно — то, что пришлось выяснить, а не предположить.** Первый вариант
проверки читал боксы слов через `pdftotext -bbox` и показывал «переполнение» на
знаке минус в `10⁻¹⁴`. Я вырезал эту строку из растрированной страницы и
посмотрел: **чернила внутри рамки**, а бокс шире — `pdftotext` отдаёт бокс
**продвижения** глифа, а не его чернил. То есть проверка по боксам давала
**ложное** переполнение того же размера, что и настоящее. Поэтому проверка
переписана на **чернила**: страница растеризуется при 200 dpi, и ищутся тёмные
пиксели вне рамки. Геометрия берётся из самого преамбула (a4paper, поля 2,4 см),
а не из догадки.

**И одна типографская правка, которую эта проверка потребовала.** `microtype`
по умолчанию включает **выступ** (protrusion): глиф законно выходит на ~2 pt в
поле. Это стандартная типографика и глазу не видно, но делает вопрос «все ли
слова внутри рамки?» неразрешимым без допуска, а допуск того же размера скрыл бы
и настоящее переполнение. Выступ **выключен** (`-V microtypeoptions=protrusion=false`),
растяжение и кернинг оставлены. Теперь проверка **точна**, а не с допуском.

**Проверка стала инструментом и порчей.** `verify_frames.py` — 4 проверки
(обе статьи, все страницы). В селфтест добавлена **C10b**: рамка расширяется до
0,4 см, PDF пересобирается, и проверка обязана покраснеть — **краснеет**
(`FRAME_BAD`). То есть доказано, что она **может** упасть, а не просто зелёная.

Итог: **12 страниц EN, 13 RU, ни одного пикселя чернил вне рамки** (худшее
превышение 0 px при допуске 1 px). Селфтест: **29/29**.

## Ход 230, вторая часть: четыре самопротиворечия, найденные владельцем

Владелец (op_e3a156582467) прочитал уже собранный PDF и нашёл четыре места, где
работа противоречит **своей же** главной оговорке. Все четыре подтвердились
чтением файлов, и все четыре исправлены.

1. **Подпись Fig. 5 против абзаца над ней.** Подпись называла профиль по τ₂
   «первым ограничением, которое дают сами данные», а текст рядом прямо говорил
   «это оценка силы ограничения сверху, а не измеренная граница». Подпись
   переписана: теперь она повторяет оговорку, а не отменяет её.
2. **«ι ≈ 59.2° измерено».** Это не измерение, а **подгонка**: наклонение, при
   котором оптимальное С/Ш шаблона сравнялось бы с опубликованным 19,5. Слово
   «измеренная» убрано, добавлено «подогнанное наклонение», и это вынесено
   отдельным пунктом 9 в §9 (в EN-версии пункт 8 был **дважды продублирован** —
   один и тот же текст про конвенцию π-границы; дубликат убран).
3. **Порядок величины.** Было «на три порядка выше» для 69 мс против 0,61 мс.
   Измерено: 68,78 / 0,6109 = **112,6**, то есть **два** порядка. Исправлено на
   «два», и добавлено, с чем именно идёт сравнение (с собственным масштабом
   механизма R_s/c, а не с заявленными 1,2 мс, которые в 57 раз меньше 69 мс).
   Рядом добавлен абзац **«Какое С/Ш какое»**: 0,087 — аналитическая
   откалиброванная PSD; 0,174 / 0,127 — измеренная PSD H1/L1; 0,13–0,17 —
   диапазон этой пары; 22,93 / 23,07 — **сетевые** величины из §6, сравнимы
   только с сетевым 24. Различие — в PSD, а не в физике.
4. **Дубликат в §9.** В EN-версии пункт 8 стоял дважды с одним и тем же текстом.
   Убран.

**Все четыре стали машинными проверками** (`verify_publication.py`, блок
self-contradiction checks) и **все четыре — порчами в селфтесте** (C11–C15), так
что теперь доказано, что каждая проверка **может** упасть. Селфтест: **28/28**.

Это тот же класс, что и всё, что владелец ловил в этом проекте: **утверждение
сильнее содержания**, причём на этот раз — внутри одного абзаца, а не между
разделами. И заметьте: все четыре нашёл читатель по формулировкам, ни одну — по
числам. Числа были верны; неверными были слова вокруг них.


---

# Ход 231 — два дефекта в собранном PDF, найденные владельцем

Вход владельца (msg_00231): прочитать собранный PDF и проверить два места.

## Дефект 1 — §2 противоречил §8

§2 (Метод) говорил: «в заголовок допускается только утверждение из регистра 1».
§8 говорил: заголовок держится на **формальных** A1–A3 **и** на **символьных**
A4–A6. A4–A6 — регистр 2, не 1. То есть §2 запрещал то, на что §8 опирался.

В `work/REPORT.md` этого противоречия **нет**: там §1 сформулирован мягче — «ни
одно утверждение категории гипотеза в заголовке не используется». То есть ошибка
появилась **при сборке PDF**: при переводе в статью ограничение было усилено с
«гипотеза исключена» до «только регистр 1», и §8 не был перепроверен.

**Исправлено.** §2 (EN и RU) теперь повторяет §8 и реестр (`CLAIMS.md` §F):
заголовок держится на **регистрах 1 и 2**; исключена категория **гипотеза**.

## Дефект 2 — «Two» вместо трёх

§3.2: «Two further formal results…», а под ним **три** пункта (`RelativityOfNow`,
`PastHypothesis`, `HomogeneousUniverse`). Числительное не совпадало со списком.

**Исправлено.** «Three further» / «Три дальнейших».

## Как это проверено машинно, а не памятью

Оба дефекта — класса «заголовок/подпись сильнее содержания», который владелец ловит
в проекте ход за ходом. Поэтому каждый превращён в **проверку** в
`verify_publication.py`, а каждая проверка — в **порчу** в `selftest_publication.py`:

* **C16** — §2 снова ограничивает заголовок регистром 1 → проверка краснеет.
* **C17** — числительное в §3.2 снова не совпадает с числом пунктов → краснеет.

Проверка C17 считает пункты в блоке §3.2 и требует, чтобы числительное им
соответствовало (EN: one/two/three…; RU: один/два/три…). Проверка C16 требует,
чтобы §2 называл исключённую категорию и регистры 1 и 2, и чтобы старой
формулировки «только регистр 1» не было.

## Правка §0 по приоритету владельца (op_3df3ebdb1f95)

Владелец: убрать «a formula running past the page margin» из списка ошибок — это
вёрстка, не работа; пересчитать разделение на **11/36, три ревью**; согласовать
§0, `NOTES.md` и `count_errors.py`.

**Что оказалось при пересчёте.** Дефект вёрстки **не входил** в число 47 (его нет
ни в одной из COUNTED-секций `count_errors.py`), но **входил** в список OWNER — то
есть вычитался из доли агента. Это и давало неверные 12/35. Убрав его из OWNER,
получаю ровно **11/36 в трёх вычитках** — то, что просил владелец.

**Согласовано в пяти местах** (все правки сделаны, все проверены машиной):
`paper_en.md` §0, `paper_ru.md` §0, `README.md`, `.zenodo.json` (описание депозита)
и `NOTES.md` (раздел «Ход 230, третья часть» помечен как вёрсточный и исключён из
47; устаревшее «45/7/38» в разделе «Что сделано» исправлено на 47/11/36).

**Машинная проверка согласия (C18, C19).** `verify_publication.py` теперь
**запускает `count_errors.py` как подпроцесс** и требует, чтобы числа в §0 статьи
(EN и RU) и в описании `.zenodo.json` **равнялись** его выводу. C18 портит §0
(11→12), C19 портит `.zenodo.json` — обе краснеют. Число в бумаге больше не может
разойтись с записью молча.

## Проверка

* `verify_publication.py` — **204 ok, 0 failed** → `PUBLICATION_CONSISTENT`.
* `verify_pdfs.py` — 77 ok; `verify_frames.py` — 4 ok (ни одного пикселя вне рамки).
* `selftest_publication.py` — **36/36**; `selftest_pdfs.py` — **7/7**.
* `count_errors.py` — `COUNT_CONSISTENT` (47 = 11 + 36).
* Оба PDF пересобраны; оба дефекта исправлены **в самом PDF** (проверено
  `pdftotext`): §2 «rests on registers 1 and 2 only», §3.2 «Three further».
* Артефакты измерений **не тронуты**: все 16 совпадают с `ARTIFACTS.md` побайтово.
* Регрессия: `mutation_control.py` exit 1 (предсуществующее, 15/20) и
  `sky_prep.py` exit 1 в системном интерпретаторе — **не регрессия**: `healpy`
  стоит в локальном venv, и `env/venv/bin/python sky_prep.py` даёт **exit 0**.

## Что осталось открытым

Не изменилось: **τ(t) из теории не выведен** — статья его не задаёт. Главный
открытый пункт цикла.

# Ход 232 — владелец: убрать из статьи фразу про «четыре самопротиворечия»

**Замечание владельца (msg_00232).** «прочитав собранный PDF и найдя, что работа
противоречит собственной главной оговорке в четырёх местах. убери это пожалуйста
потому что это ничего не меняет, проста не нужная информация для читателя и после
того как увидит проста закроет pdf.»

**Что сделано.** Фраза убрана из **всех четырёх читательских документов**, где она
стояла: `paper_en.md` §0, `paper_ru.md` §0, `README.md`, `.zenodo.json` (описание
депозита). Формулировка заменена на нейтральную «и **четыре** — при вычитке
собранного PDF» / «and **four** in a review of the built PDF» — **число 11/36
сохранено**, убрана только фраза, рекламирующая самопротиворечие.

**Что НЕ тронуто и почему.** Внутренняя запись — `NOTES.md` и `count_errors.py` —
сохраняет все четыре ошибки как таковые: это честный журнал, и владелец просил
убрать фразу из **статьи**, а не из записи. Раздел «Ход 230, вторая часть: четыре
самопротиворечия» и строка в `count_errors.py` остаются; `count_errors.py` даёт
прежние **11 / 36 / 47**.

**Машинная проверка (C20).** Новое замечание владельца превращено в проверку, а не
только в правку текста (standing intention). `verify_publication.py` теперь требует,
чтобы ни в одном читательском документе (EN, RU, README, `.zenodo.json`) не
встречались строки `contradicting its own`, `противоречит собственной`,
`in four places`, `в четырёх местах`. `selftest_publication.py` C20 возвращает
фразу в EN-статью — проверка обязана покраснеть: **краснеет**.

**Проверка.** `verify_publication.py` — **208 ok, 0 failed** → `PUBLICATION_CONSISTENT`;
`selftest_publication.py` — **38/38**; `selftest_pdfs.py` — 7/7; `verify_pdfs.py` —
77 ok; `verify_frames.py` — 4 ok; `count_errors.py` — `COUNT_CONSISTENT` (47 = 11 + 36).
Оба PDF пересобраны (`art_1edc98486b30`, `art_f43db7864982`); фраза отсутствует в
собранном PDF (проверено `pdftotext`, не по исходникам). Артефакты измерений не
тронуты.

# Ход 233 — владелец: не писать об ошибках, найденных при сборке PDF

**Замечание владельца (msg_00233).** «четыре—при вычитке собранного PDF. ВООБЩЕ НЕ
НУЖЕН ПИСАТ ОШИБКИ КОТОРИЙ СДЕЛАНО и найдено ВО ВРЕМЯ КОГДА СОБРАЛ pdf. это будет
смешно даже когда увитят.»

**Это второй заход.** Ход 232 я понял уже́, чем он был: владелец просил убрать
**фразу**, рекламирующую самопротиворечие, а я убрал только часть — само слово
«вычитке собранного PDF» оставил, заменив «прочитав собранный PDF и найдя…» на
нейтральное «при вычитке собранного PDF». То есть я подчинил букву и не расслышал
смысл: **читателю не нужен ни факт вычитки, ни то, сколько ревью было и на каком
этапе** они случились. Владелец повторил, и он прав дважды.

## Что сделано

**Из четырёх читательских документов** (`paper_en.md` §0, `paper_ru.md` §0,
`README.md`, `.zenodo.json`) убраны:
* «при вычитке собранного PDF» / «in a review of the built PDF»;
* «в трёх вычитках» / «in three reviews» — число ревью тоже не нужно читателю;
* абзац про вёрстку («дефект вёрстки PDF… в число 47 не входит») — целиком;
* «четыре — при вычитке собранного PDF» → «**четыре** дальнейшие ошибки
  формулировок»; источник четырёх больше не называется.

**Число сохранено: 11 / 36 / 47.** Убраны обстоятельства находки, не счёт.
Соавтор по-прежнему назван с одиннадцатью ошибками — это часть утверждения об
авторстве, и она осталась.

**Что НЕ тронуто.** Внутренняя запись: `NOTES.md` (разделы «Ход 230, вторая и
третья часть», «Ход 231») и `count_errors.py` сохраняют и этап, и обстоятельства.
Владелец просил убрать это из **подачи**, а не из журнала; стереть этап из журнала
значило бы потерять, **где** ошибка была поймана, — а это и есть самое полезное в
ней (ошибка сборки показывает, что проверки статьи не делились на текст и вёрстку).

## Машинная проверка (C21 + блок READER)

Замечание превращено в проверку, а не только в правку — это standing intention, и
на этот раз он сработал: **то же замечание пришло дважды**, значит проверки не было.

* `verify_publication.py` — блок `_PDF_STAGE_CLAUSE`: ни в одном читательском
  документе не должно быть `review of the built PDF`, `in three reviews`,
  `вычитке собранного PDF`, `в трёх вычитках`, `typesetting`, `дефект вёрстки`,
  `not among the 47`. Плюс три негативные проверки: §0 EN, §0 RU и описание
  депозита **не должны** аргументировать исключение вёрстки.
* `selftest_publication.py` **C21** возвращает абзац про вёрстку в §0 — проверка
  краснеет. **C20** перенастроена на новую формулировку «четыре дальнейшие
  ошибки формулировок» и тоже краснеет.

## Регрессия, которую я внёс сам, и как она поймана

Укоротив §0, я **сдвинул разбиение на страницы**: статья EN стала 11 страниц вместо
12. При отложенном размещении плавающих объектов LaTeX сгрудил две фигуры на одну
страницу, не смог их уложить и выдал `Overfull \vbox (289.93pt too high) ... while
\output is active` — подпись Figure 2 ушла **за нижний край листа**, то есть чернила
вне блока.

**Как измерено, а не угадано.** Собрал старый текст заново в отдельной копии:
старая сборка — 12 страниц, `Overfull \vbox` = 0, чернил вне блока нет; новая — 11
страниц, два `Overfull \vbox`, чернила на физическом краю страниц 5 и 6. То есть
регрессия моя, а не предсуществующая.

**Исправлено в преамбуле:** `\floatplacement{figure}{H}` — фигуры ставятся **точно
там, где написаны текстом**. Проверено на матрице вариантов (placeins, запрет
нижних флоатов, щедрые float-страницы, уменьшение фигур): ни один не убрал
переполнение; `[H]` убрал — 0 `Overfull`, 12 страниц EN и 13 RU, ни одного пикселя
вне блока. Цена названа: фигура может оставить пустое место внизу страницы.

## Серия сборки, в правильном порядке

`make_figures.py` → `build_paper.sh` → `make_manifest.py` → `verify_publication.py`
(так записано в докстринге `make_manifest.py`, потому что манифест покрывает и
сами верификаторы).

## Проверка

| проверка | итог |
|---|---|
| `verify_publication.py` | **212 ok, 0 failed** → `PUBLICATION_CONSISTENT` |
| `selftest_publication.py` | **40/40** (C20 перенастроена, C21 добавлена) |
| `verify_pdfs.py` / `selftest_pdfs.py` | 77 ok / **7/7** |
| `verify_frames.py` | 4 ok, 0 failed → `FRAME_OK` |
| `count_errors.py` | `COUNT_CONSISTENT` (47 = 11 + 36) |

Фраза отсутствует **в собранных PDF** (проверено `pdftotext`, 0 совпадений по
`proofread`, `built PDF`, `four places`, `PDF assembly`, `not among the 47`,
`вёрстк`, `вычитк`, `собранн`). Артефакты измерений не тронуты.

**Ещё один дефект в своих же проверках, найденный этой работой.** Порча **P3** в
`selftest_pdfs.py` искала строку «The co-author found **eleven** of them, in
three» — формулировку, которой после правки §0 больше нет. Порча стала тихо
**вакуумной**: `baseline green` оставался зелёным, а P3 краснела по чужой причине.
Это ровно тот класс, который владелец ловит в проекте ход за ходом, и он снова
вылез из правки текста, читаемого проверками. Исправлено; проверено, что P3
краснеет **чистым** FAIL («states the co-author found 11 errors: no 'found eleven'
phrase»), а не падением.

## Что осталось открытым

Не изменилось: **τ(t) из теории не выведен** — статья его не задаёт.
