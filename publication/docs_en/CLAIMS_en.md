# Claim ledger — time_flow_research cycle (revision of round 228)

Owner requirement (msg_00227): **a separate list of claims tagged
“proven / derived / measured under such assumptions / hypothesis”**, so that the
headline is not stronger than the content. Here every claim of the report is
assigned to one of **four** categories (the fourth was added in round 228 — see
below), and for each one the thing that holds it is named.

Notation (four categories, not three — revision of round 228):

* **[PROVEN]** — a **formal** result: Lean 4, no `sorry`, no custom axioms (the
  axiom audit yields only `propext, Classical.choice, Quot.sound`). Independent of
  data, noise, waveform. The tag is applied **only** where there is a Lean file
  and a theorem name.
* **[DERIVED symbolically]** — an analytic result in closed form, verified with
  symbolic algebra (sympy) and a negative control. Independent of data, **but
  there is no formal proof in Lean**. The category was introduced because A4–A6
  stood under the tag “proven” without having a Lean file.
* **[MEASURED]** — a numerical result on the real GW150914 data under named
  assumptions. The assumptions are listed in the row; removing them changes the
  number.
* **[HYPOTHESIS]** — a claim the project did not test or could not test.

---

## A. The core of the result

| # | claim | status | what holds it |
|---|---|---|---|
| A1 | In the additive model `M(θ,τ) = Aθ + Bτ`, `τ` is identifiable ⟺ `B` is injective and `range A ⊓ range B = ⊥` | **[PROVEN]** | `lean/Identifiability.lean`, `identifiable_iff` |
| A2 | If `range B ≤ range A`, then `τ` is not identifiable | **[PROVEN]** | `not_identifiable_of_range_le` |
| A3 | The criterion is not vacuous: both cases exist in one formal field | **[PROVEN]** | `lean/Instances.lean`, `demo_identifiable` / `demo_not_identifiable` |
| A4 | A constant delay `τ₀` is exactly equivalent to a shift of `t_c`; a linear delay `τ̇` to a rescaling of `M_c` | **[DERIVED symbolically]** for the leading phase: `ΔΨ = 2x^{3/2}τ₀/M_c ∝ ∂Ψ/∂t_c` and `−(5/128)τ₁x^{−5/2} ∝ ∂Ψ/∂M_c` (sympy, self-check S6 in `analysis.py`); **[MEASURED]** numerically: absorbed `1 − 1.4·10⁻¹⁵` and `1 − 1.7·10⁻¹⁵`. **No Lean file** | assumptions: linear model, leading order, analytic PSD calibrated to SNR 20 |
| A5 | The first non-degenerate term is quadratic, suppressed as `(v/c)⁸` | **[DERIVED symbolically]** via the exponent `x^{−13/2}` versus `x^{−5/2}` (sympy); **[MEASURED]** SNR 0.087 at 1.2 ms. **No Lean file** | assumptions: leading order |
| A6 | The degeneracy is NOT an artifact of the leading order: the identities `∂Ψ/∂t_c = 2πf` and `f·∂Ψ_chirp/∂f = M_c·∂Ψ_chirp/∂M_c` hold at **every** PN order | **[DERIVED symbolically]** (Euler’s theorem in `v`; no PN coefficient enters), verified with sympy at PN orders 0…3.5 with a negative control (`pn_orders.py`). **No Lean file** | — |
| A7 | Relativity-of-simultaneity objection: there is a boost that reverses the sign of `Δt'` when `|Δt| < |Δx|` | **[PROVEN]** | `lean/RelativityOfNow.lean`, `exists_boost_reversing_time_order` |
| A8 | With a free base expansion (dark-energy fit) the parameter `δ` is not identifiable; with it fixed, it is | **[PROVEN]** | `lean/HomogeneousUniverse.lean`; 3 negative controls go red |
| A9 | A non-zero arrow of time requires a boundary that breaks T-symmetry (on a family closed under reversal the arrow sums to zero) | **[PROVEN]** | `lean/PastHypothesis.lean`, `arrow_sum_zero`; non-vacuity: `paramArrow_ne_zero`, witness `demo_closure_is_the_line` (0 versus 2) |

## B. H3 on the real GW150914 data

| # | claim | status | what holds it |
|---|---|---|---|
| B1 | The joint H1+L1 fit with **one** `τ₂` and free `(M_c, t_c, φ_c)` gives, on the ±38.4 ms grid, a **one-sided bound** on `τ₂`: the profile drops by 3.50σ (pn kernel) and 2.35σ (imr kernel). **This is an upper estimate of the strength of the constraint, not a measured bound.** Free `η` and spins would absorb more of the quadratic term (§4b: with `η` fixed and zero spins, 88 % of the term is already absorbed by the source parameters), so the true constraint is **weaker** | **[MEASURED]** as an upper estimate | `h3_joint_fit.py`, `art_32983e9ae39a`; assumptions: spins 0, phase-only delay, `η = 29/36`, **sky position fixed at the maximum-likelihood point of the published map** (90 % area 616 deg² — the position is NOT fitted jointly with `τ₂`) |
| B2 | On the recorded ±9.6 ms grid the profile is **not constrained** (drop 1.38σ / 0.85σ) | **[MEASURED]** | same |
| B1b | What the number “+3.6 ms” means: it is the **right crossing of the level `peak − 1σ`** (the first grid point to the right of the peak where the profile falls 1.0 below the peak), not the grid edge and not a confidence interval. On the ±9.6 ms grid the peak sits at the **left edge** (`−9.6 ms`, where the profile is still rising), so the right crossing is **+5.4 ms**, not +3.6; on ±38.4 ms the peak is interior (`−18.0 ms`) and the right crossing is **+3.6 ms**. There is no left crossing on either grid — the profile does not drop by 1σ on the left. So the “bound of +3.6 ms” is a right crossing at an **interior** peak; it lands inside ±9.6 ms only because on the narrow grid the peak runs to the edge and the crossing moves to +5.4 ms | **[MEASURED]** (definition of the number) | `h3_joint_fit.py`, field `one_sigma_high_amplitude_ms`; verified by recomputing the profile from the artifact |
| B3 | The same fit **recovers** a 20 ms injection (20.3 ± 5.1 ms) and **does not resolve** the paper’s claimed 1.2 ms (spread 4.6 ms at an amplitude of 1.2 ms) | **[MEASURED]** | NC7/NC8; without this, “not constrained” is indistinguishable from “the fit is blind” |
| B4 | `φ_c` is **not identifiable** in a coherent network fit with a fixed sky position: it multiplies every detector’s coefficient by the same phase | **[DERIVED symbolically]** structurally (from the form `C_d = F₊A − iF×B`; no Lean file), **[MEASURED]** as the identity `hc₀ = −i·hp₀` (deviation 3.3·10⁻⁴) | — |
| B5 | The whole detector response transformation collapses to **one complex coefficient** `C_d = F₊A(ι) − iF×B(ι)` | **[MEASURED]** | phenomxpy convention: `hp = A e^{2iφ}hp₀`, `hc = B e^{2iφ}hc₀` |
| B6 | The optimal SNR at `ι = 0` is 25.31 (H1) / 17.16 (L1); their ratio is 1.4749 | **[MEASURED]** | assumptions: spins 0, leading mode only, sky position from the map |
| B7 | The ratio of optimal SNRs 1.4749 against the published single-detector 19.5/13.3 = 1.4662 — **agreement at the 0.6 % level, but this is NOT a 0.6 %-precision test.** The published 19.5 and 13.3 are observed quantities with noise of order ±1, so the ratio itself carries a tolerance of about **7 %** (dominated by the 13.3 ± 1 contribution; quadrature ≈ 9 %). The 0.6 % agreement lies deep inside that tolerance and therefore does not constrain the position more strongly than ±7 %. The real checks on the position are the three quantities named by the owner: the **delay** `τ_H1 − τ_L1 = −6.898 ms` against the published 6.9 (+0.5/−0.4) ms, the **area** 616.4 deg² against 610 deg², and the **shuffle control** (the accepted map reading gives 25× over a shuffled map) | **[MEASURED]**, precision named honestly | NC5; `sky_prep.py` |
| B8 | The coherent network SNR at `τ₂ = 0` is 22.17, at the peak 22.93 (pn kernel) / 23.07 (imr kernel) | **[MEASURED]** | `verify_h3_joint.py`, V3/V4, independent path |
| B9 | The coherent SNR ≤ the incoherent SNR for every template (Cauchy–Schwarz) | **[MEASURED]** | NC9 |
| B10 | The sign of the light-travel delay was determined **by measurement**: `exp(+2πifτ_d)` gives 21.08, the opposite sign 16.62 against an incoherent 21.10 | **[MEASURED]** | the sign is not postulated |
| B11 | `τ_H1 − τ_L1 = −6.898 ms` from the LAL geometry and the sky position | **[MEASURED]** | consistent with the published 6.9 (+0.5/−0.4) ms |
| B12 | The π-bound on the quadratic-delay amplitude is **1.667 ms** (τ₂ = 0.0023370 s⁻²), and the paper’s 1.2 ms uses **72 %** of that budget. **The bound is convention-dependent:** it is computed for the paper’s picture (delay zero at the band start, growing toward merger). Placing the parabola’s vertex at coalescence gives **24.96 ms** under the same π-criterion, fifteen times looser | **[MEASURED]** under the named convention; the convention dependence is named as a number | `h3_grid_bounds.py`, `art_2343189a11d1`; both values re-derived independently in `publication/verify_publication.py` |

## C. Sky position (new in this round)

| # | claim | status | what holds it |
|---|---|---|---|
| C1 | Maximum-likelihood position: RA = 134.80°, Dec = −69.79° | **[MEASURED]** | `sky_prep.py`; two independent HEALPix readers (healpy and astropy_healpix) agree to 1.8·10⁻¹⁵ rad |
| C2 | 90 % credible area = 616.4 deg², 50 % = 149.1 deg² | **[MEASURED]** | published 610 deg² (90 %) and 150 deg² (50 %) — arXiv:1602.03840 |
| C3 | The correct FITS-table reading was chosen **by measurement**, not by trusting a library: `data['PROB'].reshape(-1)` as NESTED places 94 % of the probability on the 6.9 ms annulus (concentration 0.164), whereas the `healpy.read_map` array does not (0.023) | **[MEASURED]** | control: a shuffled map gives 0.0061 on average, 0.0065 maximum (20 runs); the accepted reading is **25×** the null |
| C4 | `F₊² + F×²` does not depend on the polarization angle `ψ` | **[MEASURED]** | NC4, spread < 10⁻⁹ |
| C5 | The projection `(r_H1 − r_L1)·n/c` requires a GMST rotation: celestial coordinates of the map against the Earth system of the detector vertices. Without the rotation the peak sits at −3.4 ms with concentration 0.026; with it, at −6.9 ms with 0.164 | **[MEASURED]** | — |

## D. H4 (the paper’s arithmetic)

| # | claim | status | what holds it |
|---|---|---|---|
| D1 | The paper’s two “independent” numbers (0.6 and 1.2 ms) differ by exactly the geometric factor 2.0083 | **[MEASURED]** | `h4_paper_arithmetic.py` v2 |
| D2 | The volume-weighted mean of `(1+z)` over the stated shell `[1.5, 4] Rₛ` is 1.25034 | **[MEASURED]** | assumption: weight “time ∝ volume”; weighting by created volume `V(M)−V(0)` gives 1.28596, flat gives 1.24173 |
| D3 | The agreement of 1.250 with the paper’s text **depends on the choice of weight**; the paper gives no physical derivation of the weight | **[MEASURED]** | `h4_weight_audit.py` |
| D4 | The claim “the multiplier 1.25 cannot be produced” is **refuted** (my own error of round 217) | **[MEASURED]** | root 2.7778 inside `[1.5, 4]`; v2 scans the whole interval |

## E. Still open / hypothesis

| # | claim | status | why |
|---|---|---|---|
| E1 | The specific delay function `τ(t)` is **not derived** from the “now”-theory | **[HYPOTHESIS]** | the paper does not specify it; therefore what is proven is degeneracy *within the adopted linear model*, not the impossibility of testing the whole theory |
| E2 | Aligned spins only (`IMRPhenomT`); precession `χ_p < 0.71` is not covered | **[HYPOTHESIS]** | a named simplification |
| E3 | The delay is phase-only; the amplitude is not delayed | **[HYPOTHESIS]** | a named model choice |
| E4 | `η` is fixed at the published 29/36; a joint fit over `(M_c, η, spins, ι)` would widen the profile | **[HYPOTHESIS]** | **Not a “safe side”.** For the claim “no bound exists” fewer free parameters is safe, but B1 claims the opposite — that a bound **does** exist. For that claim, fixing `η` and spins works in the **unsafe** direction: free `η` and spins would absorb more of the quadratic term and weaken the bound. That is why B1 is tagged as an upper estimate |
| E5 | The sign of the ET/CE benefit from a wide band is undetermined without `τ(t)` | **[HYPOTHESIS]** | with a fixed amplitude it falls (0.087→0.0047), with a fixed `τ₂` it rises (0.087→7.625) |
| E6 | L1 does not reach the published 13.3 even at `ι = 90°` (minimum 13.58) | **[MEASURED]**, but the residue is not explained | the residue is attributed to antenna response, which the model lacks |
| E7 | The factor-~3 residue in the 37.3↔7.27 gap is **CLOSED** by the real IMR template (`h3_imr_check.py`, round 224): more than half of the deficit was a leading-order defect; the residue is explained by orientation (`ι ≈ 59.2°`, round 226) | **[MEASURED]** | the earlier phrasing “not separable without the IMR form” is **obsolete** and was withdrawn (revision of round 228). What remains open is not this residue but the L1 residue (E6) |

---

## F. What the headline is allowed to claim

The report’s headline (“the prediction is not falsifiable on GW150914 for a
structural reason”) rests on **A1–A3** — **formal** results in Lean, independent
of data — and on **A4–A6** — analytic identities derived symbolically (sympy)
with negative controls, but **not** formalized in Lean. The phrasing I am obliged
to use: **in Lean, A1–A3 are proven; A4–A6 are derived symbolically.** The earlier
phrase “the headline rests on the formal results A1–A6” was stronger than the
content — exactly the error the ledger was introduced to prevent.

Everything concerning numerical bounds on `τ₂` (**B1–B3**) is **[MEASURED]** under
the listed assumptions and **will change** if the spins, `η` or the phase-only
nature of the delay are relaxed. **B1 is an upper estimate of the strength of the
constraint**, not a measured bound.

No claim in the **[HYPOTHESIS]** category is used in the headline.

---

## G. Revisions of round 228 (per the owner’s review, msg_00228)

Five remarks, all accepted and incorporated. None was an error of physics — all
five were errors of **formulation**, and the first was found by the owner from the
ledger’s legend, not from the numbers.

| # | remark | what was done |
|---|---|---|
| 1 | A4–A6 are tagged **[PROVEN]** but have no Lean file — there it is analytics and sympy | A **fourth tag [DERIVED symbolically]** was introduced; A4–A6 were moved into it; §F was rewritten: “in Lean, **A1–A3** are proven”. Verified by `verify_claims_228.py` V3: `lean/*.lean` contains neither the coefficients `τ₀/τ₁/τ₂`, nor the chirp phase, nor Euler’s theorem |
| 2 | B1 is argued too confidently; E4 calls fixing `η`/spins the “safe side”, but for the claim “a bound exists” that is not safe | B1 was tagged as an **upper estimate of the strength of the constraint**; the assumption that the **sky position is fixed at the maximum-likelihood point** (area 616 deg²) was added to B1; E4 was rewritten: for “no bound” fixing is safe, for “a bound exists” it is unsafe |
| 3 | “+3.6 ms” appears only on the ±38.4 ms grid, though 3.6 lies inside ±9.6 ms, where “none” is stated | A definition was introduced: **+3.6 ms is the right crossing of the level `peak − 1σ`** at an **interior** peak. On ±9.6 ms the peak sits at the **left edge**, and the right crossing is **+5.4 ms**; on ±38.4 ms the peak is interior (`−18.0 ms`) and the crossing is **+3.6 ms**. Verified by `verify_claims_228.py` V1/V2 |
| 4 | “0.6 %” in B7 reads more precisely than it is: the ratio 19.5/13.3 carries a ~7 % tolerance | B7 was rewritten: 0.6 % is agreement, but **not** a 0.6 %-precision test; the real checks on the position are the delay, the area and the shuffle control |
| 5 | §1 places “an amplitude of 69 ms is needed” next to “a 20 ms injection is recovered (20.3 ± 5.1)” | The current fit resolution **≈ 5 ms** was stated (spread 4.6 ms at 1.2 ms and 5.1 ms at 20 ms); 69 ms was tagged as a **leading-order** result |
| 6 | §9 says “the factor-~3 residue is closed with the real IMR”, while E7 says “not separable without the IMR form” | The contradiction was removed: one phrasing kept — the residue is **closed** by the full IMR (`h3_imr_check.py`); the obsolete “not separable” was deleted from E7 |

**What this says about the ledger.** The ledger I wrote myself contained item 1 —
and the owner found it from the **legend**, not from the numbers. The tag
“proven” stood on claims that had no proof, and that is exactly the error
“headline stronger than content” for which the ledger was introduced. The lesson:
the category of a claim must be checked **by a machine** (is there a Lean file
with that theorem name?), not by the author’s memory. `verify_claims_228.py` V3 is
the first such machine check of a tag.
