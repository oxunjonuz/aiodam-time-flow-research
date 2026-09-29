# Pre-registration — cycle 2026-09-27

**This file is frozen before the decisive run.** Everything below was written
BEFORE the experiment, the formalization and the analysis. Edits after the run are
allowed only as a separate dated amendment at the end of the file.

Task folder: `/work/time_flow_research_20260927/`
Working folder for results: `/work/time_flow_research_20260927/work/`

> English translation of `work/PREREGISTRATION.md`. The Russian original is kept
> alongside; the English text is the primary one for publication.

---

## 0. What was already done BEFORE the freeze (disclosure)

During reconnaissance (before this file) I:

1. Read `README.md`, the Müller–Maguire paper (arXiv:1606.07975v1) in full, and
   the video transcript.
2. Checked the sha256 of the three local sources — they match those stated in the
   README.
3. **Reproduced equation (4.3) of the paper** — the volume estimate ΔV for
   GW150914. Result: `(ΔV_total)^(1/3)/c = 1.227 ms`, ratio
   `(ΔV)^(1/3)/Rs(62) = 2.008`. This is a reproduction of someone else’s
   calculation, not my test.
4. Established that `gwosc.org` and `dcc.ligo.org` are unreachable from the
   container (connection timeout), i.e. **the raw LIGO data could not be
   obtained**. This is a blocker, not a result; it remains a blocker.
5. Checked the tools: Lean 4.19.0 + Mathlib (prebuilt oleans,
   `/work/Shopify/audit-work/mathlib`), z3 4.13.3, numpy/scipy/sympy. No `h5py`.

None of items 3–5 is a test of the hypotheses below.

---

## 1. Research vector

Of the three directions proposed by the owner, I choose the **logical-mathematical
audit of falsifiability** of the “now”-theory prediction, with a numerical check
on a gravitational-wave signal model.

Reason for the choice: the README itself poses the question “is the prediction
specified precisely enough to be distinguished from GR and from a variation of the
source parameters”. That is a question of **identifiability**, not of detector
sensitivity. It can be answered without the raw LIGO data — so the blocker of
item 4 does not kill it. And it is formally checkable.

What I do NOT do, and why:
- I do not test the GW150914 data (no data access; see item 4).
- I do not build a new theory of gravity.
- I do not derive consequences for consciousness and free will from gravitational
  waves.

---

## 2. The model everything is computed on

The phase term of the merger frequency sweep, leading order:

```
Ψ(f) = 2π f t_c − φ_c − π/4 + (3/128)(π M_c f)^(−5/3)
```

where `M_c` is the chirp mass, `t_c` the coalescence time, `φ_c` the coalescence
phase.

The “now”-theory prediction: the radiation is delayed by `τ(t)`, where `t` is the
emission time. In the frequency domain this gives an additive phase correction

```
ΔΨ(f) = 2π f τ(t(f)),      t(f) = t_c − (5/256) M_c^(−5/3) (π f)^(−8/3)
```

Key point: **all parameters (M_c, t_c, φ_c) are fitted freely in the LIGO
analysis.** The identifiability question is whether the vector `ΔΨ` lies in the
linear span of the vectors `∂Ψ/∂M_c`, `∂Ψ/∂t_c`, `∂Ψ/∂φ_c`.

---

## 3. Hypotheses (frozen)

**H1 — exact degeneracy of a constant delay.**
A constant delay `τ(t) ≡ τ_0` gives `ΔΨ = 2π f τ_0`, which is **identically** a
shift `t_c → t_c + τ_0`. Hence `τ_0` is not identifiable at any signal-to-noise
ratio and for any noise curve.
*Falsifier:* exhibit a functional of `h(f)` sensitive to `τ_0` and insensitive to
`t_c`. Prediction: no such functional exists.

**H2 — exact degeneracy of a linearly growing delay.**
If `τ(t) = τ_0 + τ̇(t − t_ref)`, the `∝ τ̇` term gives
`ΔΨ = −(5/3) τ̇ · Ψ_chirp + (a term ∝ f)`, i.e. **identically** equivalent to a
rescaling of the chirp mass
`M_c → M_c (1 − (5/3)τ̇)^(3/5) ≈ M_c (1 − τ̇)`.
Hence `τ̇` is not identifiable separately from `M_c`.
*Falsifier:* show that the residual after the best fit over `M_c` with an injected
`τ̇ ≠ 0` differs from zero by more than machine precision.

**H3 — only the second and higher orders are identifiable.**
For a smooth `τ(t)` the phase contribution expands in powers of `(π M_c f)^(2/3)`;
the zeroth- and first-order terms in `t` are degenerate (H1, H2), and the residue
is suppressed relative to the leading phase by a factor of order `(v/c)²`.
*Falsifier:* show that the quadratic term gives a residual SNR ≥ 1 on GW150914
with aLIGO noise and the effect size claimed by the paper.

**H4 — the paper’s quantitative estimate does not have the stated precision.**
The paper gives 1–1.5 ms. I check: (a) the “dimensional analysis” (0.6 ms) and the
“Schwarzschild computation” (1.2 ms) are not independent estimates but differ by
exactly the geometric factor `(ΔV)^(1/3)/Rs`; (b) the multiplier 1.5/1.2 = 1.25 is
introduced by the words “with the additional gravitational redshift” without
derivation.
*Falsifier:* find a derivation of the multiplier 1.25 in the paper, or show that
`(ΔV)^(1/3)/Rs(62) = 1`.

**H5 — the relativity-of-simultaneity objection to the “growing block”.**
If “existing now” is a spacelike hypersurface (the edge of growth), then in 1+1
Minkowski geometry, for any two events with `|Δt| < |Δx|/c` there exist inertial
frames in which the sign of `Δt'` is opposite. Hence the predicate “exists” cannot
be simultaneously (a) invariant and (b) equal to “lies in the past of the edge”.
*Falsifier:* exhibit an invariant predicate on events that coincides with “in the
past of the edge” in one frame and is frame-independent.

---

## 4. Thresholds (assigned BEFORE the run)

| # | Quantity | Threshold | Meaning |
|---|---|---|---|
| T1 | max\|ΔΨ residual\| after the fit, constant delay | ≤ 1e-10 rad | H1 confirmed |
| T2 | max\|ΔΨ residual\| after the fit, linear delay | ≤ 1e-10 rad | H2 confirmed |
| T3 | \|ρ(t_c, τ_0)\| in the Fisher matrix | ≥ 0.999999 | exact degeneracy |
| T4 | \|ρ(M_c, τ̇)\| | ≥ 0.999999 | same for M_c |
| T5 | SNR of the quadratic residual at the paper’s effect | < 1 | H3 confirmed |
| T6 | SNR of the same residual at ×100 the effect | ≥ 5 | **power control** |
| T7 | SNR of the residual at zero delay | ≤ 1e-9 | zero control |
| T8 | False-alarm fraction on pure noise (100 realizations) | ≤ 5 % | calibration |

Threshold T6 is the most important: without it, “did not detect” is
indistinguishable from “cannot detect”.

---

## 5. Negative controls (each must go red)

1. **K1 — zero injection.** τ = 0 → residual exactly 0.
2. **K2 — revert H1.** Fix `t_c` → the residual from a constant delay is non-zero.
3. **K3 — revert H2.** Fix `M_c` → the residual from a linear delay is non-zero.
4. **K4 — scrambled phase.** A random smooth addition → the residual is large.
5. **K5 — ×100 amplification.** The quadratic effect ×100 → detected (T6).
6. **K6 — sign flip.** τ̇ → −τ̇ → the residual is equally zero.
7. **K7 — a foreign noise curve.** A different S_n(f) → the correlations remain 1.

---

## 6. Formalization (Lean 4 + Mathlib)

`work/lean/Identifiability.lean`:

- **T_ident.** For `M(θ,τ) = A θ + B τ`: `τ` is identifiable
  (⟺ `M(θ,τ) = M(θ',τ') ⟹ τ = τ'`) **if and only if** `B` is injective and
  `range A ⊓ range B = ⊥`.
- **C1.** If `range B ≤ range A`, then `τ` is not identifiable.
- **C2.** If `B` is injective and `range A ⊓ range B = ⊥`, then `τ` is
  identifiable.

`work/lean/RelativityOfNow.lean`:

- **T_sim.** In 1+1 Minkowski: if `|Δt| < |Δx|/c`, there is a boost (|v| < c)
  that reverses the sign of `Δt' = γ(Δt − vΔx/c²)`.

Both theorems — without `sorry` and without `axiom`. The presence of `sorry` is a
failure of the formalization, not a partial success.

---

## 7. Independent check

`work/verify_independent.py` — **imports** no module of the main analysis, reads
only the artifacts from disk, recomputes the key numbers with its own code. A
discrepancy > 1e-9 is a failure.

`work/verify_lean.sh` — rebuilds both `.lean` files from scratch + greps for
`sorry`/`axiom`.

---

## 8. What will count as a result

- **H1/H2 confirmed** → the prediction is not falsifiable on GW150914 for a
  structural reason (degeneracy), independently of detector sensitivity. This
  translates the paper’s own caveat §5 from prose into a theorem.
- **H1/H2 refuted** → a measurable combination exists; I compute its SNR.
- **Blocker** → if Lean/Mathlib does not build or the procedure does not converge,
  I record it as a blocker with the exact command and error text.

---

## 9. Boundaries declared in advance

1. The waveform is leading order, not full IMRPHENOM/NR. The H1/H2 conclusion does
   not depend on this (it is algebraic); the H3 conclusion does.
2. The noise curve is an analytic approximation, not a measured PSD.
3. The GW150914 data are not used; the parameters are taken from arXiv:1602.03837.
4. No claims about consciousness, free will, or the “reality of the present”.
5. The video transcript is not a source of physical claims.
