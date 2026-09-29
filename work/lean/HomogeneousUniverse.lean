import Mathlib

/-!
# Homogeneous universe: "time creation" is not identifiable from the expansion

The "now" theory's own §6 concedes that in a homogeneous universe its effect is
indistinguishable from dark energy. This file turns that concession into a
theorem about *what the data can identify*, in the same additive form as
`Identifiability.lean`.

Model.  The observable is the expansion history: the measured Hubble rate as a
function of scale factor, `obs H δ a = H a + δ a`, where `H` is the baseline
expansion the theory would predict and `δ` is the extra contribution attributed
to created time. Only `obs` is measured.

Two statements, and they are a dichotomy:

* `creation_not_identifiable_unpinned` — if the baseline `H` is free (as it is
  when dark energy is fitted, which is exactly what cosmology does), then `δ` is
  NOT identifiable: every `δ` is reproduced by a shifted baseline with a
  different `δ`. So no measurement of the expansion can single out a
  created-time contribution.
* `creation_identifiable_pinned` — if the baseline `H` is independently fixed,
  `δ` IS identifiable. This is the anti-vacuum side: the degeneracy is a
  statement about which quantity is free, not a statement that nothing can ever
  be measured.

Anti-vacuity witnesses are included: a concrete pair showing the unpinned
degeneracy is real (`demo_unpinned_pair`), and the corresponding pinned instance
(`demo_pinned_identifiable`). Without them the first theorem could be satisfied
by a one-element parameter space.
-/

namespace TimeFlowCosmo

/-- The observable: measured expansion rate at scale factor `a`, given a baseline
`H` and a created-time contribution `δ`. -/
def obs (H δ : ℝ → ℝ) (a : ℝ) : ℝ := H a + δ a

/-- `δ` is identifiable from the observable when the baseline is FIXED: distinct
`δ` give distinct observables. -/
def IdentifiablePinned (H : ℝ → ℝ) (δ : ℝ → ℝ) : Prop :=
  ∀ δ', (∀ a, obs H δ a = obs H δ' a) → δ = δ'

/-- `δ` is identifiable from the observable when the baseline is FREE: no other
pair `(H', δ')` produces the same observable unless `δ' = δ`. -/
def IdentifiableFree (H : ℝ → ℝ) (δ : ℝ → ℝ) : Prop :=
  ∀ δ' H', (∀ a, obs H' δ' a = obs H δ a) → δ = δ'

/-- **Non-identifiability (free baseline).** If the baseline expansion is fitted
freely — as it is in any cosmological parameter estimation that fits dark energy
— then the created-time contribution `δ` cannot be recovered from the expansion
history. Every `δ` is mimicked exactly by a shifted baseline with a different
`δ'`. -/
theorem creation_not_identifiable_unpinned (H δ : ℝ → ℝ) : ¬ IdentifiableFree H δ := by
  intro hident
  -- the baseline H + δ with NO creation reproduces the same observable
  have h0 : δ = fun _ => 0 :=
    hident (fun _ => 0) (fun a => H a + δ a) (by intro a; simp [obs])
  -- the baseline H + δ - 1 with creation ≡ 1 also reproduces it
  have h1 : δ = fun _ => 1 :=
    hident (fun _ => 1) (fun a => H a + δ a - 1) (by intro a; simp only [obs]; ring)
  have e0 : δ 0 = 0 := congrFun h0 0
  have e1 : δ 0 = 1 := congrFun h1 0
  linarith

/-- **Identifiability (pinned baseline).** If the baseline expansion `H` is fixed
independently of the fit, then `δ` is identifiable: the observable determines it
exactly. This is the anti-vacuum branch — the degeneracy above is about which
quantity is free, not about an impossibility of measuring anything. -/
theorem creation_identifiable_pinned (H : ℝ → ℝ) (δ : ℝ → ℝ) :
    IdentifiablePinned H δ := by
  intro δ' h
  funext a
  have ha := h a
  simp only [obs] at ha
  linarith

/-- **The degeneracy is real, not vacuous.** A concrete created-time
contribution (constant 1) and a genuinely DIFFERENT parameter (constant 0) with a
different baseline produce the SAME observable. So
`creation_not_identifiable_unpinned` is not talking about an empty space of
alternatives. -/
theorem demo_unpinned_pair :
    ∃ (δ δ' H' : ℝ → ℝ), δ ≠ δ' ∧ (∀ a, obs H' δ' a = obs (fun _ => 1) δ a) := by
  refine ⟨fun _ => 1, fun _ => 0, fun _ => 2, ?_, ?_⟩
  · intro h
    have := congrFun h 0
    norm_num at this
  · intro a
    simp only [obs]
    norm_num

/-- **The pinned branch is inhabited.** With the baseline fixed at the constant
`1`, distinct `δ` give distinct observables, so `IdentifiablePinned` is not a
vacuous statement about a one-element parameter space. -/
theorem demo_pinned_identifiable :
    IdentifiablePinned (fun _ => 1) (fun _ => 0) :=
  creation_identifiable_pinned (fun _ => 1) (fun _ => 0)

/-- **The pinning is exactly the line that separates the two branches.** With the
baseline pinned to `H`, two observables that agree force `δ' = δ`; the
hypothesis that makes the unpinned theorem fail is precisely "some other baseline
exists", which is what a dark-energy fit supplies. Stated as: if the baseline is
pinned to `H`, then agreement of observables forces agreement of `δ`. -/
theorem baseline_unique_when_pinned (H δ H' δ' : ℝ → ℝ)
    (h : ∀ a, obs H' δ' a = obs H δ a) (hH : ∀ a, H' a = H a) :
    ∀ a, δ' a = δ a := by
  intro a
  have ha := h a
  simp only [obs] at ha
  rw [hH a] at ha
  linarith

end TimeFlowCosmo
