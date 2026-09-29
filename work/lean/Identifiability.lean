import Mathlib

/-!
# Identifiability of an additive delay in a linear model

Model: observations `M(θ,τ) = A θ + B τ`, where `θ` are nuisance parameters
(fitted freely) and `τ` is the parameter of interest (here: the "now"-theory
delay). This file proves the exact criterion for when `τ` can be recovered from
the data at all — independently of noise, signal-to-noise ratio or the detector.

Physics application: `A` are the source-parameter directions (`M_c`, `t_c`,
`φ_c`) and `B` is the delay direction. `τ` is identifiable iff the delay
direction is not a linear combination of the source-parameter directions.
-/

namespace TimeFlow

variable {Θ Τ H : Type*}
variable [AddCommGroup Θ] [Module ℝ Θ]
variable [AddCommGroup Τ] [Module ℝ Τ]
variable [AddCommGroup H] [Module ℝ H]

/-- `τ` is identifiable in the additive model `A θ + B τ`: two parameter pairs
that produce the same observation must have the same `τ`. -/
def Identifiable (A : Θ →ₗ[ℝ] H) (B : Τ →ₗ[ℝ] H) : Prop :=
  ∀ θ θ' τ τ', A θ + B τ = A θ' + B τ' → τ = τ'

/-- **Main criterion.** In the additive model `A θ + B τ`, the parameter `τ` is
identifiable exactly when `B` is injective and the ranges of `A` and `B` meet
only at zero. -/
theorem identifiable_iff (A : Θ →ₗ[ℝ] H) (B : Τ →ₗ[ℝ] H) :
    Identifiable A B ↔
      Function.Injective B ∧ Disjoint (LinearMap.range A) (LinearMap.range B) := by
  constructor
  · intro hident
    refine ⟨?_, ?_⟩
    · intro τ τ' hτ
      exact hident 0 0 τ τ' (by simpa using hτ)
    · rw [disjoint_iff_inf_le]
      rintro x ⟨⟨θ, hθ⟩, ⟨τ, hτ⟩⟩
      have heq : A θ + B 0 = A 0 + B τ := by
        rw [map_zero, map_zero, add_zero, zero_add]
        exact hθ.trans hτ.symm
      have hτ0 : τ = 0 := (hident θ 0 0 τ heq).symm
      rw [Submodule.mem_bot, ← hτ, hτ0, map_zero]
  · rintro ⟨hinj, hdisj⟩ θ θ' τ τ' heq
    have h3 : B (τ - τ') = A (θ' - θ) := by
      rw [map_sub, map_sub]
      exact (sub_eq_sub_iff_add_eq_add).mpr (by rw [add_comm (B τ), heq])
    have hmemA : B (τ - τ') ∈ LinearMap.range A := ⟨θ' - θ, h3.symm⟩
    have hmemB : B (τ - τ') ∈ LinearMap.range B := ⟨τ - τ', rfl⟩
    have hzero : B (τ - τ') = 0 := hdisj.le_bot ⟨hmemA, hmemB⟩
    have hsub : τ - τ' = 0 := hinj (by rw [map_zero]; exact hzero)
    exact sub_eq_zero.mp hsub

/-- **Collapse (degeneracy).** If the delay direction `B` lies inside the range
of the nuisance directions `A`, then `τ` is *not* identifiable: every effect of
`τ` can be absorbed by a shift of `θ`. -/
theorem not_identifiable_of_range_le {A : Θ →ₗ[ℝ] H} {B : Τ →ₗ[ℝ] H}
    (hle : LinearMap.range B ≤ LinearMap.range A) (hB : B ≠ 0) :
    ¬ Identifiable A B := by
  intro hident
  rw [identifiable_iff] at hident
  obtain ⟨hinj, hdisj⟩ := hident
  obtain ⟨τ, hτ⟩ : ∃ τ, B τ ≠ 0 := by
    by_contra hcon
    push_neg at hcon
    exact hB (LinearMap.ext hcon)
  have hmem : B τ ∈ LinearMap.range A := hle ⟨τ, rfl⟩
  have hzero : B τ = 0 := hdisj.le_bot ⟨hmem, ⟨τ, rfl⟩⟩
  exact hτ hzero

/-- **Zero direction.** A delay direction with no effect (`B = 0`) is not
identifiable whenever the delay parameter space is nontrivial. -/
theorem not_identifiable_of_zero [Nontrivial Τ] (A : Θ →ₗ[ℝ] H) :
    ¬ Identifiable A (0 : Τ →ₗ[ℝ] H) := by
  intro hident
  rw [identifiable_iff] at hident
  obtain ⟨hinj, _⟩ := hident
  obtain ⟨τ, hτ⟩ := exists_ne (0 : Τ)
  exact hτ (hinj (by simp))

end TimeFlow