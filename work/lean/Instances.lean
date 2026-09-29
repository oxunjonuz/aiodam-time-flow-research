import Mathlib
import Identifiability

/-!
# Anti-vacuity instances for the identifiability criterion

A theorem about identifiability is worthless if no model ever satisfies it.
This file exhibits, in the *same* formal setting:

* `demo_identifiable` — a concrete model where the delay IS identifiable
  (the criterion's positive branch is inhabited);
* `demo_not_identifiable` — a concrete model where the delay is NOT
  identifiable, matching the physics case (delay direction inside the range of
  the nuisance directions).

Together they show the criterion `identifiable_iff` is a real dichotomy, not a
statement that is trivially true or trivially false.
-/

namespace TimeFlow

/-- Concrete model: observation space `ℝ × ℝ`, nuisance parameter `ℝ` acting on
the first coordinate, delay parameter `ℝ` acting on the second coordinate.
The two directions are disjoint, so the delay is identifiable. -/
def A1 : ℝ →ₗ[ℝ] (ℝ × ℝ) where
  toFun a := (a, 0)
  map_add' a b := by ext <;> simp
  map_smul' c a := by ext <;> simp

def B1 : ℝ →ₗ[ℝ] (ℝ × ℝ) where
  toFun b := (0, b)
  map_add' a b := by ext <;> simp
  map_smul' c a := by ext <;> simp

theorem B1_injective : Function.Injective B1 := by
  intro a b h
  simpa [B1] using congrArg Prod.snd h

theorem ranges_disjoint : Disjoint (LinearMap.range A1) (LinearMap.range B1) := by
  rw [disjoint_iff_inf_le]
  rintro x ⟨⟨a, ha⟩, ⟨b, hb⟩⟩
  rw [Submodule.mem_bot]
  have h1 : x.1 = 0 := by rw [← hb]; rfl
  have h2 : x.2 = 0 := by rw [← ha]; rfl
  exact Prod.ext h1 h2

/-- Positive instance: the criterion's hypothesis holds, hence the delay is
identifiable. Without this, `identifiable_iff` could be vacuously one-sided. -/
theorem demo_identifiable : Identifiable A1 B1 :=
  identifiable_iff A1 B1 |>.mpr ⟨B1_injective, ranges_disjoint⟩

/-- Concrete model matching the physics case: the delay direction `B2` is the
*first* coordinate, which already lies in the range of the nuisance map `A2`.
This is the algebraic image of "a constant delay is just a shift of `t_c`". -/
def A2 : ℝ →ₗ[ℝ] (ℝ × ℝ) where
  toFun a := (a, 0)
  map_add' a b := by ext <;> simp
  map_smul' c a := by ext <;> simp

def B2 : ℝ →ₗ[ℝ] (ℝ × ℝ) where
  toFun b := (b, 0)
  map_add' a b := by ext <;> simp
  map_smul' c a := by ext <;> simp

theorem range_le : LinearMap.range B2 ≤ LinearMap.range A2 := by
  rintro x ⟨b, hb⟩
  exact ⟨b, by rw [← hb]; rfl⟩

theorem B2_ne_zero : B2 ≠ 0 := by
  intro h
  have := congrArg (fun (f : ℝ →ₗ[ℝ] (ℝ × ℝ)) => f 1) h
  simpa [B2] using this

/-- Negative instance: the delay direction is inside the range of the nuisance
map, so the delay is NOT identifiable. This is the formal content of the
"constant delay = shift of `t_c`" degeneracy. -/
theorem demo_not_identifiable : ¬ Identifiable A2 B2 :=
  not_identifiable_of_range_le range_le B2_ne_zero

end TimeFlow