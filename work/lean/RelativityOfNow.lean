import Mathlib

/-!
# Relativity of simultaneity against a "growing edge of now"

A growing-block picture says: exactly one spatial slice (the "edge of now") is
real, and the block grows. If that slice is a spacelike hypersurface, then the
predicate "is real now" must be a property of the event together with the slice,
and one would like it to be observer-independent.

This file proves the obstruction in the simplest case: in 1+1 Minkowski
spacetime, for any two events with a spacelike separation and distinct times,
there is an inertial frame in which their time order is reversed. Hence
"earlier than" is not an observer-independent relation between
spacelike-separated events.
-/

namespace TimeFlow

/-- A 1+1 Minkowski boost with velocity parameter `v` (units `c = 1`),
`|v| < 1`. Returns the transformed time and space coordinates. -/
noncomputable def boost (v t x : ℝ) : ℝ × ℝ :=
  let γ := (1 - v ^ 2) ^ (-(1 / 2 : ℝ))
  (γ * (t - v * x), γ * (x - v * t))

/-- **Relativity of simultaneity.** If two events are spacelike separated
(`|Δt| < |Δx|`) and not simultaneous (`Δt ≠ 0`), then there is a boost with
`|v| < 1` that reverses their time order. -/
theorem exists_boost_reversing_time_order {t x : ℝ} (ht : t ≠ 0) (h : |t| < |x|) :
    ∃ v : ℝ, |v| < 1 ∧ (t - v * x) * t < 0 := by
  have hx : x ≠ 0 := by
    intro hx0
    rw [hx0, abs_zero] at h
    exact not_lt_of_ge (abs_nonneg t) h
  set u : ℝ := t / x with hu
  have hu_ne : u ≠ 0 := by
    rw [hu]; exact div_ne_zero ht hx
  have hu_pos : 0 < |u| := abs_pos.mpr hu_ne
  have hu_lt : |u| < 1 := by
    rw [hu, abs_div]
    rw [div_lt_one (abs_pos.mpr hx)]
    exact h
  refine ⟨u * ((1 + 1 / |u|) / 2), ?_, ?_⟩
  · rw [abs_mul]
    have hk : |(1 + 1 / |u|) / 2| = (1 + 1 / |u|) / 2 :=
      abs_of_pos (by positivity)
    rw [hk]
    have hmul : |u| * ((1 + 1 / |u|) / 2) = (|u| + 1) / 2 := by
      field_simp
      ring
    rw [hmul]
    linarith
  · have hvx : (u * ((1 + 1 / |u|) / 2)) * x = t * ((1 + 1 / |u|) / 2) := by
      rw [hu]; field_simp; ring
    rw [hvx]
    have hk1 : 1 < (1 + 1 / |u|) / 2 := by
      have : 1 < 1 / |u| := (one_lt_div hu_pos).mpr hu_lt
      linarith
    nlinarith [sq_pos_of_ne_zero ht]

end TimeFlow