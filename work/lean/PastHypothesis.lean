import Mathlib

/-!
# The arrow of time comes from the boundary, not from the laws

The thermodynamic arrow of time is usually stated as: entropy increases toward
the future. But the microscopic laws (Newton, Schrödinger, Einstein) are
time-reversal symmetric — they do not distinguish past from future. So where
does the arrow come from? The standard answer (Boltzmann, and the "past
hypothesis" of Albert, 2000) is: from a **boundary condition** — the universe
started in a low-entropy state — not from the dynamics.

This file formalizes the *logical skeleton* of that answer, so the claim is a
theorem and not a slogan. The structure:

* `IsTraj T x` — `x : ℤ → α` is a bi-infinite trajectory of the reversible
  dynamics `T` (a bijection).
* `rev R x` — the time reversal of `x` under an involution `R`.
* `TRSymmetric T R` — the dynamics is time-reversal symmetric: `R T = T⁻¹ R`.
* `IsArrow R a` — an "arrow": a function on trajectories that flips sign under
  time reversal.

Theorems:

* `rev_involutive` — reversal is an involution on trajectories.
* `rev_isTraj` — **the reversal of a trajectory is a trajectory.** The set of
  histories of a time-symmetric dynamics is closed under time reversal.
* `arrow_vanishes_on_symmetric` — a trajectory that is its own reversal carries
  **no arrow**: any arrow must assign it zero.
* `rev_allowed_of_invariant` — a time-reversal-**invariant** boundary condition
  keeps the allowed set closed under reversal.
* `arrow_sum_zero` — **conclusion.** Over any finite reversal-closed family of
  trajectories, an arrow sums to exactly zero. So a time-symmetric boundary
  cannot make time "point" anywhere; a nonzero arrow requires a boundary that
  **breaks** the time-reversal symmetry.

Anti-vacuity (`Instances`-style): `demo_arrow_vanishes` exhibits a concrete
model in which the phenomenon occurs; `Bfuture_not_invariant` and
`demo_boundary_breaks_reversal` exhibit a boundary that does break the symmetry,
so the hypothesis of the conclusion is not vacuous.

**Scope, stated honestly.** This is the logic of the past-hypothesis argument,
not the physics of entropy increase. Nothing here says the universe *does* have
a low-entropy boundary; it says that *if* the laws are time-symmetric, then the
arrow cannot come from them, and must come from a symmetry-breaking boundary.
-/

namespace TimeFlow

variable {α : Type*}

/-- A bi-infinite trajectory of the dynamics `T`. -/
def IsTraj (T : α ≃ α) (x : ℤ → α) : Prop := ∀ n, x (n + 1) = T (x n)

/-- Time reversal of a trajectory under the map `R`. -/
def rev (R : α ≃ α) (x : ℤ → α) : ℤ → α := fun n => R (x (-n))

/-- The dynamics `T` is time-reversal symmetric under `R`: `R T = T⁻¹ R`. -/
def TRSymmetric (T R : α ≃ α) : Prop := ∀ y, R (T y) = T.symm (R y)

/-- Reversal is an involution whenever `R` is. -/
lemma rev_involutive (R : α ≃ α) (hR : ∀ y, R (R y) = y) (x : ℤ → α) :
    rev R (rev R x) = x := by
  funext n
  simp only [rev]
  rw [show -(-n) = n from by ring]
  exact hR (x n)

/-- Key commutation: `R` carries the inverse dynamics to the forward dynamics. -/
lemma rev_comm (T R : α ≃ α) (h : TRSymmetric T R) (y : α) :
    R (T.symm y) = T (R y) := by
  have h1 := h (T.symm y)
  rw [T.apply_symm_apply] at h1
  rw [h1, T.apply_symm_apply]

/-- **The reversal of a trajectory is a trajectory.** The histories of a
time-reversal-symmetric dynamics are closed under time reversal. -/
lemma rev_isTraj (T R : α ≃ α) (h : TRSymmetric T R) {x : ℤ → α}
    (hx : IsTraj T x) : IsTraj T (rev R x) := by
  intro n
  simp only [rev]
  have hstep : x (-n) = T (x (-(n + 1))) := by
    have := hx (-(n + 1))
    simpa [show -(n + 1) + 1 = -n from by ring] using this
  have hstep2 : x (-(n + 1)) = T.symm (x (-n)) := by
    rw [hstep, T.symm_apply_apply]
  rw [hstep2, rev_comm T R h]

/-- An "arrow": a function on trajectories that flips sign under time reversal. -/
def IsArrow (R : α ≃ α) (a : (ℤ → α) → ℤ) : Prop := ∀ x, a (rev R x) = - a x

/-- **A self-reverse trajectory carries no arrow.** -/
theorem arrow_vanishes_on_symmetric (R : α ≃ α) (a : (ℤ → α) → ℤ)
    (ha : IsArrow R a) {x : ℤ → α} (hx : rev R x = x) : a x = 0 := by
  have h := ha x
  rw [hx] at h
  omega

/-- Boundary condition: the state at time `0` lies in `B`. -/
def Allowed (B : Set α) (x : ℤ → α) : Prop := x 0 ∈ B

/-- **A time-reversal-invariant boundary keeps the allowed set closed under
reversal.** -/
theorem rev_allowed_of_invariant (R : α ≃ α) {B : Set α}
    (hB : ∀ y, y ∈ B ↔ R y ∈ B) {x : ℤ → α} (hx : Allowed B x) :
    Allowed B (rev R x) := by
  simp only [Allowed, rev]
  rw [show -(0 : ℤ) = 0 from by ring]
  exact (hB (x 0)).mp hx

/-- **Conclusion: the arrow averages to zero over a reversal-closed family.**
If a finite set of trajectories is closed under time reversal, any arrow sums to
exactly zero on it. A time-symmetric boundary therefore cannot orient time; a
nonzero arrow requires a boundary that breaks the symmetry. -/
theorem arrow_sum_zero (R : α ≃ α) (hR : ∀ x, rev R (rev R x) = x)
    (a : (ℤ → α) → ℤ) (ha : IsArrow R a)
    (S : Finset (ℤ → α)) (hS : ∀ x ∈ S, rev R x ∈ S) :
    ∑ x ∈ S, a x = 0 := by
  have h1 : ∑ x ∈ S, a (rev R x) = ∑ x ∈ S, a x := by
    refine Finset.sum_bij (fun x _ => rev R x) ?_ ?_ ?_ ?_
    · intro x hx; exact hS x hx
    · intro x _ y _ hxy
      have := congrArg (rev R) hxy
      rwa [hR x, hR y] at this
    · intro b hb
      exact ⟨rev R b, hS b hb, hR b⟩
    · intro x _; rfl
  have h2 : ∑ x ∈ S, a (rev R x) = ∑ x ∈ S, -a x :=
    Finset.sum_congr rfl (fun x _ => ha x)
  rw [h2, Finset.sum_neg_distrib] at h1
  linarith

/-! ## Anti-vacuity: a concrete model

Dynamics = successor on `ℤ` (a bijection), reversal = negation. The trajectory
`n ↦ n` is its own reversal, so it carries no arrow; and the boundary
`{n | 0 ≤ n}` is *not* reversal-invariant, so the hypothesis of the conclusion
is genuinely breakable. -/

/-- Successor, as a bijection of `ℤ`. -/
def succ : ℤ ≃ ℤ where
  toFun n := n + 1
  invFun n := n - 1
  left_inv n := by ring
  right_inv n := by ring

/-- Negation, as an involution of `ℤ`. -/
def negZ : ℤ ≃ ℤ where
  toFun n := -n
  invFun n := -n
  left_inv n := by ring
  right_inv n := by ring

lemma negZ_involutive : ∀ n : ℤ, negZ (negZ n) = n := by
  intro n; simp [negZ]

lemma succ_TRSymmetric : TRSymmetric succ negZ := by
  intro y
  show negZ (succ y) = succ.symm (negZ y)
  change -(y + 1) = succ.symm (-y)
  have hs : succ.symm (-y) = -y - 1 := rfl
  rw [hs]
  ring

lemma id_traj : IsTraj succ (fun n : ℤ => n) := by
  intro n; simp [succ]

lemma id_rev_fixed : rev negZ (fun n : ℤ => n) = (fun n : ℤ => n) := by
  funext n
  simp [rev, negZ]

/-- Anti-vacuity: in the concrete model any arrow vanishes on the symmetric
trajectory, so the phenomenon is not empty. -/
theorem demo_arrow_vanishes (a : (ℤ → ℤ) → ℤ) (ha : IsArrow negZ a) :
    a (fun n : ℤ => n) = 0 :=
  arrow_vanishes_on_symmetric negZ a ha id_rev_fixed

/-- The reversal of a trajectory is again a trajectory, in the concrete model. -/
theorem demo_rev_is_traj : IsTraj succ (rev negZ (fun n : ℤ => n)) :=
  rev_isTraj succ negZ succ_TRSymmetric id_traj

/-- Anti-vacuity for the conclusion: a reversal-closed family exists and the
sum of any arrow over it is zero. -/
theorem demo_sum_zero (a : (ℤ → ℤ) → ℤ) (ha : IsArrow negZ a) :
    ∑ x ∈ ({(fun n : ℤ => n)} : Finset (ℤ → ℤ)), a x = 0 := by
  refine arrow_sum_zero negZ (fun x => rev_involutive negZ negZ_involutive x)
    a ha _ ?_
  intro x hx
  rw [Finset.mem_singleton] at hx
  rw [hx, id_rev_fixed]
  exact Finset.mem_singleton_self _

/-- A boundary condition that is **not** time-reversal invariant. -/
def Bfuture : Set ℤ := {n | 0 ≤ n}

/-- The hypothesis of the conclusion is breakable: `Bfuture` is not invariant. -/
theorem Bfuture_not_invariant :
    ¬ (∀ y : ℤ, y ∈ Bfuture ↔ negZ y ∈ Bfuture) := by
  intro h
  have h1 := h 1
  simp only [Bfuture, Set.mem_setOf_eq, negZ, Equiv.coe_fn_mk] at h1
  norm_num at h1

/-- With a symmetry-breaking boundary, the reversal of an allowed trajectory
need not be allowed — so `arrow_sum_zero`'s hypothesis genuinely fails here. -/
theorem demo_boundary_breaks_reversal :
    Allowed Bfuture (fun n : ℤ => n + 1) ∧
      ¬ Allowed Bfuture (rev negZ (fun n : ℤ => n + 1)) := by
  constructor
  · simp only [Allowed, Bfuture, Set.mem_setOf_eq]; norm_num
  · simp only [Allowed, Bfuture, Set.mem_setOf_eq, rev, negZ, Equiv.coe_fn_mk]
    norm_num

/-! ### Non-vacuity of the arrow predicate itself

`arrow_sum_zero` would be vacuous if no arrow existed. Here is a concrete
**nonzero** arrow: the value at time `0`. It flips sign under reversal purely
because the reversal map is negation. -/

/-- The value of a trajectory at time `0`. -/
def paramArrow : (ℤ → ℤ) → ℤ := fun x => x 0

/-- It is a genuine arrow (and not the zero arrow). -/
theorem paramArrow_isArrow : IsArrow negZ paramArrow := by
  intro x
  simp only [paramArrow, rev, negZ, Equiv.coe_fn_mk]
  ring

/-- The arrow is not identically zero. -/
theorem paramArrow_ne_zero : paramArrow ≠ 0 := by
  intro h
  have := congrFun h (fun n : ℤ => n + 1)
  simp only [paramArrow, Pi.zero_apply] at this
  norm_num at this

/-- **Non-vacuous instance of the conclusion.** A nonzero arrow takes opposite,
nonzero values on a trajectory and its reversal — the cancellation is real, not
trivial. (Stated without a `Finset` because `ℤ → ℤ` has no `DecidableEq`.) -/
theorem demo_nontrivial_arrow_flips :
    paramArrow (fun n : ℤ => n + 1) = 1 ∧
      paramArrow (rev negZ (fun n : ℤ => n + 1)) = -1 ∧
      paramArrow (fun n : ℤ => n + 1) +
        paramArrow (rev negZ (fun n : ℤ => n + 1)) = 0 := by
  have hrev : rev negZ (fun n : ℤ => n + 1) = (fun n : ℤ => n - 1) := by
    funext n; simp only [rev, negZ, Equiv.coe_fn_mk]; ring
  refine ⟨?_, ?_, ?_⟩
  · simp only [paramArrow]; norm_num
  · rw [hrev]; simp only [paramArrow]; norm_num
  · rw [hrev]; simp only [paramArrow]; norm_num

/-- **The closure hypothesis is exactly the dividing line.** The same nonzero
arrow on two two-element families that differ only in whether the family is
closed under reversal: closed → the sum is `0`; not closed → the sum is `2`.
So `arrow_sum_zero` is not a theorem about a trivial arrow. -/
theorem demo_closure_is_the_line :
    (∑ i : Fin 2, paramArrow (if i = 0 then (fun n : ℤ => n + 1)
                            else rev negZ (fun n : ℤ => n + 1))) = 0 ∧
    (∑ _i : Fin 2, paramArrow (fun n : ℤ => n + 1)) = 2 := by
  have hrev : rev negZ (fun n : ℤ => n + 1) = (fun n : ℤ => n - 1) := by
    funext n; simp only [rev, negZ, Equiv.coe_fn_mk]; ring
  constructor
  · rw [Fin.sum_univ_two]
    simp only [Fin.isValue, Fin.zero_eta, ↓reduceIte]
    rw [hrev]
    simp only [paramArrow]
    norm_num
  · rw [Fin.sum_univ_two]
    simp only [paramArrow]
    norm_num

end TimeFlow
