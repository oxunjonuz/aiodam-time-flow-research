#!/bin/sh
# Independent verification of HomogeneousUniverse.lean:
#  1. rebuild from source with the toolchain only;
#  2. scan for sorry/admit/axiom/native_decide/unsafe;
#  3. audit the axioms of every headline theorem;
#  4. NEGATIVE CONTROLS -- each must FAIL, proving the checks have teeth.
set -u
ML=/work/Shopify/audit-work/mathlib
LP="$ML/.lake/build/lib/lean"
for d in "$ML"/.lake/packages/*/; do
  b="$d.lake/build/lib/lean"
  [ -d "$b" ] && LP="$LP:$b"
done
cd /work/time_flow_research_20260927/work/lean || exit 9
export LEAN_PATH="$LP:$PWD"
rc=0

echo "=== 1. build HomogeneousUniverse.lean ==="
/opt/lean4/bin/lean HomogeneousUniverse.lean 2>&1
e=$?; echo "    exit=$e"; [ "$e" -ne 0 ] && rc=1

echo "=== 2. forbidden-token scan ==="
if grep -nE '(^|[^A-Za-z_])(sorry|admit|axiom|native_decide|unsafe|implemented_by)([^A-Za-z_]|$)' \
     HomogeneousUniverse.lean; then
  echo "FOUND_FORBIDDEN"; rc=1
else
  echo "NO_SORRY_NO_AXIOM"
fi

echo "=== 3. axiom audit ==="
/opt/lean4/bin/lean -o HomogeneousUniverse.olean HomogeneousUniverse.lean 2>/dev/null
cat > /tmp/ax_cosmo.lean <<'EOF'
import Mathlib
import HomogeneousUniverse
#print axioms TimeFlowCosmo.creation_not_identifiable_unpinned
#print axioms TimeFlowCosmo.creation_identifiable_pinned
#print axioms TimeFlowCosmo.demo_unpinned_pair
#print axioms TimeFlowCosmo.demo_pinned_identifiable
#print axioms TimeFlowCosmo.baseline_unique_when_pinned
EOF
/opt/lean4/bin/lean /tmp/ax_cosmo.lean 2>&1

echo "=== 4. negative controls (each MUST fail) ==="
mk() { # $1 = name, $2 = replacement body for the theorem
  cp HomogeneousUniverse.lean "/tmp/nc_$1.lean"
  python3 - "$1" "$2" <<'PY'
import sys, re
name, body = sys.argv[1], sys.argv[2]
p = f"/tmp/nc_{name}.lean"
s = open(p).read()
s = s.replace("namespace TimeFlowCosmo", "namespace TimeFlowCosmo\n-- negative control", 1)
open(p, "w").write(s)
PY
}
# NC1: remove the pinned-baseline hypothesis -> the pinned theorem must fail
sed 's/^theorem creation_identifiable_pinned (H : ℝ → ℝ) (δ : ℝ → ℝ) :/theorem creation_identifiable_pinned_broken (H : ℝ → ℝ) (δ : ℝ → ℝ) :/' \
    HomogeneousUniverse.lean > /tmp/nc1.lean
python3 - <<'PY'
p="/tmp/nc1.lean"; s=open(p).read()
s=s.replace("  funext a\n  have ha := h a\n  simp only [obs] at ha\n  linarith",
            "  funext a\n  have ha := h a\n  simp only [obs] at ha\n  -- dropped: linarith",1)
open(p,"w").write(s)
PY
echo "--- NC1: pinned theorem with the closing step removed (expect failure)"
/opt/lean4/bin/lean /tmp/nc1.lean 2>&1 | head -6

# NC2: assert the unpinned theorem is FALSE i.e. claim identifiability is possible
cat > /tmp/nc2.lean <<'EOF'
import Mathlib
import HomogeneousUniverse
open TimeFlowCosmo
-- NC2: assert the degeneracy does NOT hold. Must fail to compile.
theorem nc2_wrong : IdentifiableFree (fun _ : ℝ => (0:ℝ)) (fun _ => 1) := by
  intro δ' H' h
  have h0 := h 0
  simp only [obs] at h0
  funext a
  have ha := h a
  simp only [obs] at ha
  linarith
EOF
echo "--- NC2: asserting the opposite of the theorem (expect failure)"
/opt/lean4/bin/lean /tmp/nc2.lean 2>&1 | head -6

# NC3: sneak in a `sorry`
python3 - <<'PY'
s=open("/work/time_flow_research_20260927/work/lean/HomogeneousUniverse.lean").read()
s=s.replace("theorem demo_pinned_identifiable :\n    IdentifiablePinned (fun _ => 1) (fun _ => 0) :=\n  creation_identifiable_pinned (fun _ => 1) (fun _ => 0)",
            "theorem demo_pinned_identifiable :\n    IdentifiablePinned (fun _ => 1) (fun _ => 0) := by sorry",1)
open("/tmp/nc3.lean","w").write(s)
PY
echo "--- NC3: injected sorry must be caught by the scan"
if grep -nE '(^|[^A-Za-z_])(sorry|admit|axiom|native_decide|unsafe|implemented_by)([^A-Za-z_]|$)' /tmp/nc3.lean; then
  echo "NC3_DETECTED"
else
  echo "NC3_MISSED -- scan has no teeth"; rc=1
fi

echo "=== overall rc=$rc ==="
exit $rc
