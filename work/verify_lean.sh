#!/bin/sh
# Build the Lean formalization from scratch and check for sorry/axiom.
# Independent path: uses only the toolchain and the .lean sources on disk.
set -u
ML=/work/Shopify/audit-work/mathlib
LP="$ML/.lake/build/lib/lean"
for d in "$ML"/.lake/packages/*/; do
  b="$d.lake/build/lib/lean"
  [ -d "$b" ] && LP="$LP:$b"
done
cd /work/time_flow_research_20260927/work/lean || exit 9
export LEAN_PATH="$LP:$PWD"
echo "=== LEAN_PATH entries: $(echo "$LP" | tr ':' '\n' | wc -l) ==="
rc=0
for f in Identifiability.lean RelativityOfNow.lean Instances.lean PastHypothesis.lean HomogeneousUniverse.lean; do
  echo "--- building $f"
  /opt/lean4/bin/lean "$f" 2>&1
  e=$?
  echo "    exit=$e"
  [ "$e" -ne 0 ] && rc=1
done
echo "=== sorry/axiom scan ==="
if grep -nE '(^|[^A-Za-z_])(sorry|admit|axiom)([^A-Za-z_]|$)' \
     Identifiability.lean RelativityOfNow.lean Instances.lean PastHypothesis.lean HomogeneousUniverse.lean; then
  echo "FOUND_FORBIDDEN"; rc=1
else
  echo "NO_SORRY_NO_AXIOM"
fi
echo "=== axiom audit of the headline theorems ==="
# `lean <file>.lean` type-checks but does not emit oleans; emit them so the
# audit can import the modules.
for f in Identifiability RelativityOfNow Instances PastHypothesis HomogeneousUniverse; do
  /opt/lean4/bin/lean -o "$f.olean" "$f.lean" 2>/dev/null
done
cat > /tmp/ax_audit.lean <<'EOF'
import Mathlib
import Identifiability
import RelativityOfNow
import PastHypothesis
import HomogeneousUniverse
#print axioms TimeFlow.identifiable_iff
#print axioms TimeFlow.not_identifiable_of_range_le
#print axioms TimeFlow.exists_boost_reversing_time_order
#print axioms TimeFlow.rev_isTraj
#print axioms TimeFlow.arrow_vanishes_on_symmetric
#print axioms TimeFlow.arrow_sum_zero
#print axioms TimeFlow.demo_arrow_vanishes
#print axioms TimeFlow.demo_boundary_breaks_reversal
#print axioms TimeFlow.paramArrow_isArrow
#print axioms TimeFlow.paramArrow_ne_zero
#print axioms TimeFlow.demo_nontrivial_arrow_flips
#print axioms TimeFlow.demo_closure_is_the_line
#print axioms TimeFlowCosmo.creation_not_identifiable_unpinned
#print axioms TimeFlowCosmo.creation_identifiable_pinned
#print axioms TimeFlowCosmo.demo_unpinned_pair
#print axioms TimeFlowCosmo.demo_pinned_identifiable
#print axioms TimeFlowCosmo.baseline_unique_when_pinned
EOF
/opt/lean4/bin/lean /tmp/ax_audit.lean 2>&1
echo "=== overall rc=$rc ==="
exit $rc