#!/bin/sh
# Reproduces every machine check behind:
#   "No finite level of the NPA hierarchy is exact for the doubly-tilted
#    CHSH functional near the critical tilt" (npa_nonexactness_proof.tex)
# and its companion note. Exits 0 iff every verifier passes.
# Requirements: Python 3 (>=3.10) with sympy and numpy. Only the ladder
# verifiers (verify_v*_rational_base.py) additionally exercise fractions-
# based exact linear algebra; no floating point enters any load check.
set -e
for v in verify_pipeline_laws verify_r1_induction verify_witness_symbolic \
         verify_bridge_lemma verify_small_lambda verify_framework_lemma \
         verify_sstar_rational verify_v2_rational_base \
         verify_v3_rational_base verify_v4_rational_base cleanroom_check \
         audit_spot; do
  echo "== $v =="
  python3 "$v.py"
done
echo "ALL VERIFIERS PASSED"
