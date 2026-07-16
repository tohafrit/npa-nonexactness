# No finite level of the NPA hierarchy is exact for the doubly-tilted CHSH functional

**Anton Pakhunov** — pakhunov.anton.n@gmail.com

This repository accompanies the paper

> *No finite level of the NPA hierarchy is exact for the doubly-tilted CHSH
> functional near the critical tilt*, arXiv:2607.13762,
> `paper/npa_nonexactness_proof.pdf`

which resolves the open question of Gigena, Panwar, Scala, Araújo, Farkas and
Chaturvedi [npj Quantum Inf. **11**, 82 (2025)]: for the critical family
B_s = (1−s/2)(⟨A0⟩+⟨B0⟩) + CHSH, **no finite level of the
Navascués–Pironio–Acín hierarchy is exact on any neighbourhood of the
critical point**. For every level k ≥ 2 there are certified rational g_k > 0
and s*_k > 0 with c_k(s) ≥ 4 − s + g_k s² on (0, s*_k], while the quantum
value is cubically flat. Unconditionally a_2 > 1/39, a_3 > 1/188, a_4 > 1/641.

The companion note (arXiv:2607.13774, `companion/`) develops the phase-transition geometry the
proof is built on (supercritical exactness theorem, Motzkin scaling limit,
the exact almost-quantum coefficient 3/64, the certified non-exactness
ladder).

## Reproducing every machine check

The proof is computer-assisted in the strict sense: finite exact-integer or
exact-rational verifications are constituent parts of the argument. No
floating point enters any proof-critical check.

```sh
sh run_all_checks.sh        # all 12 verifiers, ~40 s; exits 0 iff all pass
```

Requirements: Python ≥ 3.10 with `sympy` and `numpy` (see `VERSIONS.txt` for
the environment used for the reported runs).

| verifier | what it proves / checks |
|---|---|
| `verify_pipeline_laws.py` | structural laws: dual certificate, face/product law, kernel law, uniform Slater point (k ≤ 8) |
| `verify_r1_induction.py` | the witness identity N^T Γ(y*) N = u u^T: exhaustive exact-integer grid, 1,413,721 kernel pairs (side ≤ 20; margin ≤ 24/≤ 26) plus per-level checks k = 2..8 |
| `verify_witness_symbolic.py` | the same identity proven **symbolically**: all 58,081 fine regimes, E ≡ 0 as a polynomial — no interpolation argument needed |
| `cleanroom_check.py` | independent clean-room implementation written from the paper text alone (stdlib only), 55 checks incl. the witness identity k = 2..7 |
| `verify_bridge_lemma.py` | the swap-split bridge lemma (k = 2..5) |
| `verify_small_lambda.py` | the small-λ feasible point and its positive value (k = 2, 3, 4) |
| `verify_framework_lemma.py` | exact arc feasibility (no-repair lemma): exact rational PD along the arc (k = 2, 3, 4) |
| `verify_sstar_rational.py` | certified **rational** interval endpoints s*_2 ≥ 1.64e-3, s*_3 ≥ 9.29e-3, s*_4 ≥ 7.55e-3 |
| `verify_v{2,3,4}_rational_base.py` | the certified ladder: a_2 > 1/39, a_3 > 1/188, a_4 > 1/641 (exact rational certificates) |
| `audit_spot.py` | independent re-verification pass (fresh implementations of the core objects; composition checks) |
| `companion/anc/` | the companion note's certificates: supercritical exactness, non-exactness ladder (Sturm-verified), a_{1+AB} = 3/64 exact |

`scenario_ext.py`, `npa_general.py` are the shared word/moment-matrix library.

## Citing

    A. Pakhunov, "No finite level of the NPA hierarchy is exact for the
    doubly-tilted CHSH functional near the critical tilt", arXiv:2607.13762.
    A. Pakhunov, "A phase transition in the exactness of the NPA hierarchy
    at the critical doubly-tilted CHSH functional", arXiv:2607.13774.
