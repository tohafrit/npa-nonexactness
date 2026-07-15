"""Float-free verifier: NPA level 3 is NOT exact at s=1/20 for the doubly-tilted
CHSH functional. Extends the certified ladder (level 1+AB at s=1/10, level 2 at
s=1/5) to level 3 -- concrete evidence that "no finite NPA level is exact" does not
stop at level 2.

Certificate (cert_subcrit_L3.json): a rational moment vector y over the level-3
word set (ScenarioL(3), 25x25 moment matrix) with
  (i)   Gamma_3(y) exactly positive semidefinite (rational LDL) and y_1 = 1
        -- a valid level-3 pseudo-moment (feasible for NPA level 3);
  (ii)  objective B_s.y (exact rational) strictly exceeds the quantum value c_Q(s),
        the largest real root of Gigena et al.'s exact sextic, verified by Sturm.
Hence c_3(1/20) >= B_s.y > c_Q(1/20): level 3 overshoots. QED.
Requires sympy (exact Sturm root count).
"""
import json
from fractions import Fraction
import sympy as sp
from scenario_ext import ScenarioL

def psd_exact(M):
    n = len(M); A = [r[:] for r in M]
    for k in range(n):
        if A[k][k] < 0: return False
        if A[k][k] == 0:
            if any(A[k][j] != 0 for j in range(k, n)): return False
            continue
        for i in range(k+1, n):
            f = A[i][k]/A[k][k]
            for j in range(k, n): A[i][j] -= f*A[k][j]
    return True

def main():
    s3 = ScenarioL(3); N, idx, words = s3.N, s3.idx, s3.words
    def W(a, b): return words[(a, b)]
    al = sp.Rational(39, 40)  # 1 - s/2, s = 1/20
    d = json.load(open("cert_subcrit_L3.json"))
    yr = [Fraction(n, dd) for n, dd in zip(d["y_num"], d["y_den"])]
    # (i) feasibility at level 3
    G = [[yr[idx[i, j]] for j in range(N)] for i in range(N)]
    feas = psd_exact(G) and yr[W((), ())] == 1
    # (ii) objective (exact) and c_Q comparison via Sturm
    ys = [sp.Rational(x.numerator, x.denominator) for x in yr]
    obj = (al*ys[W((0,), ())] + al*ys[W((), (0,))] + ys[W((0,), (0,))]
           + ys[W((0,), (1,))] + ys[W((1,), (0,))] - ys[W((1,), (1,))])
    lam = sp.Symbol('lam'); a = b = al
    CO = [sp.Integer(4), -4*a*b, 11*a*a*b*b-16*a*a-16*b*b-64,
          8*a**3*b**3-24*a**3*b-24*a*b**3+96*a*b,
          -6*a**4*b*b-6*a*a*b**4-64*a*a*b*b+20*a**4+96*a*a+20*b**4+96*b*b+320,
          -168*a**3*b**3+60*a**5*b+160*a**3*b+60*a*b**5+160*a*b**3-576*a*b,
          a**6*(27*b*b-8)+a**4*(-54*b**4+48*b*b+32)
          + a*a*(27*b**6+48*b**4-400*b*b+128)-8*b**6+32*b**4+128*b*b-512]
    S = sp.Poly(sp.expand(sum(c*lam**(6-k) for k, c in enumerate(CO))), lam)
    roots_above = S.count_roots(obj, sp.oo)  # exact Sturm
    beats = (roots_above == 0)
    print(f"(i) Gamma_3(y) exact PSD (25x25) + normalized: {feas}")
    print(f"(ii) objective (exact) = {float(obj):.12f}")
    print(f"     c_Q(1/20) largest sextic root ~ {float(max(r for r in sp.Poly(S,lam).nroots() if abs(sp.im(r))<1e-9).as_real_imag()[0]):.12f}")
    print(f"     sextic real roots above objective (Sturm): {roots_above}  (0 => obj > c_Q)")
    ok = feas and beats
    print("\nLEVEL 3 NOT EXACT AT s=1/20 VERIFIED (float-free)" if ok
          else "\nVERIFICATION FAILED")
    return 0 if ok else 1

if __name__ == "__main__":
    raise SystemExit(main())
