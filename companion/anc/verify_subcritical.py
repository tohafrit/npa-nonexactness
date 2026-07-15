"""Float-free verifier: NPA level 1+AB is NOT exact at s=1/10 for the
doubly-tilted CHSH functional. Complements the exact supercritical theorem.

Certificate (cert_subcrit_1AB.json): a rational moment vector y with
  (i)   Gamma(y) exactly positive semidefinite (rational LDL) and y_1 = 1
        -- so y is a valid level-1+AB pseudo-moment (feasible for NPA);
  (ii)  objective B_{s}.y (exact rational) strictly exceeds the quantum value
        c_Q(s), the largest real root of Gigena et al.'s exact sextic, verified
        by Sturm's theorem (no root of the sextic lies above the objective).
Hence c_{1+AB}(1/10) >= B_s.y > c_Q(1/10): level 1+AB overshoots. QED.
Requires sympy (exact Sturm root count).
"""
import json
from fractions import Fraction
import sympy as sp
from npa_general import Scenario

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
    s1 = Scenario(2, 2); N, idx = s1.N, s1.words
    al = sp.Rational(19, 20)  # 1 - s/2, s = 1/10
    d = json.load(open("cert_subcrit_1AB.json"))
    yr = [Fraction(n, dd) for n, dd in zip(d["y_num"], d["y_den"])]
    # (i) feasibility
    G = [[yr[s1.idx[i, j]] for j in range(N)] for i in range(N)]
    feas = psd_exact(G) and yr[idx[((), ())]] == 1
    # (ii) objective (exact) and c_Q comparison via Sturm
    ys = [sp.Rational(x.numerator, x.denominator) for x in yr]
    obj = (al*ys[idx[((0,),())]] + al*ys[idx[((),(0,))]] + ys[idx[((0,),(0,))]]
           + ys[idx[((0,),(1,))]] + ys[idx[((1,),(0,))]] - ys[idx[((1,),(1,))]])
    lam = sp.Symbol('lam'); a = b = al
    CO = [sp.Integer(4), -4*a*b, 11*a*a*b*b-16*a*a-16*b*b-64,
          8*a**3*b**3-24*a**3*b-24*a*b**3+96*a*b,
          -6*a**4*b*b-6*a*a*b**4-64*a*a*b*b+20*a**4+96*a*a+20*b**4+96*b*b+320,
          -168*a**3*b**3+60*a**5*b+160*a**3*b+60*a*b**5+160*a*b**3-576*a*b,
          a**6*(27*b*b-8)+a**4*(-54*b**4+48*b*b+32)
          + a*a*(27*b**6+48*b**4-400*b*b+128)-8*b**6+32*b**4+128*b*b-512]
    S = sp.Poly(sp.expand(sum(c*lam**(6-k) for k, c in enumerate(CO))), lam)
    roots_above = S.count_roots(obj, sp.oo)  # exact Sturm
    beats = (roots_above == 0)  # obj exceeds the largest root c_Q
    print(f"(i) Gamma(y) exact PSD + normalized: {feas}")
    print(f"(ii) objective (exact) = {obj} = {float(obj):.10f}")
    print(f"     sextic real roots above objective (Sturm): {roots_above}  (0 => obj > c_Q)")
    ok = feas and beats
    print("\nLEVEL 1+AB NOT EXACT AT s=1/10 VERIFIED (float-free)" if ok
          else "\nVERIFICATION FAILED")
    return 0 if ok else 1

if __name__ == "__main__":
    raise SystemExit(main())
