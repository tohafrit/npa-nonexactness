"""Independent float-free verifier for the SUPERCRITICAL EXACTNESS THEOREM
(npa_overshoot_note.tex, Theorem: phase transition).

Checks, in exact rational arithmetic, that Z0, Z1 (level 1+AB, 9x9) satisfy:
  (a) Z0, Z1 symmetric and PSD (rational LDL);
  (b) affine identity coefficients:
        4 - B_11 . y = <Z0, Gamma(y)>   for all y   (B_11 = <A0>+<B0>+CHSH-part)
        2 - <A0> - <B0> = <Z1, Gamma(y)> for all y
      i.e. for every word class w: class-sum(Z0) = target0_w, class-sum(Z1)=target1_w.
Consequently Z(nu)=Z0+nu Z1 >= 0 certifies B_{1+nu,1+nu}.y <= 4+2nu at level
1+AB, hence at every NPA level: the hierarchy is EXACT on the whole
supercritical ray, at all levels.
"""
import json
from fractions import Fraction
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
            for j in range(k, n):
                A[i][j] -= f*A[k][j]
    return True

def load(fn):
    d = json.load(open(fn))
    N = len(d["Z_num"])
    return [[Fraction(d["Z_num"][p][q], d["Z_den"][p][q]) for q in range(N)]
            for p in range(N)]

def main():
    s1 = Scenario(2, 2); N, idx = s1.N, s1.words
    cls = {}
    for p in range(N):
        for q in range(N):
            cls.setdefault(s1.idx[p, q], []).append((p, q))
    # B_11 . y = <A0>+<B0> + C00+C01+C10-C11 ; identity 4 - B11.y = <Z0,Gamma>
    B0v = {idx[((0,),())]: 1, idx[((),(0,))]: 1, idx[((0,),(0,))]: 1,
           idx[((0,),(1,))]: 1, idx[((1,),(0,))]: 1, idx[((1,),(1,))]: -1}
    tgt0 = {w: (Fraction(4) if w == idx[((),())] else Fraction(-B0v.get(w, 0)))
            for w in range(s1.NV)}
    tgt1 = {w: Fraction(0) for w in range(s1.NV)}
    tgt1[idx[((),())]] = Fraction(2)
    tgt1[idx[((0,),())]] = Fraction(-1)
    tgt1[idx[((),(0,))]] = Fraction(-1)

    # Z_A: 1 - <A0> = <Z_A, Gamma(y)>;  Z_B: 1 - <B0> = <Z_B, Gamma(y)>
    tgtA = {w: Fraction(0) for w in range(s1.NV)}
    tgtA[idx[((), ())]] = Fraction(1); tgtA[idx[((0,), ())]] = Fraction(-1)
    tgtB = {w: Fraction(0) for w in range(s1.NV)}
    tgtB[idx[((), ())]] = Fraction(1); tgtB[idx[((), (0,))]] = Fraction(-1)

    ok = True
    for fn, tgt, name in (("cert_Z0_critical.json", tgt0, "Z0"),
                          ("cert_Z1_supercritical.json", tgt1, "Z1"),
                          ("cert_ZA.json", tgtA, "ZA"),
                          ("cert_ZB.json", tgtB, "ZB")):
        Z = load(fn)
        sym = all(Z[p][q] == Z[q][p] for p in range(N) for q in range(N))
        sums = all(sum(Z[p][q] for p, q in cls[w]) == tgt[w] for w in range(s1.NV))
        psd = psd_exact(Z)
        print(f"{name}: symmetric {sym} | affine-identity class-sums {sums} | PSD {psd}")
        ok &= sym and sums and psd
    print("\n=> Z(alpha,beta) = Z0 + (alpha-1) ZA + (beta-1) ZB satisfies")
    print("   (2+alpha+beta) - B_{alpha,beta}.y = <Z(alpha,beta), Gamma(y)>  (by linearity),")
    print("   and for alpha,beta >= 1 it is a NONNEGATIVE combination of the exactly-PSD")
    print("   Z0, ZA, ZB, hence PSD. Therefore the NPA hierarchy is exact at EVERY level")
    print("   for all alpha,beta >= 1 (the symmetric ray alpha=beta=1+nu is the diagonal).")
    print("\nSUPERCRITICAL EXACTNESS (quadrant alpha,beta>=1) VERIFIED (float-free)" if ok
          else "\nVERIFICATION FAILED")
    return 0 if ok else 1

if __name__ == "__main__":
    raise SystemExit(main())
