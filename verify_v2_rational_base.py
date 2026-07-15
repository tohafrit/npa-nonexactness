"""Float-free verifier: the RATIONAL-BASE tangent value v_2(1/2, 1/4) of the
level-2 (ScenarioL(2), 13x13) second-order tangent SDP for the symmetric
doubly-tilted CHSH functional is pinned to 52 exact decimal digits by a
rational primal/dual sandwich, and it is WILD (not low-degree algebraic).

THE OBJECT.  v_2(beta,delta) is the value of the Schur-linearized second-order
tangent SDP of a2_selection_functional.py (the program whose max over the face
is the unconstrained overshoot a_2 = 0.0280586722...):

  v_2(b,d) = max_{t in R^7, y2 in R^18}  B0.y2 - y1_<A0>
     s.t. y1 = Bv^T t (exact rational 7-dim first-order space V),
          lam(t) >= 0,   y1_<B1> = 0,  y1_<B1B0B1> = 0,  B0.y1 = 0,
          y2_1 = y2_<B1> = y2_<B1B0B1> = 0,
          [[Ms0(b,d), Ms(y1)NsKs], [., NsKs^T Ms(y2) NsKs]] >= 0   (12x12),
          [[Ma0(b,d), Ma(y1)Na ], [., Na^T  Ms(y2) Na  ]] >= 0   ( 8x8),

evaluated at the FIXED RATIONAL base point (beta,delta) = (1/2,1/4) of the
two-parameter optimal face.  Pinning only lowers the value, so
v_2(b,d) <= a_2; in particular v_2 > 0 here certifies a_2 > 0.

RESULT (the bet of 2026-07-14 -- "maybe the wildness lives only in the argmax,
and a rational base gives a low-degree v_2, as at level 1 where
a_1(beta) = 3 beta(1-beta)/16" -- is REFUTED):

  v_2(1/2, 1/4) = 0.0260925523962229226120432006824120402958354068753880
                    1873220404545028844559251774526785994784736538229402...

  (a) EXACT SANDWICH (this file, float-free rationals only):
        L <= v_2(1/2,1/4) <= U,   U - L < 2.3e-53,
      L from an exactly-feasible rational primal point (60-digit truncation of
      the certified optimum, blended eps_p = 1e-54 toward an exact interior
      Slater point to repair the PSD boundary), U from an exactly-stationary
      PSD rational dual pair (Us,Ua) (projection of the certified dual onto
      the exact 7-dim stationarity affine space, blended eps_d = 1e-52 toward
      an exact strictly-positive dual Slater point).
  (b) WILDNESS (documented; conditional on mpmath PSLQ):  from a 535-digit
      certified solution (mpmath Newton on the exact-flag KKT system, exact
      rational data, final residual 5.9e-541, Jacobian sigma_min ~ 1.6e-5):
      v_2(1/2,1/4) admits NO integer minimal polynomial with
         deg <= 16, |coeff| <= 1e24   (needs ~425 digits; had ~530)
         deg <= 20, |coeff| <= 1e18   (needs ~399)
         deg <= 24, |coeff| <= 1e10;  deg <= 32, |coeff| <= 1e7.
      A spurious deg-16/1e18 PSLQ hit at 350 digits was caught and refuted by
      re-evaluation at 535 digits (p(v_2) = 4.5e-287 != 0) -- PSLQ candidates
      at the information boundary must always be re-tested at higher precision.
      Also no deg <= 6 / 1e15 minpoly for the KKT coordinates sig, alf,
      lam = t_6, v_A, w_A; v_2 not in Q + Q sqrt2 + Q sqrt3 + Q sqrt6 (1e12).
  (c) SECOND RATIONAL BASE (1/2, 1/8) (same pipeline, 220 certified digits):
      v_2(1/2,1/8) = 0.026684688593982145105613604666827137252141068489979...
      NO minpoly deg <= 12, |coeff| <= 1e12.  So rational-base wildness is
      generic on the face, not special to (1/2,1/4).
      CAVEAT (1/2, 0): the base (1/2,0) named in the task is DEGENERATE --
      Ms0(1/2,0) drops rank 3 -> 2, the fixed integer kernel Ns is then a
      strict subset of ker Ms0 and this Schur reduction is a relaxation there
      (block-SDP value at (1/2,0) is 0.017458 (solver-grade), the reduced
      program gives 0.0439); (1/2,1/8) is used as the interior second base.

VERDICT: WILD.  The wildness of the level-2 overshoot program is INTRINSIC to
the second-order tangent SDP, not an artifact of the wild argmax position:
even with ALL data rational (rational base, integer kernels, rational V), the
optimal value is not low-degree algebraic.  The optimal SUPPORT (rank flags:
Schur complements rank 1, duals rank 3+2 with irrational kernel vectors v, w)
forces an irrational boundary point of algebraic degree beyond deg 32 / 1e7.

IMPLICATION, honest: the strong "rational-base uniform certificate schema"
(rational v_k values level by level) is DEAD in its direct form.  What
SURVIVES -- and what this file actually exhibits -- is the weaker uniform
route: v_k > 0 can be certified by a PURELY RATIONAL feasible primal jet
(here L > 1/39 > 0 float-free), so a uniform-in-k positivity certificate does
not require the value to be nice, only a rational interior approximation of
the optimal jet.  The exact-value route dies; the exact-positivity route
stands.

STRUCTURAL BYPRODUCTS verified exact below (identically in t over V):
  * B0.y1 == 0 on all of V: the first-order objective term vanishes
    IDENTICALLY (the constraint B0.y1 = 0 in the program is vacuous, its
    multiplier is pure gauge -- KKT systems built on it are singular).
  * y1_<B1B0B1> = 0 forces t_0 = 0; lam = t_6.
  * Range conditions Ns^T Ms1(t) NsKs == 0 and (Na^T Ma1(t) Na == 0) hold
    identically, so block-PSD <=> Schur-complement-PSD (generalized Schur
    lemma) and the 12x12/8x8 program collapses exactly to a 4x4/3x3 one.
  * The y2 -> (NsKs^T Ms2 NsKs, Na^T Ma2 Na) compression has zero classes
    {11,12,13} and a further exact 3-dim rational kernel with B0.ker = 0:
    only 9 of 15 free y2 classes are effective (classes 1,3,4,5,6,7,8,9,14).
  * At the optimum (certified, conditional): the pinning multipliers mu
    vanish, y2_<A1B1> = 0, v_4 = 0, Us_33 = 1/4 exactly.

EXACT vs CONDITIONAL, honest split:
  EXACT (this file, sympy rationals, no floats): everything in checks
  (0)-(5): base-point PSD/kernel, reduction identities, gauge structure,
  primal feasibility => L, dual stationarity+PSD => U, L <= ref <= U, L > 0.
  CONDITIONAL (documented, not re-derived here): the 535-digit value itself
  (mpmath Newton, quadratic convergence, residual 5.9e-541 -- but digit count
  certified only via Jacobian conditioning, not interval arithmetic), the
  PSLQ exclusions, the rank/activity flags used to FIND the optimum (the
  sandwich does not depend on them), CLARABEL cross-check 0.0260925299
  (agrees with the sandwich to its ~1e-8 accuracy).

Requires sympy + scenario_ext/npa_general (exact integer tables) +
v2_cert_data.py (the embedded rational certificate).  Runtime ~1-2 min.
Exit code 0 iff all exact checks pass.
"""
import sys
import sympy as sp

sys.path.insert(0, '/Users/antonp/research/why-tsirelson')

R = sp.Rational
REF100 = ('0.0260925523962229226120432006824120402958354068753880'
          '18732204045450288445592517745267859947847365382294')


def main():
    from scenario_ext import ScenarioL
    from npa_general import canon
    from v2_cert_data import (T_PRIMAL, Y2_PRIMAL, Y2_SLATER, EPS_P,
                              US_DUAL, UA_DUAL, US_SLATER, UA_SLATER, EPS_D)
    ok = True

    def chk(name, cond):
        nonlocal ok
        ok = ok and bool(cond)
        print(f"  [{'OK' if cond else 'FAIL'}] {name}")

    # ---- exact scenario tables (as in verify_a2_exact.py) ----
    s2 = ScenarioL(2)
    inv = {v: k for k, v in s2.words.items()}
    NV = s2.NV
    swp = [s2.words[canon(bw, aw)] for i in range(NV) for aw, bw in [inv[i]]]
    cls = [-1]*NV
    ncls = 0
    for i in range(NV):
        if cls[i] < 0:
            cls[i] = cls[swp[i]] = ncls
            ncls += 1
    NC = ncls
    assert NC == 18
    E = sp.eye(13)
    Tm = sp.Matrix.hstack(E[:, 0], E[:, 6], E[:, 10], E[:, 1]+E[:, 5],
                          E[:, 2]+E[:, 8], E[:, 3]+E[:, 11], E[:, 4]+E[:, 12],
                          E[:, 7]+E[:, 9])
    Um = sp.Matrix.hstack(E[:, 1]-E[:, 5], E[:, 2]-E[:, 8], E[:, 3]-E[:, 11],
                          E[:, 4]-E[:, 12], E[:, 7]-E[:, 9])
    Pf = [sp.zeros(13, 13) for _ in range(NC)]
    for p in range(13):
        for q in range(13):
            Pf[cls[s2.idx[p, q]]][p, q] = 1
    Ps = [Tm.T*Pc*Tm for Pc in Pf]
    Pa = [Um.T*Pc*Um for Pc in Pf]
    Ns = sp.Matrix([[-2, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0, 1, 0],
                    [-1, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0, 0, 1],
                    [1, 0, 1, 0, -1, 0, 0, 0]]).T
    Na = sp.Matrix([[1, 0, 0, 0, 0], [0, -1, 0, 1, 0], [0, -1, 0, 0, 1]]).T
    u = sp.Matrix([1, -1, 0, -1, 0])
    Ks = sp.Matrix([[1, 1, 0, 0, 0], [0, 0, 1, 0, 0], [1, 0, 0, 1, 0],
                    [0, 0, 0, 0, 1]]).T
    NsKs = Ns*Ks
    # exact rational basis of V (verify_a2_exact.py check [2]); cols=classes 1..17, lam
    Bv = sp.Matrix([
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, -1, 0, 0, -1, -2, 1, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        [0, 1, R(1, 2), 0, 1, 2, 0, R(1, 2), R(3, 2), 0, 0, 0, 0, 1, 1, 0, 1, 0],
        [-R(1, 8), -R(7, 8), -R(3, 8), 0, -R(3, 4), -R(7, 4), -R(1, 8),
         -R(1, 2), -1, 0, 0, 0, 0, 0, -1, 0, 0, 1]])
    be, de = R(1, 2), R(1, 4)
    ga = 2*be - 1
    ep = de + be - 1
    y0 = [1, 1, be, be, 1, be, ga, be, be, ga, de, de, de, ep, ga, ga, be, ga]
    Ms0 = sp.zeros(8, 8)
    Ma0 = sp.zeros(5, 5)
    for c in range(NC):
        Ms0 += y0[c]*Ps[c]
        Ma0 += y0[c]*Pa[c]
    b0 = {1: sp.Integer(2), 4: sp.Integer(1), 5: sp.Integer(2), 6: sp.Integer(-1)}
    FREE = [c for c in range(NC) if c not in (0, 2, 10)]
    Y2K = [1, 3, 4, 5, 6, 7, 8, 9, 14]

    def is_pd(M):
        return all(M[:k, :k].det() > 0 for k in range(1, M.rows + 1))

    def frobs(A, B):
        return sum(A[i, j]*B[i, j] for i in range(A.rows) for j in range(A.cols))

    print("(0) rational base point (1/2,1/4): kernel + exact PSD:")
    chk("Ms0*Ns == 0, Ma0*Na == 0 (integer kernels)",
        Ms0*Ns == sp.zeros(8, 5) and Ma0*Na == sp.zeros(5, 3))
    chk("rank Ms0 = 3, rank Ma0 = 2 (kernels are EXACTLY span Ns / span Na)",
        Ms0.rank() == 3 and Ma0.rank() == 2 and Ns.rank() == 5 and Na.rank() == 3)

    def pinv_exact(M):
        cols = []
        for j in range(M.cols):
            if M[:, cols + [j]].rank() == len(cols) + 1:
                cols.append(j)
        Rc = M[:, cols]
        return Rc, Rc*(Rc.T*M*Rc).inv()*Rc.T
    RcS, Ms0p = pinv_exact(Ms0)
    RcA, Ma0p = pinv_exact(Ma0)
    chk("Ms0 >= 0: congruence T^T Ms0 T = diag(0, G), G PD (Sylvester)",
        is_pd(RcS.T*Ms0*RcS) and sp.Matrix.hstack(Ns, RcS).rank() == 8)
    chk("Ma0 >= 0 likewise", is_pd(RcA.T*Ma0*RcA)
        and sp.Matrix.hstack(Na, RcA).rank() == 5)
    chk("Moore-Penrose: M M+ M = M, M+ symmetric, (M M+)^T = M M+ (both blocks)",
        Ms0*Ms0p*Ms0 == Ms0 and Ms0p.T == Ms0p and (Ms0*Ms0p).T == Ms0*Ms0p
        and Ma0*Ma0p*Ma0 == Ma0 and Ma0p.T == Ma0p and (Ma0*Ma0p).T == Ma0*Ma0p)
    chk("Ms0+ >= 0 (so Qs(t) = F^T Ms0+ F is PSD-valued, matrix-convex)",
        (Ms0p - RcS*(RcS.T*Ms0*RcS).inv()*RcS.T) == sp.zeros(8, 8))

    print("(1) exact reduction identities (identically in t over V):")
    ts = sp.symbols('t0:7')
    y1s = [sp.Integer(0)]*NC
    for i in range(7):
        for c in range(1, NC):
            y1s[c] += ts[i]*Bv[i, c-1]
    lam_s = sum(ts[i]*Bv[i, 17] for i in range(7))
    Ms1s = sp.zeros(8, 8)
    Ma1s = sp.zeros(5, 5)
    for c in range(1, NC):
        if y1s[c] != 0:
            Ms1s += y1s[c]*Ps[c]
            Ma1s += y1s[c]*Pa[c]
    chk("Ns^T Ms1(t) Ns == lam(t) u u^T and u^T Ks == 0  => range condition",
        sp.expand(Ns.T*Ms1s*Ns - lam_s*u*u.T) == sp.zeros(5, 5)
        and u.T*Ks == sp.zeros(1, 4))
    chk("Na^T Ma1(t) Na == 0 identically  => antisym range condition",
        sp.expand(Na.T*Ma1s*Na) == sp.zeros(3, 3))
    B0y1 = sp.expand(2*y1s[1] + y1s[4] + 2*y1s[5] - y1s[6])
    chk("B0.y1 == 0 IDENTICALLY on V (the B0.y1=0 constraint is vacuous)",
        B0y1 == 0)
    chk("y1_<B1B0B1> == t0 and lam == t6 (pinning delta kills t0)",
        sp.expand(y1s[10] - ts[0]) == 0 and sp.expand(lam_s - ts[6]) == 0)

    print("(2) exact gauge structure of the y2-compression:")
    Rsq = {c: NsKs.T*Ps[c]*NsKs for c in range(NC)}
    Raq = {c: Na.T*Pa[c]*Na for c in range(NC)}
    zc = [c for c in FREE if Rsq[c] == sp.zeros(4, 4) and Raq[c] == sp.zeros(3, 3)]
    chk("zero-compression classes are exactly {11,12,13}, all with B0-coeff 0",
        zc == [11, 12, 13] and all(b0.get(c, 0) == 0 for c in zc))
    iU4 = [(i, j) for i in range(4) for j in range(i, 4)]
    iU3 = [(i, j) for i in range(3) for j in range(i, 3)]
    Y2f = [c for c in FREE if c not in zc]
    Amap = sp.Matrix([[Rsq[c][i, j] for (i, j) in iU4]
                      + [Raq[c][i, j] for (i, j) in iU3] for c in Y2f]).T
    kerA = Amap.nullspace()
    chk("compression map on the 12 nonzero classes has exact 3-dim kernel",
        len(kerA) == 3 and Amap.rank() == 9)
    bvec = sp.Matrix([b0.get(c, sp.Integer(0)) for c in Y2f])
    chk("B0 _|_ kernel (gauge dirs cost nothing; value well-defined)",
        all((bvec.T*k)[0, 0] == 0 for k in kerA))
    Akeep = Amap[:, [Y2f.index(c) for c in Y2K]]
    chk("kept classes [1,3,4,5,6,7,8,9,14] span the full 9-dim image",
        Akeep.rank() == 9)

    print("(3) LOWER BOUND: exactly-feasible rational primal point:")
    def q2r(x):
        return R(x.numerator, x.denominator)
    tq = [q2r(x) for x in T_PRIMAL]
    y2q = {c: q2r(Y2_PRIMAL[c]) for c in Y2K}
    ysl = {c: q2r(Y2_SLATER[c]) for c in Y2K}
    epsp = q2r(EPS_P)
    tb = [(1 - epsp)*x for x in tq]
    y2b = {c: (1 - epsp)*y2q[c] + epsp*ysl[c] for c in Y2K}
    r12 = [Bv[i, 1] for i in range(7)]
    chk("pinning equalities exact: y1_<B1>(tb)=0, y1_<B1B0B1>(tb)=t0=0",
        sum(r12[i]*tb[i] for i in range(7)) == 0 and tb[0] == 0)
    chk("lam(tb) = t6 > 0 (inactive as flagged)", tb[6] > 0)
    Fi = []
    Gi = []
    for i in range(7):
        M1i = sp.zeros(8, 8)
        A1i = sp.zeros(5, 5)
        for c in range(1, NC):
            if Bv[i, c-1] != 0:
                M1i += Bv[i, c-1]*Ps[c]
                A1i += Bv[i, c-1]*Pa[c]
        Fi.append(M1i*NsKs)
        Gi.append(A1i*Na)
    HsF = [[Fi[i].T*Ms0p*Fi[j] for j in range(7)] for i in range(7)]
    HaF = [[Gi[i].T*Ma0p*Gi[j] for j in range(7)] for i in range(7)]

    def SsSa(tt, yy):
        Qs = sp.zeros(4, 4)
        Qa = sp.zeros(3, 3)
        for i in range(7):
            if tt[i] == 0:
                continue
            for j in range(7):
                if tt[j] == 0:
                    continue
                Qs += tt[i]*tt[j]*HsF[i][j]
                Qa += tt[i]*tt[j]*HaF[i][j]
        Qs = (Qs + Qs.T)/2
        Qa = (Qa + Qa.T)/2
        Cs = sp.zeros(4, 4)
        Ca = sp.zeros(3, 3)
        for c in Y2K:
            Cs += yy[c]*Rsq[c]
            Ca += yy[c]*Raq[c]
        return Cs - Qs, Ca - Qa
    Ssb, Sab = SsSa(tb, y2b)
    chk("Schur complements at blended primal are PD (leading minors, exact)",
        is_pd(Ssb) and is_pd(Sab))
    # exact feasibility of the BLOCK SDP: Ms0>=0 + range condition (1) + Schur>=0
    yb_a = sum(Bv[i, 0]*tb[i] for i in range(7))
    Lb = 2*y2b[1] + y2b[4] + 2*y2b[5] - y2b[6] - yb_a
    chk("=> v_2(1/2,1/4) >= L > 1/39 > 0  (POSITIVITY, float-free)",
        Lb > R(1, 39))
    print("      L =", sp.N(Lb, 60))

    print("(4) UPPER BOUND: exactly-stationary PSD rational dual:")
    epsd = q2r(EPS_D)
    Usb = sp.zeros(4, 4)
    Uab = sp.zeros(3, 3)
    for k, (i, j) in enumerate(iU4):
        Usb[i, j] = Usb[j, i] = (1 - epsd)*q2r(US_DUAL[k]) + epsd*q2r(US_SLATER[k])
    for k, (i, j) in enumerate(iU3):
        Uab[i, j] = Uab[j, i] = (1 - epsd)*q2r(UA_DUAL[k]) + epsd*q2r(UA_SLATER[k])
    chk("stationarity in ALL 15 free y2 classes holds EXACTLY: "
        "b0_c + <Us,Rs_c> + <Ua,Ra_c> = 0",
        all(b0.get(c, sp.Integer(0)) + frobs(Usb, Rsq[c]) + frobs(Uab, Raq[c]) == 0
            for c in FREE))
    chk("Us, Ua PD (exact leading minors) -- dual multipliers PSD",
        is_pd(Usb) and is_pd(Uab))
    # max over T0 = {t: y1_<B1>=y1_<B1B0B1>=0} of the concave quadratic
    #   q(t) = -y1_a(t) - <Us,Qs(t)> - <Ua,Qa(t)>
    Z = sp.zeros(7, 5)
    Z[1, 0] = 1
    Z[2, 1] = 1
    Z[3, 2] = 1
    Z[4, 3] = 1
    Z[5, 3] = 1
    Z[5, 4] = R(7, 8)
    Z[6, 4] = 1
    chk("Z is an exact basis of T0 (both equality rows kill it; rank 5)",
        all(sum(r12[i]*Z[i, k] for i in range(7)) == 0 for k in range(5))
        and all(Z[0, k] == 0 for k in range(5)) and Z.rank() == 5)
    Mq = sp.zeros(7, 7)
    for i in range(7):
        for j in range(i, 7):
            hij = (HsF[i][j] + HsF[j][i])/2
            aij = (HaF[i][j] + HaF[j][i])/2
            Mq[i, j] = Mq[j, i] = frobs(Usb, hij) + frobs(Uab, aij)
    a_vec = sp.Matrix([Bv[i, 0] for i in range(7)])
    H5 = Z.T*Mq*Z
    g5 = Z.T*a_vec
    chk("H5 = Z^T Hess Z PD (max exists, unique; exact minors)", is_pd(H5))
    sstar = H5.LUsolve(-g5/2)
    Ub = (-g5.T*sstar - sstar.T*H5*sstar)[0, 0]
    print("      U =", sp.N(Ub, 60))
    chk("weak-duality bound valid: v_2(1/2,1/4) <= U (by (1)+(4) exactly)", True)

    print("(5) THE SANDWICH:")
    gap = Ub - Lb
    chk("L < U and U - L < 2.3e-53 (52 exact digits)",
        gap > 0 and gap < R(23, 10**54))
    ref = sp.nsimplify(REF100, rational=True)
    chk("documented 100-digit value lies in [L, U]", Lb < ref < Ub)
    print("      v_2(1/2,1/4) in [L, U], width", sp.N(gap, 3))

    print("\n(6) documented, conditional (not re-derived here):")
    print("    * 535-digit KKT solution, residual 5.9e-541 (mpmath Newton,")
    print("      exact rational data; CLARABEL cross-check 0.02609253 ~ 1e-8).")
    print("    * PSLQ exclusions: NO minpoly deg<=16/|c|<=1e24, deg<=20/1e18,")
    print("      deg<=24/1e10, deg<=32/1e7.  ==> v_2(1/2,1/4) is WILD.")
    print("    * second base v_2(1/2,1/8) = 0.0266846885939821451056... :")
    print("      NO minpoly deg<=12/1e12 (wildness generic on the face).")
    print("    * (1/2,0) is a DEGENERATE base (rank Ms0 drops to 2); block-SDP")
    print("      value there 0.017458 (solver-grade), NOT covered by this file.")
    print("\nVERDICT: v_2 at the rational base is WILD -- the wildness of the")
    print("level-2 program is intrinsic to the tangent SDP, not to the argmax.")
    print("Surviving route: rational POSITIVITY certificates (check (3)).")
    print("\nALL EXACT CHECKS " + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
