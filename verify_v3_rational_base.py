"""Float-free verifier: the LEVEL-3 rational-base tangent value
v_3(1/2, 1/4, 1/2) of the ScenarioL(3) (25x25) second-order tangent SDP for
the symmetric doubly-tilted CHSH functional B_s = (1-s/2)(<A0>+<B0>) + CHSH
is POSITIVE, certified in purely rational arithmetic:

    v_3(1/2, 1/4, 1/2)  >=  L  >  1/188  >  0        (all checks exact)

with L = 0.005336346... (exact rational printed below), CLARABEL reference
value v_3 = 0.00534002157 (solver-grade).  This extends the level-2 rational
positivity certificate (verify_v2_rational_base.py: v_2(1/2,1/4) > 1/39) one
level up the NPA hierarchy -- the second rung of the conjectured
level-uniform ladder "v_k > eps_k > 0 for all k ==> a_k > 0 for all k".

THE LEVEL-3 FACE (new exact results, checks (0)):
  The s->0 optimal face of {Gamma_3(y) >= 0, y_1 = 1, B_0.y = 4} in the
  party-symmetrized class space (34 classes) is EXACTLY 3-dimensional and
  AFFINE, coordinates beta = <A1>, delta = <B1B0B1>, phi = <(B1B0)^2 B1>
  (phi is the new level-3 soft word), with the forced relations
     <A1B1>            = 2 beta - 1          (level-1 diagonal relation)
     <A1 B1B0B1>       = beta + delta - 1    (level-2 cross relation)
     <A1A0A1 B1B0B1>   = 2 delta - 1         (NEW level-3 diagonal relation)
     <A1 (B1B0)^2 B1>  = beta + phi - 1      (NEW level-3 cross relation)
  and B0.y == 4 IDENTICALLY on the face.  (This explains why the naive
  "mixed = 1/8" pinning of campaign N+49a was infeasible: the mixed moment
  is slaved, m = 2 delta - 1.)  Pattern: level-k face = (beta, delta_2,
  ..., delta_k), delta_j = <(B1B0)^{j-1} B1>, dimension k (level 1: 1-dim
  beta landscape; level 2: (beta, delta); level 3: (beta, delta, phi)).
  The kernel of Gamma_3 on the face is the (beta,delta,phi)-INDEPENDENT
  integer 18-dim space = 15 "trivial" one-sided A0/B0-reduction relations
  + the 3 dressings {R1, A0 R1, B0 R1} of the level-2 self-testing relation
  R1 = (1 - A1 - B1 + A1B1)|psi> = (1-A1)(1-B1)|psi> = 0.  Sym/antisym
  split 10 + 8; base ranks 4 + 3 (level 2: 5 + 3 and 3 + 2).

THE TANGENT PROGRAM (level-3 instance of the a2_selection_functional
pattern; v_3 <= a_3 since pinning the base and restricting the first-order
support only lower the tangent value):

  v_3(b,d,f) = max_{t in R^9, y2}   B0.y2 - y1_<A0>
     s.t. y1 = Bv3^T t  (exact 9-dim first-order space V3),
          lam(t) >= 0,  y1_beta = y1_delta = y1_phi = 0   (face pins),
          y2_1 = y2_beta = y2_delta = y2_phi = 0,
          [[Ms0, Ms(y1) Ns3 Ks3], [., (Ns3Ks3)^T Ms(y2) Ns3Ks3]] >= 0 (23x23)
          [[Ma0, Ma(y1) Na3   ], [., Na3^T     Ms(y2) Na3    ]] >= 0 (19x19)

  where V3 is the exact solution space (dim 9, rank check included) of
     Ns3^T M(y1) Ns3 == lam * u3 u3^T,   Na3^T M(y1) Na3 == 0
  with the INTEGER level-3 lifting direction u3 (path-flag from Richardson-
  extrapolated finite differences of the s>0 optimizer, residual ~1e-3):
     u3 = (1,-1,1,0,-1,1,0,1,0,-2)  in the Ns3 kernel coordinates.
  In basis space, Ns3*u3 is NOT the padded level-2 direction u2: it equals
  (padded u2) - 2*R1 + explicit support on the new dressed kernel vectors.

LEVEL-UNIFORM STRUCTURE VERIFIED EXACTLY HERE (same shape as level 2):
  * B0.y1 == 0 IDENTICALLY on V3 (first-order objective term vacuous).
  * y1_<A0> == -lam/8 IDENTICALLY on V3, hence objective = B0.y2 + lam/8;
    at the optimum B0.y2 = -lam/16, so v_k = lam*_k / 16 at BOTH levels
    (level 2: lam* = 0.41748, v_2 = 0.02609; level 3: lam* = 0.08546,
    v_3 = 0.00534): the gain is literally 1/16 of the affordable
    first-order unlock.
  * Range conditions Ns3^T Ms1(t) Ns3Ks3 == 0, Na3^T Ma1(t) Na3 == 0 hold
    identically => block-PSD <=> Schur-complement-PSD (generalized Schur).
  * The face pins kill t0 (phi-pin), as delta-pin killed t0 at level 2.
  * y2-compression gauge: zero classes exactly {24,31,32} (the
    (B1B0)^2B1-padding family; level 2: {11,12,13} = B1B0B1-paddings),
    all with B0-coefficient 0; compression map rank 23 with 4-dim kernel
    orthogonal to B0 (gauge directions cost nothing).

CERTIFICATE METHOD (v2_cert_data pattern): rationalize the CLARABEL-optimal
jet (denominator 1e12), blend eps_p = 1e-7 toward an exact rational Slater
point of the y2-compression (strictly PD compressed pair, margin ~12), and
verify feasibility by exact rational leading minors.  Then
   L = B0.y2_blend - y1_<A0>(t_blend)  <=  v_3  <=  a_3.

EXACT vs CONDITIONAL, honest split:
  EXACT (this file, sympy rationals, no floats): all checks (0)-(4) --
  face structure + kernels + identities + V3/Z3/Ks3 rank certificates +
  Slater PD + blended-primal feasibility => L > 1/188 > 0.
  CONDITIONAL (documented, not re-derived): the CLARABEL reference value
  0.00534002157 and its optimality (L is a lower bound regardless); the
  path flags used to FIND u3 (the certificate does not depend on how u3
  was found -- V3 is *defined* by u3 and the restriction only lowers v_3);
  the tangent-program framework v_3 <= a_3 (Bonnans-Shapiro second-order
  expansion of c_3(s) = 4 - s + a_3 s^2 + o(s^2), as at levels 1+AB, 2);
  no exact dual upper bound is attempted here (level-2 experience: the
  value is expected WILD; positivity, not the value, is the target).

Requires sympy + scenario_ext/npa_general + v3_cert_data.py.
Runtime ~2-3 min.  Exit code 0 iff all exact checks pass.
"""
import sys
import sympy as sp

sys.path.insert(0, '/Users/antonp/research/why-tsirelson')

R = sp.Rational


def main():
    from scenario_ext import ScenarioL
    from npa_general import canon
    from v3_cert_data import (NS3, NA3, BV3, U3, KS3, Z3,
                              S_PRIMAL, Y2_PRIMAL, Y2_SLATER, EPS_P)
    ok = True

    def chk(name, cond):
        nonlocal ok
        ok = ok and bool(cond)
        print(f"  [{'OK' if cond else 'FAIL'}] {name}")

    # ---- exact scenario tables ----
    s3 = ScenarioL(3)
    N = s3.N
    inv = {v: k for k, v in s3.words.items()}
    NV = s3.NV
    swp = [s3.words[canon(bw, aw)] for i in range(NV) for aw, bw in [inv[i]]]
    cls = [-1]*NV
    ncls = 0
    for i in range(NV):
        if cls[i] < 0:
            cls[i] = cls[swp[i]] = ncls
            ncls += 1
    NC = ncls
    assert N == 25 and NC == 34
    bi = {b: i for i, b in enumerate(s3.basis)}

    def W(a, b):
        return cls[s3.words[canon(a, b)]]

    # party swap on the 25-word basis -> sym (14) / antisym (11) compressions
    sig = [bi[(bw, aw)] for (aw, bw) in s3.basis]
    fixed = [i for i in range(N) if sig[i] == i]
    pairs = sorted({tuple(sorted((i, sig[i]))) for i in range(N) if sig[i] != i})
    E = sp.eye(N)
    Tm = sp.Matrix.hstack(*([E[:, i] for i in fixed]
                            + [E[:, i]+E[:, j] for i, j in pairs]))
    Um = sp.Matrix.hstack(*[E[:, i]-E[:, j] for i, j in pairs])
    Pf = [sp.zeros(N, N) for _ in range(NC)]
    for p in range(N):
        for q in range(N):
            Pf[cls[s3.idx[p, q]]][p, q] += 1
    Ps = [Tm.T*P*Tm for P in Pf]
    Pa = [Um.T*P*Um for P in Pf]

    # ---- (0) the level-3 face: 18-dim integer kernel, affine dim 3 ----
    print("(0) level-3 s->0 face: integer kernel + exact 3-dim affine hull:")

    def colreduce(aw, bw):
        aw, bw = list(aw), list(bw)
        ch = True
        while ch:
            ch = False
            if bw and bw[-1] == 0:
                bw.pop()
                ch = True
            if aw and aw[-1] == 0:
                aw.pop()
                ch = True
        return (tuple(aw), tuple(bw))
    kervecs = []
    for q, (aw, bw) in enumerate(s3.basis):
        r = colreduce(aw, bw)
        if r != (aw, bw):
            v = [0]*N
            v[q] = 1
            v[bi[r]] = -1
            kervecs.append(v)
    chk("15 trivial A0/B0-reduction relations", len(kervecs) == 15)

    def evec(*pl):
        v = [0]*N
        for coef, w in pl:
            v[bi[w]] += coef
        return v
    kervecs.append(evec((1, ((), ())), (-1, ((1,), ())), (-1, ((), (1,))),
                        (1, ((1,), (1,)))))                       # R1
    kervecs.append(evec((1, ((), ())), (-1, ((0, 1), ())), (-1, ((), (1,))),
                        (1, ((0, 1), (1,)))))                     # A0 R1
    kervecs.append(evec((1, ((), ())), (-1, ((1,), ())), (-1, ((), (0, 1))),
                        (1, ((1,), (0, 1)))))                     # B0 R1
    Nk = sp.Matrix(kervecs).T
    chk("kernel matrix 25x18 has rank 18", Nk.shape == (25, 18) and Nk.rank() == 18)

    be, de, ph = sp.symbols('beta delta phi')
    ga = 2*be - 1
    ep = de + be - 1
    rho = ph + be - 1
    m = 2*de - 1
    YF = {0: 1, 1: 1, 2: be, 3: be, 4: be, 5: de, 6: 1, 7: be, 8: be, 9: ga,
          10: ga, 11: de, 12: be, 13: ga, 14: ga, 15: ga, 16: de, 17: ep,
          18: ph, 19: de, 20: ep, 21: ep, 22: ga, 23: de, 24: ph, 25: de,
          26: ep, 27: ep, 28: ep, 29: ga, 30: ep, 31: ph, 32: rho, 33: m}
    Gam = sp.zeros(N, N)
    for p in range(N):
        for q in range(N):
            Gam[p, q] = YF[cls[s3.idx[p, q]]]
    chk("Gamma(beta,delta,phi) * Nk == 0 IDENTICALLY (face parametrization)",
        sp.expand(Gam*Nk) == sp.zeros(N, 18))
    # affine hull is EXACTLY 3-dim: unknown y (34 classes), Gamma(y)Nk=0, y_1=1
    ysym = [sp.Symbol(f'q{c}') for c in range(NC)]
    Gs = sp.zeros(N, N)
    for p in range(N):
        for q in range(N):
            Gs[p, q] = ysym[cls[s3.idx[p, q]]]
    GN = Gs*Nk
    eqs = [GN[i, j] for i in range(N) for j in range(18) if GN[i, j] != 0]
    A_, _ = sp.linear_eq_to_matrix(eqs + [ysym[0] - 1], ysym)
    chk("face affine hull has dim EXACTLY 3 (rank 31 of 34)", A_.rank() == 31)
    b0c = {W((0,), ()): sp.Integer(2), W((0,), (0,)): sp.Integer(1),
           W((0,), (1,)): sp.Integer(2), W((1,), (1,)): sp.Integer(-1)}
    B0y0 = sum(co*YF[c] for c, co in b0c.items())
    chk("B0.y == 4 IDENTICALLY on the face", sp.expand(B0y0 - 4) == 0)
    chk("new relations: <A1A0A1 B1B0B1> = 2 delta - 1 and "
        "<A1 (B1B0)^2 B1> = beta + phi - 1 (by parametrization)",
        YF[33] == m and YF[32] == rho)

    # ---- (1) rational base point (1/2, 1/4, 1/2) ----
    print("(1) rational base (1/2,1/4,1/2): kernels + exact PSD:")
    Ns3 = sp.Matrix(NS3)
    Na3 = sp.Matrix(NA3)
    Ms0s = Tm.T*Gam*Tm
    Ma0s = Um.T*Gam*Um
    chk("Ms0*Ns3 == 0 and Ma0*Na3 == 0 IDENTICALLY in (beta,delta,phi)",
        sp.expand(Ms0s*Ns3) == sp.zeros(14, 10)
        and sp.expand(Ma0s*Na3) == sp.zeros(11, 8))
    sub = {be: R(1, 2), de: R(1, 4), ph: R(1, 2)}
    Ms0 = Ms0s.subs(sub)
    Ma0 = Ma0s.subs(sub)
    chk("rank Ms0 = 4, rank Ma0 = 3; kernels EXACTLY span Ns3 (10) / Na3 (8)",
        Ms0.rank() == 4 and Ma0.rank() == 3
        and Ns3.rank() == 10 and Na3.rank() == 8)

    def is_pd(M):
        return all(M[:k, :k].det() > 0 for k in range(1, M.rows + 1))

    def pinv_exact(M):
        cols = []
        for j in range(M.cols):
            if M[:, cols + [j]].rank() == len(cols) + 1:
                cols.append(j)
        Rc = M[:, cols]
        return Rc, Rc*(Rc.T*M*Rc).inv()*Rc.T
    RcS, Ms0p = pinv_exact(Ms0)
    RcA, Ma0p = pinv_exact(Ma0)
    chk("Ms0 >= 0, Ma0 >= 0 (Sylvester on independent-column congruence)",
        is_pd(RcS.T*Ms0*RcS) and sp.Matrix.hstack(Ns3, RcS).rank() == 14
        and is_pd(RcA.T*Ma0*RcA) and sp.Matrix.hstack(Na3, RcA).rank() == 11)
    chk("Moore-Penrose exact: M M+ M = M, M+ sym, (M M+)^T = M M+ (both)",
        Ms0*Ms0p*Ms0 == Ms0 and Ms0p.T == Ms0p and (Ms0*Ms0p).T == Ms0*Ms0p
        and Ma0*Ma0p*Ma0 == Ma0 and Ma0p.T == Ma0p and (Ma0*Ma0p).T == Ma0*Ma0p)

    # ---- (2) the first-order space V3 and its exact identities ----
    print("(2) V3 (dim 9), lifting direction u3, pinned space Z3:")
    u3 = sp.Matrix(U3)
    Ks3 = sp.Matrix(KS3)
    Bv3 = sp.Matrix(BV3)          # 9 x 34: classes 1..33, then lam
    chk("u3^T Ks3 == 0 and rank Ks3 = 9 (exact basis of u3-perp)",
        u3.T*Ks3 == sp.zeros(1, 9) and Ks3.rank() == 9)
    ts = sp.symbols('t0:9')
    y1s = [sp.Integer(0)]*NC
    for i in range(9):
        for c in range(1, NC):
            y1s[c] += ts[i]*Bv3[i, c-1]
    lam_s = sum(ts[i]*Bv3[i, NC-1] for i in range(9))
    Ms1s = sp.zeros(14, 14)
    Ma1s = sp.zeros(11, 11)
    for c in range(1, NC):
        if y1s[c] != 0:
            Ms1s += y1s[c]*Ps[c]
            Ma1s += y1s[c]*Pa[c]
    chk("Ns3^T M(y1) Ns3 == lam(t) u3 u3^T IDENTICALLY in t",
        sp.expand(Ns3.T*Ms1s*Ns3 - lam_s*u3*u3.T) == sp.zeros(10, 10))
    chk("Na3^T M(y1) Na3 == 0 IDENTICALLY in t",
        sp.expand(Na3.T*Ma1s*Na3) == sp.zeros(8, 8))
    # V3 is EXACTLY the solution space (dim certificate)
    yy = [sp.Integer(0)] + [sp.Symbol(f'z{c}') for c in range(1, NC)]
    lmg = sp.Symbol('lmg')
    Mzs = sp.zeros(14, 14)
    Mza = sp.zeros(11, 11)
    for c in range(1, NC):
        Mzs += yy[c]*Ps[c]
        Mza += yy[c]*Pa[c]
    G1 = sp.expand(Ns3.T*Mzs*Ns3 - lmg*u3*u3.T)
    G2 = sp.expand(Na3.T*Mza*Na3)
    eqs = ([G1[i, j] for i in range(10) for j in range(i, 10)]
           + [G2[i, j] for i in range(8) for j in range(i, 8)])
    A2_, _ = sp.linear_eq_to_matrix(eqs, yy[1:]+[lmg])
    chk("solution space of the 91 kernel conditions has dim EXACTLY 9",
        34 - A2_.rank() == 9)
    B0y1 = sp.expand(sum(co*y1s[c] for c, co in b0c.items()))
    chk("B0.y1 == 0 IDENTICALLY on V3 (first-order term vacuous)", B0y1 == 0)
    chk("y1_<A0> == -lam/8 IDENTICALLY on V3 (objective = B0.y2 + lam/8)",
        sp.expand(y1s[W((0,), ())] + lam_s/8) == 0)
    Z3m = sp.Matrix(Z3)
    cb, cd, cf = W((1,), ()), W((), (1, 0, 1)), W((), (1, 0, 1, 0, 1))
    pins = [sum(sp.Symbol(f't{i}')*Bv3[i, c-1] for i in range(9))
            for c in (cb, cd, cf)]
    Ap, _ = sp.linear_eq_to_matrix(pins, list(sp.symbols('t0:9')))
    chk("Z3 exact basis of pinned t-space {y1_beta=y1_delta=y1_phi=0} (dim 6)",
        Ap*Z3m == sp.zeros(3, 6) and Z3m.rank() == 6 and Ap.rank() == 3)

    # ---- (3) y2-compression gauge structure ----
    print("(3) exact gauge structure of the y2-compression:")
    NsK = Ns3*Ks3
    Rsq = {c: NsK.T*Ps[c]*NsK for c in range(NC)}
    Raq = {c: Na3.T*Pa[c]*Na3 for c in range(NC)}
    PINS2 = (0, cb, cd, cf)
    FREE = [c for c in range(NC) if c not in PINS2]
    zc = [c for c in FREE if Rsq[c] == sp.zeros(9, 9)
          and Raq[c] == sp.zeros(8, 8)]
    chk("zero-compression classes exactly {24,31,32}, all with B0-coeff 0",
        zc == [24, 31, 32] and all(b0c.get(c, 0) == 0 for c in zc))
    iU9 = [(i, j) for i in range(9) for j in range(i, 9)]
    iU8 = [(i, j) for i in range(8) for j in range(i, 8)]
    Y2f = [c for c in FREE if c not in zc]
    Amap = sp.Matrix([[Rsq[c][i, j] for (i, j) in iU9]
                      + [Raq[c][i, j] for (i, j) in iU8] for c in Y2f]).T
    kerA = Amap.nullspace()
    chk("compression map on the 27 nonzero classes has rank 23, kernel dim 4",
        Amap.rank() == 23 and len(kerA) == 4)
    bvec = sp.Matrix([b0c.get(c, sp.Integer(0)) for c in Y2f])
    chk("B0 _|_ compression kernel (gauge dirs cost nothing)",
        all((bvec.T*k)[0, 0] == 0 for k in kerA))

    # ---- (4) THE CERTIFICATE: blended rational primal => L > 1/188 ----
    print("(4) LOWER BOUND: exactly-feasible rational primal jet:")

    def q2r(x):
        return R(x.numerator, x.denominator)
    sq = [q2r(x) for x in S_PRIMAL]
    y2q = {c: q2r(v) for c, v in Y2_PRIMAL.items()}
    ysl = {c: q2r(v) for c, v in Y2_SLATER.items()}
    epsp = q2r(EPS_P)
    chk("primal/Slater y2 supported off the pinned classes",
        all(c in FREE for c in y2q) and all(c in FREE for c in ysl))
    sb = [(1 - epsp)*x for x in sq]
    y2b = {c: (1 - epsp)*y2q.get(c, sp.Integer(0)) + epsp*ysl.get(c, sp.Integer(0))
           for c in set(y2q) | set(ysl)}
    tb = Z3m*sp.Matrix(sb)
    lam_b = sum(tb[i]*Bv3[i, NC-1] for i in range(9))
    y1a_b = sum(tb[i]*Bv3[i, W((0,), ())-1] for i in range(9))
    chk("pins exact: y1_beta = y1_delta = y1_phi = 0 at tb (by Z3), lam(tb) > 0",
        all(sum(tb[i]*Bv3[i, c-1] for i in range(9)) == 0
            for c in (cb, cd, cf)) and lam_b > 0)
    # Slater point strictly PD (compressed pair)
    Csl = sp.zeros(9, 9)
    Cal = sp.zeros(8, 8)
    for c, v in ysl.items():
        if v != 0:
            Csl += v*Rsq[c]
            Cal += v*Raq[c]
    chk("Slater compressed pair strictly PD (exact leading minors)",
        is_pd(Csl) and is_pd(Cal))
    # Schur complements at the blended point
    FZ = []
    GZ = []
    for k in range(6):
        tcol = Z3m[:, k]
        M1 = sp.zeros(14, 14)
        A1 = sp.zeros(11, 11)
        for c in range(1, NC):
            v = sum(tcol[i]*Bv3[i, c-1] for i in range(9))
            if v != 0:
                M1 += v*Ps[c]
                A1 += v*Pa[c]
        FZ.append(M1*NsK)
        GZ.append(A1*Na3)
    Qs = sp.zeros(9, 9)
    Qa = sp.zeros(8, 8)
    for i in range(6):
        if sb[i] == 0:
            continue
        for j in range(6):
            if sb[j] == 0:
                continue
            Qs += sb[i]*sb[j]*(FZ[i].T*Ms0p*FZ[j])
            Qa += sb[i]*sb[j]*(GZ[i].T*Ma0p*GZ[j])
    Qs = (Qs + Qs.T)/2
    Qa = (Qa + Qa.T)/2
    Cs = sp.zeros(9, 9)
    Ca = sp.zeros(8, 8)
    for c, v in y2b.items():
        if v != 0:
            Cs += v*Rsq[c]
            Ca += v*Raq[c]
    chk("Schur complements at blended primal are PD (exact leading minors)",
        is_pd(Cs - Qs) and is_pd(Ca - Qa))
    Lb = sum(co*y2b.get(c, sp.Integer(0)) for c, co in b0c.items()) - y1a_b
    chk("=> v_3(1/2,1/4,1/2) >= L > 1/188 > 0  (POSITIVITY, float-free)",
        Lb > R(1, 188))
    print("      L =", sp.N(Lb, 30), " (= %s)" % Lb)
    print("      lam(tb) =", sp.N(lam_b, 12), "  [v_k = lam*/16 pattern]")

    print("\n(5) documented, conditional (not re-derived here):")
    print("    * CLARABEL reference v_3(1/2,1/4,1/2) = 0.00534002157 (~1e-9),")
    print("      lam* = 0.08545812; L captures 99.93% of it.")
    print("    * u3 path-flag: Richardson FD of the s>0 optimizer; first-order")
    print("      compressed form rank 1 (top eig 1.664, next 0.0012), antisym")
    print("      form 0 at first order; sym kernel unlocks 1 dir at O(s),")
    print("      2 at O(s^2); antisym 1 at O(s^2) -- same shape as level 2.")
    print("    * v_3 <= a_3 (~0.0062-0.0096 solver-grade): tangent-program")
    print("      pinning inequality, same framework as levels 1+AB and 2.")
    print("\nVERDICT: the rational-jet positivity ladder CLIMBS to level 3:")
    print("v_3(1/2,1/4,1/2) > 1/188 float-free.  Ladder: v_2 > 1/39, v_3 > 1/188.")
    print("\nALL EXACT CHECKS " + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
