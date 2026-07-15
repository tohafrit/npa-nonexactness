"""Float-free verifier: the LEVEL-4 rational-base tangent value
v_4(1/2, 1/4, 1/2, 1/2) of the ScenarioL(4) (41x41) second-order tangent SDP
for the symmetric doubly-tilted CHSH functional B_s = (1-s/2)(<A0>+<B0>)+CHSH
is POSITIVE, certified in purely rational arithmetic:

    v_4(1/2, 1/4, 1/2, 1/2)  >=  L  >  1/641  >  0      (all checks exact)

with L = 1952448434305137/1250000000000000000 = 0.00156195874744...,
CLARABEL reference value v_4 = 0.00156556 (solver-grade; L captures 99.77%).
THIRD rung of the rational-positivity ladder:
    v_2(1/2,1/4) > 1/39,   v_3(1/2,1/4,1/2) > 1/188,   v_4 > 1/641.

THE LEVEL-4 FACE (checks (0)): the s->0 optimal face of {Gamma_4(y) >= 0,
y_1 = 1, B0.y = 4} in the party-symmetrized class space (55 classes) is
EXACTLY 4-dimensional and AFFINE, coordinates (beta, delta2, delta3, delta4)
with delta_j = <(B1B0)^{j-1} B1> -- the face law of N+50 continues.  New
forced relations at level 4 include the first conjugate-x-conjugate cross
   <A1 (B1B0)^3 B1>          = beta   + delta4 - 1
   <A1A0A1 (B1B0)^2 B1>      = delta2 + delta3 - 1   (NEW shape: d_j + d_k - 1)
and B0.y == 4 IDENTICALLY on the face.  The kernel of Gamma_4 on the face is
the face-independent INTEGER 32-dim space = 26 one-sided A0/B0-reductions
+ 6 dressings of R1 = (1-A1)(1-B1)|psi>:
   {R1, A0R1, B0R1, A1A0R1, B1B0R1, A0B0R1}
-- level 4 picks up the DOUBLE dressing A0B0R1 (absent at level 3); the
dressing count follows |w| <= k-2.  Sym/antisym split 18 + 14, base ranks
5 + 4 (levels 2,3: ranks 3+2, 4+3 -- the (k+1, k) pattern).

u_4 AND THE UNLOCK-RAY STRUCTURE (checks (2), the campaign deliverable):
  In the CANONICAL structured kernel basis [embedded level-3 sym kernel (10)
  | new trivial sym-reductions (6) | dressA1A0R1+B1B0R1 | dressA0B0R1]:

     u_4 = ( u_3 | -1, 0, -1, -1, -1, 0 | 2 | 1 )      -- INTEGER,

  i.e. u_4 EXTENDS u_3 EXACTLY (old 10 components identical), with new
  components (-1,0,-1,-1,-1,0) on the new trivial-reduction vectors and
  (2, 1) on the two new dressings.  (Extension chain: u_2 = (1,-1,0,-1,0),
  u_3 = u_2 + (1,1,0,1 | -2), u_4 = u_3 + new part above.)

  NEW STRUCTURAL DISCOVERY (exact, verified here): the rank-1 unlock
  direction is NOT unique.  The PSD first-order unlock forms
  {N^T M(y1) N : antisym form = 0} contain a PROJECTIVE LINE of rank-1
  rays u_4(t) = u_4 + t*r_4,  r_4 = e_R1 - 2 e_dress3 + 2 e_dress4s
  + e_dress4ab (= e_R1 + dressing-part of u_4); u_4 is the t = 0 member.
  The SAME structure exists at level 3 (u_3(t) with r_3 = e_R1 - 2 e_dress3,
  u_3 = t=0) -- found retroactively; the previously reported "the" u_3 is
  the canonical member of a line of rays.  VERIFIED EXACTLY here:
  vec(u4(t) u4(t)^T) lies in the exact image space S for all t (checked at
  t = 0, 1, 2 -- 3 points suffice for a quadratic curve in a subspace),
  while the PADDED u_3 (zero new components) is NOT in S: it admits NO
  first-order unlock at level 4 (its V-space forces lam == 0).  The naive
  "pad with zeros" extension is dead; the dressing components are FORCED.

LEVEL-UNIFORM STRUCTURE VERIFIED EXACTLY HERE (same shape as levels 2, 3):
  * V4 = solution space of {Ns^T M(y1) Ns == lam u4 u4^T, Na^T M(y1) Na == 0}
    has dim EXACTLY 11  (levels 2,3: 7, 9 -- the 2k+3 pattern).
  * B0.y1 == 0 IDENTICALLY on V4; y1_<A0> == -lam/8 IDENTICALLY, hence
    objective = B0.y2 + lam/8 and at the optimum B0.y2 = -lam/16:
    v_k = lam*_k/16 at ALL THREE levels
    (0.41748/16 = 0.02609; 0.08546/16 = 0.00534; 0.02509/16 = 0.00157).
  * The 4 face pins leave a 7-dim t-space (Z4) on which lam survives.
  * y2-compression gauge: zero classes exactly {42,43,44} (the (B1B0)^3B1-
    padding family; levels 2,3: {11,12,13}, {24,31,32}), all with
    B0-coefficient 0; the exact compression kernel (dim 5) is ORTHOGONAL
    to B0 (gauge directions cost nothing => the tangent SDP is bounded).
    CAUTION (new level-4 lesson, documented): replacing the rank-1 ray by
    the full 4-dim unlock FACE makes the compression coarser and B0 is
    then NOT orthogonal to the gauge -- that "face program" is unbounded
    (overclaims); the RAY-restricted program used here is the sound one.

CERTIFICATE METHOD (identical to v2/v3): rationalize the CLARABEL-optimal
jet (denominator 1e12), blend eps_p = 1e-7 toward an exact rational Slater
point of the y2-compression (strictly PD pair, margin ~9.5), verify
feasibility by exact rational elimination pivots (Sylvester).  Then
   L = B0.y2_blend - y1_<A0>(t_blend)  <=  v_4  <=  a_4.

EXACT vs CONDITIONAL, honest split:
  EXACT (this file, sympy rationals, no floats): all checks (0)-(4).
  CONDITIONAL (documented, not re-derived): the CLARABEL reference value
  and its optimality (L is a lower bound of the tangent value regardless);
  the tangent-program framework inequality v_4 <= a_4 (Bonnans-Shapiro
  second-order expansion, as at levels 1+AB, 2, 3).  NEW HONEST FLAG: at
  level 3 the pinned-path extrapolation of (c(s)-4+s)/s^2 matched the
  tangent value (0.0053-0.0054 vs 0.00534); at level 4 the same-method fit
  over s in [0.02, 0.08] is unstable (0.0002-0.0016 depending on fit
  order) and does not yet cleanly confirm 0.00157 -- higher-precision
  path solves are needed; the discrepancy is flagged, not resolved.

Requires sympy + scenario_ext/npa_general + v3_cert_data + v4_cert_data.
Runtime ~1-2 min.  Exit code 0 iff all exact checks pass.
"""
import sys
import time
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy import QQ

sys.path.insert(0, '/Users/antonp/research/why-tsirelson')

R = sp.Rational


def ddim(A):
    """exact nullity of a sympy Matrix via DomainMatrix (fast)."""
    dm = DomainMatrix.from_Matrix(A).convert_to(QQ)
    return A.cols - dm.rank()


def drank(A):
    return DomainMatrix.from_Matrix(A).convert_to(QQ).rank()


def main():
    t00 = time.time()
    from scenario_ext import ScenarioL
    from npa_general import canon
    from v3_cert_data import NS3
    from v4_cert_data import (NS4, NA4, BV4, U4, R4, U3PAD, KS4, Z4,
                              S_PRIMAL, Y2_PRIMAL, Y2_SLATER, EPS_P)
    ok = True

    def chk(name, cond):
        nonlocal ok
        ok = ok and bool(cond)
        print(f"  [{'OK' if cond else 'FAIL'}] {name}")

    # ---- exact scenario tables ----
    s4 = ScenarioL(4)
    N = s4.N
    inv = {v: k for k, v in s4.words.items()}
    NV = s4.NV
    swp = [s4.words[canon(bw, aw)] for i in range(NV) for aw, bw in [inv[i]]]
    cls = [-1]*NV
    ncls = 0
    for i in range(NV):
        if cls[i] < 0:
            cls[i] = cls[swp[i]] = ncls
            ncls += 1
    NC = ncls
    assert N == 41 and NC == 55
    bi = {b: i for i, b in enumerate(s4.basis)}

    def W(a, b):
        return cls[s4.words[canon(a, b)]]

    sig = [bi[(bw, aw)] for (aw, bw) in s4.basis]
    fixed = [i for i in range(N) if sig[i] == i]
    pairs = sorted({tuple(sorted((i, sig[i]))) for i in range(N) if sig[i] != i})
    E = sp.eye(N)
    Tm = sp.Matrix.hstack(*([E[:, i] for i in fixed]
                            + [E[:, i]+E[:, j] for i, j in pairs]))
    Um = sp.Matrix.hstack(*[E[:, i]-E[:, j] for i, j in pairs])
    Pf = [sp.zeros(N, N) for _ in range(NC)]
    for p in range(N):
        for q in range(N):
            Pf[cls[s4.idx[p, q]]][p, q] += 1
    Ps = [Tm.T*P*Tm for P in Pf]
    Pa = [Um.T*P*Um for P in Pf]

    # ---- (0) the level-4 face ----
    print("(0) level-4 s->0 face: 32-dim integer kernel + exact 4-dim affine hull:")

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
    for q, (aw, bw) in enumerate(s4.basis):
        r = colreduce(aw, bw)
        if r != (aw, bw):
            v = [0]*N
            v[q] = 1
            v[bi[r]] = -1
            kervecs.append(v)
    chk("26 trivial A0/B0-reduction relations", len(kervecs) == 26)

    def dressvec(wa, wb):
        v = [0]*N
        for da, db, sg in [((), (), 1), ((1,), (), -1),
                           ((), (1,), -1), ((1,), (1,), 1)]:
            a, b = colreduce(tuple(wa)+da, tuple(wb)+db)
            v[bi[(a, b)]] += sg
        return v
    DRESS = [((), ()), ((0,), ()), ((), (0,)),
             ((1, 0), ()), ((), (1, 0)), ((0,), (0,))]
    for wa, wb in DRESS:
        kervecs.append(dressvec(wa, wb))
    Nk = sp.Matrix(kervecs).T
    chk("kernel = 26 reductions + 6 dressings {R1,A0R1,B0R1,A1A0R1,B1B0R1,"
        "A0B0R1}: 41x32, rank 32", Nk.shape == (41, 32) and drank(Nk) == 32)

    be, d2, d3, d4 = sp.symbols('beta delta2 delta3 delta4')
    ysym = [sp.Symbol(f'q{c}') for c in range(NC)]
    Gs = sp.zeros(N, N)
    for p in range(N):
        for q in range(N):
            Gs[p, q] = ysym[cls[s4.idx[p, q]]]
    GN = Gs*Nk
    eqs = [GN[i, j] for i in range(N) for j in range(32) if GN[i, j] != 0]
    A_, bb_ = sp.linear_eq_to_matrix(eqs + [ysym[0] - 1], ysym)
    chk("face affine hull has dim EXACTLY 4 (rank 51 of 55)", drank(A_) == 51)
    # solve for the face parametrization YF in (beta, delta2, delta3, delta4)
    sol = sp.linsolve((A_, bb_), ysym)
    sol = list(sol)[0]
    cb = W((1,), ())
    cd2 = W((), (1, 0, 1))
    cd3 = W((), (1, 0, 1, 0, 1))
    cd4 = W((), (1, 0, 1, 0, 1, 0, 1))
    free = sorted({s for e in sol for s in e.free_symbols}, key=lambda s: s.name)
    subsol = sp.solve([sol[cb]-be, sol[cd2]-d2, sol[cd3]-d3, sol[cd4]-d4],
                      free, dict=True)
    chk("face coords = (beta, delta2, delta3, delta4), delta_j = <(B1B0)^{j-1}B1>",
        len(subsol) == 1)
    YF = [sp.expand(e.subs(subsol[0])) for e in sol]
    Gam = sp.zeros(N, N)
    for p in range(N):
        for q in range(N):
            Gam[p, q] = YF[cls[s4.idx[p, q]]]
    chk("Gamma(b,d2,d3,d4)*Nk == 0 IDENTICALLY (face parametrization)",
        sp.expand(Gam*Nk) == sp.zeros(N, 32))
    b0c = {W((0,), ()): sp.Integer(2), W((0,), (0,)): sp.Integer(1),
           W((0,), (1,)): sp.Integer(2), W((1,), (1,)): sp.Integer(-1)}
    chk("B0.y == 4 IDENTICALLY on the face",
        sp.expand(sum(co*YF[c] for c, co in b0c.items()) - 4) == 0)
    cX1 = W((1,), (1, 0, 1, 0, 1, 0, 1))       # <A1 (B1B0)^3 B1>
    cX2 = W((1, 0, 1), (1, 0, 1, 0, 1))        # <A1A0A1 (B1B0)^2 B1>
    chk("new forced relations: <A1(B1B0)^3B1> = beta+delta4-1 and "
        "<A1A0A1 (B1B0)^2B1> = delta2+delta3-1",
        sp.expand(YF[cX1] - (be + d4 - 1)) == 0
        and sp.expand(YF[cX2] - (d2 + d3 - 1)) == 0)

    # ---- (1) rational base point (1/2, 1/4, 1/2, 1/2) ----
    print("(1) rational base (1/2,1/4,1/2,1/2): kernels + exact PSD:")
    Ns4 = sp.Matrix(NS4)
    Na4 = sp.Matrix(NA4)
    Ms0s = Tm.T*Gam*Tm
    Ma0s = Um.T*Gam*Um
    chk("Ms0*NS4 == 0 and Ma0*NA4 == 0 IDENTICALLY in (b,d2,d3,d4)",
        sp.expand(Ms0s*Ns4) == sp.zeros(23, 18)
        and sp.expand(Ma0s*Na4) == sp.zeros(18, 14))
    sub = {be: R(1, 2), d2: R(1, 4), d3: R(1, 2), d4: R(1, 2)}
    Ms0 = Ms0s.subs(sub)
    Ma0 = Ma0s.subs(sub)
    chk("rank Ms0 = 5, rank Ma0 = 4; kernels EXACTLY span NS4 (18) / NA4 (14)",
        drank(Ms0) == 5 and drank(Ma0) == 4
        and drank(Ns4) == 18 and drank(Na4) == 14)
    # NS4 is the CANONICAL STRUCTURED basis: verify its first 10 columns are
    # the embedded level-3 sym kernel (NS3) -- the substrate of "u4 extends u3"
    s3 = ScenarioL(3)
    bi3 = {b: i for i, b in enumerate(s3.basis)}
    sig3 = [bi3[(bw, aw)] for (aw, bw) in s3.basis]
    fixed3 = [i for i in range(25) if sig3[i] == i]
    pairs3 = sorted({tuple(sorted((i, sig3[i]))) for i in range(25)
                     if sig3[i] != i})
    E3 = sp.eye(25)
    Tm3 = sp.Matrix.hstack(*([E3[:, i] for i in fixed3]
                             + [E3[:, i]+E3[:, j] for i, j in pairs3]))
    NS3w = Tm3*sp.Matrix(NS3)     # level-3 word space, 25 x 10
    emb = sp.zeros(N, 10)
    for i, w3 in enumerate(s3.basis):
        for j in range(10):
            if NS3w[i, j] != 0:
                emb[bi[w3], j] = NS3w[i, j]
    # column permutation NS3 -> canonical order [old2(5) | newtriv3(4) | dress3]
    PERM = [0, 1, 3, 4, 8, 2, 5, 6, 7, 9]
    chk("NS4 columns 0..9 == embedded level-3 kernel (canonical order)",
        Tm*Ns4[:, :10] == sp.Matrix.hstack(*[emb[:, j] for j in PERM]))

    def is_pd(M):
        """exact PD via elimination pivots (Sylvester)."""
        A = M.as_mutable()
        n = A.rows
        for k in range(n):
            if A[k, k] <= 0:
                return False
            for i in range(k+1, n):
                f = A[i, k]/A[k, k]
                if f != 0:
                    A[i, k:] = A[i, k:] - f*A[k, k:]
        return True

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
        is_pd(RcS.T*Ms0*RcS) and drank(sp.Matrix.hstack(Ns4, RcS)) == 23
        and is_pd(RcA.T*Ma0*RcA) and drank(sp.Matrix.hstack(Na4, RcA)) == 18)
    chk("Moore-Penrose exact: M M+ M = M, M+ sym, (M M+)^T = M M+ (both)",
        Ms0*Ms0p*Ms0 == Ms0 and Ms0p.T == Ms0p and (Ms0*Ms0p).T == Ms0*Ms0p
        and Ma0*Ma0p*Ma0 == Ma0 and Ma0p.T == Ma0p and (Ma0*Ma0p).T == Ma0*Ma0p)

    # ---- (2) V4, u4, the unlock-ray line, pinned space Z4 ----
    print("(2) V4 (dim 11), u4 = exact extension of u3, ray line, Z4:")
    u4 = sp.Matrix(U4)
    r4 = sp.Matrix(R4)
    Ks4 = sp.Matrix(KS4)
    Bv4 = sp.Matrix(BV4)          # 11 x 55: classes 1..54, then lam
    chk("u4^T Ks4 == 0 and rank Ks4 = 17 (exact basis of u4-perp)",
        u4.T*Ks4 == sp.zeros(1, 17) and drank(Ks4) == 17)
    chk("u4[:10] == u3 in canonical coords (1,-1,0,-1,0 | 1,1,0,1 | -2): "
        "u4 EXTENDS u3; new components (-1,0,-1,-1,-1,0 | 2 | 1)",
        list(u4[:10]) == [1, -1, 0, -1, 0, 1, 1, 0, 1, -2]
        and list(u4[10:]) == [-1, 0, -1, -1, -1, 0, 2, 1])
    ts = sp.symbols('t0:11')
    y1s = [sp.Integer(0)]*NC
    for i in range(11):
        for c in range(1, NC):
            y1s[c] += ts[i]*Bv4[i, c-1]
    lam_s = sum(ts[i]*Bv4[i, NC-1] for i in range(11))
    Ms1s = sp.zeros(23, 23)
    Ma1s = sp.zeros(18, 18)
    for c in range(1, NC):
        if y1s[c] != 0:
            Ms1s += y1s[c]*Ps[c]
            Ma1s += y1s[c]*Pa[c]
    chk("NS4^T M(y1) NS4 == lam(t) u4 u4^T IDENTICALLY in t",
        sp.expand(Ns4.T*Ms1s*Ns4 - lam_s*u4*u4.T) == sp.zeros(18, 18))
    chk("NA4^T M(y1) NA4 == 0 IDENTICALLY in t",
        sp.expand(Na4.T*Ma1s*Na4) == sp.zeros(14, 14))
    # V4 is EXACTLY the solution space (dim certificate) -- and lam is free
    yy = [sp.Integer(0)] + [sp.Symbol(f'z{c}') for c in range(1, NC)]
    lmg = sp.Symbol('lmg')
    Mzs = sp.zeros(23, 23)
    Mza = sp.zeros(18, 18)
    for c in range(1, NC):
        Mzs += yy[c]*Ps[c]
        Mza += yy[c]*Pa[c]
    G1 = sp.expand(Ns4.T*Mzs*Ns4 - lmg*u4*u4.T)
    G2 = sp.expand(Na4.T*Mza*Na4)
    eqs = ([G1[i, j] for i in range(18) for j in range(i, 18)]
           + [G2[i, j] for i in range(14) for j in range(i, 14)])
    A2_, _ = sp.linear_eq_to_matrix(eqs, yy[1:]+[lmg])
    chk("solution space of the 276 kernel conditions has dim EXACTLY 11 "
        "(pattern 7, 9, 11)", 55 - drank(A2_) == 11)
    chk("lam NOT forced to 0 on V4 (a basis row has lam = 8)",
        any(Bv4[i, NC-1] != 0 for i in range(11)))
    # the unlock-ray LINE and the death of the padded u3:
    # S := image space of attainable sym forms with antisym form == 0
    A3_, _ = sp.linear_eq_to_matrix(
        [G2[i, j] for i in range(14) for j in range(i, 14)], yy[1:])
    Y0 = DomainMatrix.from_Matrix(A3_).convert_to(QQ).nullspace().to_Matrix()
    iU18 = [(i, j) for i in range(18) for j in range(i, 18)]
    Simg = []
    for r_ in range(Y0.rows):
        Xf = sp.zeros(23, 23)
        for c in range(1, NC):
            if Y0[r_, c-1] != 0:
                Xf += Y0[r_, c-1]*Ps[c]
        Xc = Ns4.T*Xf*Ns4
        Simg.append([Xc[i, j] for (i, j) in iU18])
    SM = sp.Matrix(Simg)
    rkS = drank(SM)
    chk("image space S of unlock forms has dim EXACTLY 8", rkS == 8)

    def inS(w):
        X = w*w.T
        row = sp.Matrix([[X[i, j] for (i, j) in iU18]])
        return drank(sp.Matrix.vstack(SM, row)) == rkS
    chk("u4(t) = u4 + t*r4 is a LINE of rank-1 unlock rays: "
        "u4(t)u4(t)^T in S for t = 0, 1, 2 (quadratic curve => all t)",
        inS(u4) and inS(u4 + r4) and inS(u4 + 2*r4))
    chk("PADDED u3 (zero new components) admits NO unlock: u3pad u3pad^T not in S",
        not inS(sp.Matrix(U3PAD)))
    B0y1 = sp.expand(sum(co*y1s[c] for c, co in b0c.items()))
    chk("B0.y1 == 0 IDENTICALLY on V4 (first-order objective term vacuous)",
        B0y1 == 0)
    chk("y1_<A0> == -lam/8 IDENTICALLY on V4 (objective = B0.y2 + lam/8; "
        "v_k = lam*_k/16 at levels 2, 3, 4)",
        sp.expand(y1s[W((0,), ())] + lam_s/8) == 0)
    Z4m = sp.Matrix(Z4)
    pins = [sum(sp.Symbol(f't{i}')*Bv4[i, c-1] for i in range(11))
            for c in (cb, cd2, cd3, cd4)]
    Ap, _ = sp.linear_eq_to_matrix(pins, list(sp.symbols('t0:11')))
    chk("Z4 exact basis of pinned t-space {y1_b = y1_d2 = y1_d3 = y1_d4 = 0} "
        "(dim 7)", Ap*Z4m == sp.zeros(4, 7) and drank(Z4m) == 7
        and drank(Ap) == 4)

    # ---- (3) y2-compression gauge structure ----
    print("(3) exact gauge structure of the y2-compression:")
    NsK = Ns4*Ks4
    Rsq = {c: NsK.T*Ps[c]*NsK for c in range(NC)}
    Raq = {c: Na4.T*Pa[c]*Na4 for c in range(NC)}
    PINS2 = (0, cb, cd2, cd3, cd4)
    FREE = [c for c in range(NC) if c not in PINS2]
    zc = [c for c in FREE if Rsq[c] == sp.zeros(17, 17)
          and Raq[c] == sp.zeros(14, 14)]
    chk("zero-compression classes exactly {42,43,44}, all with B0-coeff 0",
        zc == [42, 43, 44] and all(b0c.get(c, 0) == 0 for c in zc))
    iU17 = [(i, j) for i in range(17) for j in range(i, 17)]
    iU14 = [(i, j) for i in range(14) for j in range(i, 14)]
    Y2f = [c for c in FREE if c not in zc]
    Amap = sp.Matrix([[Rsq[c][i, j] for (i, j) in iU17]
                      + [Raq[c][i, j] for (i, j) in iU14] for c in Y2f]).T
    kerA = DomainMatrix.from_Matrix(Amap).convert_to(QQ).nullspace().to_Matrix()
    chk("compression map on the 47 nonzero classes has rank 42, kernel dim 5",
        drank(Amap) == 42 and kerA.rows == 5)
    bvec = sp.Matrix([b0c.get(c, sp.Integer(0)) for c in Y2f])
    chk("B0 _|_ compression kernel (gauge dirs cost nothing; bounded program)",
        all(sum(kerA[i, j]*bvec[j] for j in range(len(Y2f))) == 0
            for i in range(kerA.rows)))

    # ---- (4) THE CERTIFICATE: blended rational primal => L > 1/641 ----
    print("(4) LOWER BOUND: exactly-feasible rational primal jet:")

    def q2r(x):
        return R(x.numerator, x.denominator)
    sq = [q2r(x) for x in S_PRIMAL]
    y2q = {c: q2r(v) for c, v in Y2_PRIMAL.items()}
    ysl = {c: q2r(v) for c, v in Y2_SLATER.items()}
    epsp = q2r(EPS_P)
    chk("primal/Slater y2 supported off the pinned + gauge classes",
        all(c in Y2f for c in y2q) and all(c in Y2f for c in ysl))
    sb = [(1 - epsp)*x for x in sq]
    y2b = {c: (1 - epsp)*y2q.get(c, sp.Integer(0)) + epsp*ysl.get(c, sp.Integer(0))
           for c in set(y2q) | set(ysl)}
    tb = Z4m*sp.Matrix(sb)
    lam_b = sum(tb[i]*Bv4[i, NC-1] for i in range(11))
    y1a_b = sum(tb[i]*Bv4[i, W((0,), ())-1] for i in range(11))
    chk("pins exact: y1_beta = y1_d2 = y1_d3 = y1_d4 = 0 at tb (by Z4), "
        "lam(tb) > 0",
        all(sum(tb[i]*Bv4[i, c-1] for i in range(11)) == 0
            for c in (cb, cd2, cd3, cd4)) and lam_b > 0)
    Csl = sp.zeros(17, 17)
    Cal = sp.zeros(14, 14)
    for c, v in ysl.items():
        if v != 0:
            Csl += v*Rsq[c]
            Cal += v*Raq[c]
    chk("Slater compressed pair strictly PD (exact elimination pivots)",
        is_pd(Csl) and is_pd(Cal))
    FZ = []
    GZ = []
    for k in range(7):
        tcol = Z4m[:, k]
        M1 = sp.zeros(23, 23)
        A1 = sp.zeros(18, 18)
        for c in range(1, NC):
            v = sum(tcol[i]*Bv4[i, c-1] for i in range(11))
            if v != 0:
                M1 += v*Ps[c]
                A1 += v*Pa[c]
        FZ.append(M1*NsK)
        GZ.append(A1*Na4)
    Qs = sp.zeros(17, 17)
    Qa = sp.zeros(14, 14)
    for i in range(7):
        if sb[i] == 0:
            continue
        for j in range(7):
            if sb[j] == 0:
                continue
            Qs += sb[i]*sb[j]*(FZ[i].T*Ms0p*FZ[j])
            Qa += sb[i]*sb[j]*(GZ[i].T*Ma0p*GZ[j])
    Qs = (Qs + Qs.T)/2
    Qa = (Qa + Qa.T)/2
    Cs = sp.zeros(17, 17)
    Ca = sp.zeros(14, 14)
    for c, v in y2b.items():
        if v != 0:
            Cs += v*Rsq[c]
            Ca += v*Raq[c]
    chk("Schur complements at blended primal are PD (exact pivots)",
        is_pd(Cs - Qs) and is_pd(Ca - Qa))
    Lb = sum(co*y2b.get(c, sp.Integer(0)) for c, co in b0c.items()) - y1a_b
    chk("=> v_4(1/2,1/4,1/2,1/2) >= L > 1/641 > 0  (POSITIVITY, float-free)",
        Lb > R(1, 641))
    print("      L =", sp.N(Lb, 30), " (= %s)" % Lb)
    print("      lam(tb) =", sp.N(lam_b, 12), "  [v_k = lam*/16 pattern]")

    print("\n(5) documented, conditional (not re-derived here):")
    print("    * CLARABEL reference v_4(1/2,1/4,1/2,1/2) = 0.00156556,")
    print("      lam* = 0.0250905; L captures 99.77% of it.")
    print("    * ladder ratios: v_3/v_2 = 0.205, v_4/v_3 = 0.293.")
    print("    * v_4 <= a_4: tangent-program framework (Bonnans-Shapiro), as at")
    print("      levels 1+AB, 2, 3.  FLAG: the level-4 pinned-path extrapolation")
    print("      (0.0002-0.0016, fit-order dependent) does not yet cleanly")
    print("      confirm the tangent value 0.00157 (level 3: clean match);")
    print("      higher-precision path solves needed.")
    print("    * the FACE-generalized program (4-dim unlock face instead of the")
    print("      u4 ray) is UNBOUNDED (B0 not orthogonal to its coarser gauge);")
    print("      the ray program used here is the sound one.")
    print("\nVERDICT: the rational-jet positivity ladder CLIMBS to level 4:")
    print("v_4(1/2,1/4,1/2,1/2) > 1/641 float-free.")
    print("Ladder: v_2 > 1/39, v_3 > 1/188, v_4 > 1/641.")
    print(f"\nruntime {time.time()-t00:.1f}s")
    print("\nALL EXACT CHECKS " + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
