"""Machine check of the small-lambda positivity theorem: paper
npa_nonexactness_proof.tex, Theorem 5.4 with Lemmas 5.1-5.3 and
Theorems 3.7/3.9 (development record: Theorem 7.3, pipeline_laws_note.md
sec 7; FINDINGS N+53a -> N+53).
In development-record language: R3 (lam*_k > 0) REDUCES to R1 (witness
existence) + L2 + L4, per level.

CLAIM being verified (k = 2, 3, 4, all verdicts EXACT rational):
  Let w_k be ANY pinned unit-lambda witness, i.e. a t-vector with
     lam(w) = 1,  face pins y1(w) = 0,   (=> (y1(w), 1) in V_k)
  so that  Ns^T Ms(y1(w)) Ns = u_k u_k^T,  Na^T Ma(y1(w)) Na = 0.
  Let y2* = e_{class(1)} - y0(delta=0) be the L4 uniform Slater point, and
     Cs = (NsKs)^T Ms(y2*) (NsKs),   Ca = Na^T Ma(y2*) Na      (PD, L4),
     Qs = F^T Ms0^+ F, F = Ms(y1(w)) NsKs;  Qa = G^T Ma0^+ G, G = Ma(y1(w)) Na.
  Take any rational t0 with t0*Cs - Qs > 0 and t0*Ca - Qa > 0 (exists since
  Cs, Ca PD).  Then for every lam in (0, 1/(32 t0)) the point
     y1 = lam * y1(w),   y2 = (t0 lam^2) * y2*
  is EXACTLY feasible for the level-k tangent SDP (block-PSD via generalized
  Schur: Ms0 >= 0, range condition = identity on V_k, Schur complement
  = lam^2 (t0 Cs - Qs) >= 0), and the objective is
     B0.y2 - y1_<A0> = lam/8 + (t0 lam^2) B0.y2* = lam/8 - 4 t0 lam^2,
  maximized at lam = 1/(64 t0):   v_k >= 1/(1024 t0) > 0.

  INGREDIENT CHECKS (each exact, per level):
   (a) w_k exists with small rationals (built from the cert pinned basis Z);
   (b) pins on y2 are HOMOGENEOUS (y2* vanishes at identity class + all k
       face-coordinate classes), so t*y2* is pin-feasible for all t;
   (c) B0.y2* == -4 exactly (B0.e_1 = 0, B0.y0(delta) == 4 on the face);
   (d) the range identity holds (u^T Ks == 0 => Ns^T Ms1 NsKs == 0 on V_k;
       Na^T Ma1 Na == 0 on V_k; ker Ms0 = span Ns EXACTLY, ker Ma0 = span Na);
   (e) Cs, Ca PD exactly (the L4 compressed pair at the program's blocks);
   (f) t0 Cs - Qs and t0 Ca - Qa PD exactly (t0 rational, from a numeric
       generalized-eigenvalue estimate + margin, then verified exactly);
   (g) objective identity y1_<A0>(w) == -1/8, B0.y1(w) == 0;
   (h) final exact positivity v_k >= 1/(1024 t0) > 0.

If all pass at k = 2, 3, 4: R3 follows from R1 at those levels by a uniform
mechanism whose every ingredient except R1 itself is proven for ALL k
(paper Theorem 5.4; development record: pipeline_laws_note.md sec 7).
Exit 0 iff all checks pass.
"""
import sys
import time
import sympy as sp

sys.path.insert(0, '/Users/antonp/research/why-tsirelson')
from scenario_ext import ScenarioL          # noqa: E402
from npa_general import canon               # noqa: E402

R = sp.Rational
OK = True


def chk(name, cond):
    global OK
    OK = OK and bool(cond)
    print(f"  [{'OK' if cond else 'FAIL'}] {name}")


def is_pd(M):
    """exact PD test: symmetric Gaussian elimination, all pivots > 0."""
    n = M.rows
    A = [[sp.Rational(M[i, j]) for j in range(n)] for i in range(n)]
    for k in range(n):
        piv = A[k][k]
        if piv <= 0:
            return False
        for i in range(k + 1, n):
            f = A[i][k] / piv
            if f == 0:
                continue
            for j in range(k, n):
                A[i][j] -= f * A[k][j]
    return True


def pinv_exact(M):
    cols = []
    for j in range(M.cols):
        if M[:, cols + [j]].rank() == len(cols) + 1:
            cols.append(j)
    Rc = M[:, cols]
    return Rc, Rc * (Rc.T * M * Rc).inv() * Rc.T


def face_value(a, b, delta):
    """product law: m_00 = 1, m_0j = delta_j, m_ij = delta_i + delta_j - 1."""
    i = sum(1 for x in a if x == 1)
    j = sum(1 for x in b if x == 1)
    if i == 0 and j == 0:
        return sp.Integer(1)
    if i == 0:
        return delta[j]
    if j == 0:
        return delta[i]
    return delta[i] + delta[j] - 1


class LevelData:
    def __init__(self, k):
        self.k = k
        s = ScenarioL(k)
        self.s = s
        self.N = s.N
        inv = {v: kk for kk, v in s.words.items()}
        NV = s.NV
        swp = [s.words[canon(bw, aw)] for i in range(NV)
               for aw, bw in [inv[i]]]
        cls = [-1] * NV
        nc = 0
        for i in range(NV):
            if cls[i] < 0:
                cls[i] = cls[swp[i]] = nc
                nc += 1
        self.cls, self.NC = cls, nc
        self.crep = {}
        for i in range(NV):
            self.crep.setdefault(cls[i], inv[i])
        bi = {b: i for i, b in enumerate(s.basis)}
        self.bi = bi
        sig = [bi[(bw, aw)] for (aw, bw) in s.basis]
        fixed = [i for i in range(self.N) if sig[i] == i]
        pairs = sorted({tuple(sorted((i, sig[i]))) for i in range(self.N)
                        if sig[i] != i})
        E = sp.eye(self.N)
        self.Tm = sp.Matrix.hstack(*([E[:, i] for i in fixed]
                                     + [E[:, i] + E[:, j] for i, j in pairs]))
        self.Um = sp.Matrix.hstack(*[E[:, i] - E[:, j] for i, j in pairs])
        Pf = [sp.zeros(self.N, self.N) for _ in range(self.NC)]
        for p in range(self.N):
            for q in range(self.N):
                Pf[cls[s.idx[p, q]]][p, q] += 1
        self.Ps = [self.Tm.T * P * self.Tm for P in Pf]
        self.Pa = [self.Um.T * P * self.Um for P in Pf]

    def W(self, a, b):
        return self.cls[self.s.words[canon(tuple(a), tuple(b))]]


def level2_data():
    """inline exact data from verify_v2_rational_base.py."""
    Ns = sp.Matrix([[-2, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0, 1, 0],
                    [-1, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0, 0, 1],
                    [1, 0, 1, 0, -1, 0, 0, 0]]).T
    Na = sp.Matrix([[1, 0, 0, 0, 0], [0, -1, 0, 1, 0], [0, -1, 0, 0, 1]]).T
    u = sp.Matrix([1, -1, 0, -1, 0])
    Ks = sp.Matrix([[1, 1, 0, 0, 0], [0, 0, 1, 0, 0], [1, 0, 0, 1, 0],
                    [0, 0, 0, 0, 1]]).T
    Bv = sp.Matrix([
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, -1, 0, 0, -1, -2, 1, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        [0, 1, R(1, 2), 0, 1, 2, 0, R(1, 2), R(3, 2), 0, 0, 0, 0, 1, 1, 0,
         1, 0],
        [-R(1, 8), -R(7, 8), -R(3, 8), 0, -R(3, 4), -R(7, 4), -R(1, 8),
         -R(1, 2), -1, 0, 0, 0, 0, 0, -1, 0, 0, 1]])
    Z = sp.zeros(7, 5)
    Z[1, 0] = 1
    Z[2, 1] = 1
    Z[3, 2] = 1
    Z[4, 3] = 1
    Z[5, 3] = 1
    Z[5, 4] = R(7, 8)
    Z[6, 4] = 1
    return Ns, Na, u, Ks, Bv, Z


def load_level(k):
    if k == 2:
        return level2_data()
    if k == 3:
        from v3_cert_data import NS3, NA3, BV3, U3, KS3, Z3
        q = (sp.Matrix(NS3), sp.Matrix(NA3), sp.Matrix(U3), sp.Matrix(KS3),
             sp.Matrix([[R(x.numerator, x.denominator)
                         if hasattr(x, 'numerator') else x for x in row]
                        for row in BV3]), sp.Matrix(Z3))
        return q[0], q[1], q[2], q[3], q[4], q[5]
    from v4_cert_data import NS4, NA4, BV4, U4, KS4, Z4
    return (sp.Matrix(NS4), sp.Matrix(NA4), sp.Matrix(U4), sp.Matrix(KS4),
            sp.Matrix([[R(x.numerator, x.denominator)
                        if hasattr(x, 'numerator') else x for x in row]
                       for row in BV4]), sp.Matrix(Z4))


def run_level(k):
    print(f"\n===== level k = {k} =====")
    t0c = time.time()
    L = LevelData(k)
    NC = L.NC
    Ns, Na, u, Ks, Bv, Z = load_level(k)
    nT = Bv.rows                       # dim of the t-space
    lamcol = NC - 1                    # Bv columns: classes 1..NC-1, then lam
    b0c = {L.W((0,), ()): sp.Integer(2), L.W((0,), (0,)): sp.Integer(1),
           L.W((0,), (1,)): sp.Integer(2), L.W((1,), (1,)): sp.Integer(-1)}
    ca = L.W((0,), ())
    # face-coordinate classes: beta = <A1>, delta_j = <(B1B0)^{j-1}B1>
    fc = [L.W((1,), ())] + [L.W((), tuple((i % 2) for i in range(1, 2 * j)))
                            for j in range(2, k + 1)]
    delta_base = [None, R(1, 2), R(1, 4)] + [R(1, 2)] * (k - 2)
    y0 = [face_value(*L.crep[c], delta_base) for c in range(NC)]
    Ms0 = sp.zeros(L.Ps[0].rows, L.Ps[0].rows)
    Ma0 = sp.zeros(L.Pa[0].rows, L.Pa[0].rows)
    for c in range(NC):
        if y0[c] != 0:
            Ms0 += y0[c] * L.Ps[c]
            Ma0 += y0[c] * L.Pa[c]
    ns, na = Ms0.rows, Ma0.rows

    # ---- (d1) base kernels exact: ker Ms0 = span Ns, ker Ma0 = span Na ----
    chk("base: Ms0 Ns == 0, Ma0 Na == 0, ranks exact (kernels = spans)",
        Ms0 * Ns == sp.zeros(ns, Ns.cols) and Ma0 * Na == sp.zeros(na, Na.cols)
        and Ms0.rank() == ns - Ns.cols and Ma0.rank() == na - Na.cols
        and Ns.rank() == Ns.cols and Na.rank() == Na.cols)
    RcS, Ms0p = pinv_exact(Ms0)
    RcA, Ma0p = pinv_exact(Ma0)
    chk("base: Ms0 >= 0, Ma0 >= 0 (Sylvester congruence), Moore-Penrose exact",
        is_pd(RcS.T * Ms0 * RcS) and is_pd(RcA.T * Ma0 * RcA)
        and Ms0 * Ms0p * Ms0 == Ms0 and Ms0p.T == Ms0p
        and (Ms0 * Ms0p).T == Ms0 * Ms0p
        and Ma0 * Ma0p * Ma0 == Ma0 and Ma0p.T == Ma0p
        and (Ma0 * Ma0p).T == Ma0 * Ma0p)

    # ---- (a) pinned unit-lambda witness w from the cert pinned basis Z ----
    lamZ = [sum(Z[i, j] * Bv[i, lamcol] for i in range(nT))
            for j in range(Z.cols)]
    jz = next(j for j in range(Z.cols) if lamZ[j] != 0)
    wt = [Z[i, jz] / lamZ[jz] for i in range(nT)]
    y1w = [sp.Integer(0)] * NC
    for i in range(nT):
        if wt[i] != 0:
            for c in range(1, NC):
                y1w[c] += wt[i] * Bv[i, c - 1]
    lam_w = sum(wt[i] * Bv[i, lamcol] for i in range(nT))
    chk("witness w: lam(w) == 1 and ALL face pins y1(w) == 0 (pins are "
        "homogeneous linear, so lam*w is pin-feasible for all lam)",
        lam_w == 1 and all(y1w[c] == 0 for c in fc))
    Ms1 = sp.zeros(ns, ns)
    Ma1 = sp.zeros(na, na)
    for c in range(1, NC):
        if y1w[c] != 0:
            Ms1 += y1w[c] * L.Ps[c]
            Ma1 += y1w[c] * L.Pa[c]
    chk("witness w: Ns^T Ms(y1w) Ns == u u^T and Na^T Ma(y1w) Na == 0 "
        "(R1 instance, unit lambda)",
        Ns.T * Ms1 * Ns == u * u.T
        and Na.T * Ma1 * Na == sp.zeros(Na.cols, Na.cols))

    # ---- (g) objective identities at w ----
    B0y1 = sum(co * y1w[c] for c, co in b0c.items())
    chk("witness w: y1_<A0>(w) == -1/8 and B0.y1(w) == 0 "
        "(objective = B0.y2 + lam/8)",
        y1w[ca] == -R(1, 8) and B0y1 == 0)

    # ---- (d2) range identity ----
    NsKs = Ns * Ks
    chk("range identity: u^T Ks == 0, Ns^T Ms(y1w) NsKs == 0, "
        "Na^T Ma(y1w) Na == 0  => block-PSD <=> Schur-PSD (Albert)",
        u.T * Ks == sp.zeros(1, Ks.cols)
        and Ns.T * Ms1 * NsKs == sp.zeros(Ns.cols, Ks.cols))

    # ---- (b)+(c) the uniform Slater point y2* = e_1 - y0(delta = 0) ----
    dzero = [None] + [sp.Integer(0)] * k
    y2s = [(sp.Integer(1) if c == 0 else sp.Integer(0))
           - face_value(*L.crep[c], dzero) for c in range(NC)]
    chk("y2*: ALL pins vanish (identity class + k face coordinates) -- "
        "pin set is HOMOGENEOUS, t*y2* pin-feasible for all t",
        y2s[0] == 0 and all(y2s[c] == 0 for c in fc))
    B0y2s = sum(co * y2s[c] for c, co in b0c.items())
    chk("y2*: B0.y2* == -4 exactly", B0y2s == -4)

    # ---- (e) compressed Slater pair PD at the program's blocks ----
    Cs = sp.zeros(Ks.cols, Ks.cols)
    Ca = sp.zeros(Na.cols, Na.cols)
    for c in range(1, NC):
        if y2s[c] != 0:
            Cs += NsKs.T * L.Ps[c] * NsKs * y2s[c]
            Ca += Na.T * L.Pa[c] * Na * y2s[c]
    chk("y2*: Cs = (NsKs)^T Ms(y2*) NsKs and Ca = Na^T Ma(y2*) Na are PD "
        "exactly (L4 at the program's blocks)", is_pd(Cs) and is_pd(Ca))

    # ---- (f) transfer ratio t0: numeric estimate, exact verification ----
    F = Ms1 * NsKs
    G = Ma1 * Na
    Qs = F.T * Ms0p * F
    Qa = G.T * Ma0p * G
    Qs = (Qs + Qs.T) / 2
    Qa = (Qa + Qa.T) / 2
    import numpy as np
    rho = 0.0
    for Q, C in ((Qs, Cs), (Qa, Ca)):
        Qn = np.array([[float(Q[i, j]) for j in range(Q.cols)]
                       for i in range(Q.rows)])
        Cn = np.array([[float(C[i, j]) for j in range(C.cols)]
                       for i in range(C.rows)])
        Lc = np.linalg.cholesky(Cn)
        Mn = np.linalg.solve(Lc, np.linalg.solve(Lc, Qn).T)
        rho = max(rho, float(np.max(np.linalg.eigvalsh((Mn + Mn.T) / 2))))
    t0 = None
    m = 1.02
    for _ in range(6):
        cand = R(int(rho * m * 10**6) + 1, 10**6)
        if is_pd(cand * Cs - Qs) and is_pd(cand * Ca - Qa):
            t0 = cand
            break
        m *= 1.5
    chk(f"exact rational t0 = {t0} with t0*Cs - Qs > 0 AND t0*Ca - Qa > 0 "
        "(transfer ratio certified)", t0 is not None)

    # ---- (h) the explicit feasible point and its exact objective ----
    lam = 1 / (64 * t0)
    tt = t0 * lam**2
    # Schur complements of the two blocks at (y1, y2) = (lam*y1w, tt*y2*):
    #   tt*Cs - lam^2*Qs = lam^2 (t0 Cs - Qs) > 0, same antisym -- exact by
    #   construction; re-verify literally:
    chk("feasible point: Schur complements tt*Cs - lam^2*Qs > 0 and "
        "tt*Ca - lam^2*Qa > 0 EXACTLY (with Ms0, Ma0 >= 0 + range identity "
        "=> both block matrices PSD, Albert's lemma)",
        is_pd(tt * Cs - lam**2 * Qs) and is_pd(tt * Ca - lam**2 * Qa))
    obj = sum(co * (tt * y2s[c]) for c, co in b0c.items()) - lam * y1w[ca]
    chk("objective == lam/8 - 4*t0*lam^2 == 1/(1024*t0) EXACTLY",
        obj == lam / 8 - 4 * t0 * lam**2 and obj == 1 / (1024 * t0))
    chk(f"=> v_{k} >= 1/(1024*t0) = {1/(1024*t0)} = "
        f"{float(1/(1024*t0)):.3e} > 0  (POSITIVITY from R1 alone)",
        obj > 0)
    print(f"    [t0 = {t0} ~ {float(t0):.4f}; rho_num = {rho:.4f}; "
          f"lam = {float(lam):.3e}; {time.time()-t0c:.1f}s]")
    return float(1 / (1024 * t0))


if __name__ == "__main__":
    t0 = time.time()
    bounds = {}
    for k in (2, 3, 4):
        bounds[k] = run_level(k)
    print(f"\nsmall-lambda bounds: " + ", ".join(
        f"v_{k} >= {b:.3e}" for k, b in bounds.items()))
    print("(ladder cert bounds for comparison: v_2 > 1/39, v_3 > 1/188, "
          "v_4 > 1/641 -- the small-lambda bounds are weaker but need ONLY "
          "R1 + proven laws)")
    print(f"\n[{time.time()-t0:.1f}s] "
          + ("ALL CHECKS PASS" if OK else "CHECK FAILURES"))
    raise SystemExit(0 if OK else 1)
