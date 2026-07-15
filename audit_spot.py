"""ADVERSARIAL AUDIT spot checks (independent referee, 2026-07-14).
Paper: npa_nonexactness_proof.tex (theorem numbers below are the paper's;
development-record numbers in parentheses).

Independent of the authors' verifier code paths wherever possible --
honestly per item: S1 and S2 are freshly written code (own words/
reduction/kernel/u/y* implementations, no shared imports); S3 and S4 are
independent DRIVERS that import the development stack (scenario_ext,
LevelData/load_level of verify_small_lambda); S5 is an end-to-end numeric
SDP cross-check:
  [S1] own-code witness check: fresh implementation of words/reduction/kernel/
       u-closed-form/y*; full-space N^T Gamma(y*) N == u u^T at k = 2..6 + pins.
  [S2] long-length adversarial E(K1,K2) == 0 test: random + targeted element
       pairs with side lengths FAR beyond the exhaustive <= 20/24 grid
       (up to 61), attacking the regime-completeness of paper Lemmas
       4.8-4.9 (development record: Lemmas 8.5/8.6).
  [S3] glue check: y* actually satisfies the hypothesis of paper Theorem
       5.4 in the AUTHORS' cert block bases (Ns, Na, u, Ks) at k = 2, 3, 4
       -- the composed chain Thm 4.4 -> Thm 5.4 (development record:
       8.1 -> 7.3) was never machine-checked with y* itself
       (verify_small_lambda uses cert witnesses).  Rebuild t0 and the exact
       Schur margins with the y*-witness.
  [S4] full-matrix end-to-end: with the y*-witness jet, EXACT PD of the FULL
       unsymmetrized Gamma(y(s)) at s = 1/1000 (k = 2, 3) -- bypasses the
       swap-congruence bookkeeping entirely -- and the exact objective cubic
       recomputed from the raw functional.
  [S5] numeric SDP c_k(s): level 1+AB, 2, 3 at several s; check
       c_k(s) >= 4 - s + g_k s^2 (their bound), c_3 <= c_2 <= c_{1+AB},
       and a quadratic-fit estimate of a_2 vs the certified 0.0280587.
Exit 0 iff all gates pass.
"""
import sys
import time
import random

HERE = '/Users/antonp/research/why-tsirelson'
sys.path.insert(0, HERE)

OK = True


def chk(name, cond):
    global OK
    OK = OK and bool(cond)
    print(f"  [{'OK' if cond else 'FAIL'}] {name}")


# ===================== independent core (written fresh) =====================
def red(w):
    """reduce word over involutions {0,1}: cancel adjacent equal letters."""
    out = []
    for x in w:
        if out and out[-1] == x:
            out.pop()
        else:
            out.append(x)
    return tuple(out)


def altw(maxlen):
    ws = [()]
    for L in range(1, maxlen + 1):
        for st in (0, 1):
            ws.append(tuple((st + i) % 2 for i in range(L)))
    return ws


def alt_from(first, L):
    """the alternating word of length L starting with letter `first`."""
    return tuple((first + i) % 2 for i in range(L))


def strip0(w):
    return w[:-1] if w and w[-1] == 0 else w


def feat(w):
    i = sum(1 for x in w if x == 1)
    ta = (w[0] + w[-1] - 1) if w else 0
    sa = (w[0] - w[-1]) if w else 0
    d0 = 1 if w == (0,) else 0
    return i, ta, sa, d0


def ystar4(alpha, beta):
    """4*y* (exact integer; paper eq. (2) in sec 4.2; note sec 8.1)."""
    i, ta, sa, da = feat(alpha)
    j, tb, sb, db = feat(beta)
    return (4 * i * i * j * j - 8 * i * j * (ta * tb + sa * sb)
            + 2 * (ta + tb) - 2 * (ta * ta + tb * tb)
            + 4 * ta * ta * tb * tb + 2 * (da + db))


def mom(p, q):
    return (red(tuple(reversed(p[0])) + q[0]),
            red(tuple(reversed(p[1])) + q[1]))


def sgn(w):
    return -1 if (w and w[0] == 1) else 1


def ends0(w):
    return bool(w) and w[-1] == 0


def endw0(L):
    """alternating word of length L ending in 0 (empty if L = 0)."""
    return tuple(((L - 1 - i) % 2) for i in range(L))


def T_elem(aw, bw):
    """trivial-reduction kernel element for a reducible word pair; u value."""
    r = (strip0(aw), strip0(bw))
    u = sgn(aw) * sgn(bw) if ends0(aw) != ends0(bw) else 0
    return [(1, (aw, bw)), (-1, r)], u


def D_elem(la, lb):
    """dressing W = (endw0(la), endw0(lb)); support = 4 corners; u value."""
    wa, wb = endw0(la), endw0(lb)
    sup = [(1, (strip0(wa), strip0(wb))),
           (-1, (wa + (1,), strip0(wb))),
           (-1, (strip0(wa), wb + (1,))),
           (1, (wa + (1,), wb + (1,)))]
    u = 0 if la + lb == 0 else 2 * (-1) ** (la + lb)
    return sup, u


def E_pair(K1, K2):
    """E(K1,K2) * 4 = sum c_p c_q 4*y*(m(p,q)) - 4 u1 u2 (integer)."""
    s1, u1 = K1
    s2, u2 = K2
    tot = 0
    for c1, p in s1:
        for c2, q in s2:
            tot += c1 * c2 * ystar4(*mom(p, q))
    return tot - 4 * u1 * u2


def level_kernel(k):
    """canonical level-k kernel: T reductions of basis words + dressings."""
    out = []
    for aw in altw(k):
        for bw in altw(k - len(aw)):
            if ends0(aw) or ends0(bw):
                out.append(T_elem(aw, bw))
    for la in range(0, k - 1):
        for lb in range(0, k - 1 - la):
            out.append(D_elem(la, lb))
    return out


def S1():
    print("[S1] independent full-space witness check, k = 2..6")
    for k in range(2, 7):
        ker = level_kernel(k)
        n = len(ker)
        assert n == 2 * k * k, (k, n)
        bad = 0
        for a in range(n):
            for b in range(a, n):
                if E_pair(ker[a], ker[b]) != 0:
                    bad += 1
        chk(f"k={k}: N^T Gamma(y*) N == u u^T on all {n*(n+1)//2} pairs "
            f"(own code)", bad == 0)
    # pins (own evaluation): identity, <A1>, face coords <(B1B0)^(j-1)B1>
    pins = ystar4((), ()) == 0 and ystar4((1,), ()) == 0
    for j in range(1, 9):
        pins = pins and ystar4((), alt_from(1, 2 * j - 1)) == 0
    chk("pins: y*[1] = y*[<A1>] = y*[face coords j<=8] = 0", pins)
    # convention anchors (paper Lemma 4.2 / Thm 4.4 consistency checks;
    # note sec 8.2)
    chk("y*[<A0>] == -1/2 and B0.y* == 0 (paper Lemma 4.2 anchors)",
        ystar4((0,), ()) == -2
        and 2 * ystar4((0,), ()) + ystar4((0,), (0,))
        + 2 * ystar4((0,), (1,)) - ystar4((1,), (1,)) == 0)


# ===================== S2: long-length adversarial sampling ================
def rand_word(rng, Lmax, allow_empty=True):
    L = rng.randint(0 if allow_empty else 1, Lmax)
    if L == 0:
        return ()
    return alt_from(rng.randint(0, 1), L)


def rand_T(rng, Lmax):
    while True:
        aw, bw = rand_word(rng, Lmax), rand_word(rng, Lmax)
        if ends0(aw) or ends0(bw):
            return T_elem(aw, bw)


def rand_D(rng, Lmax):
    return D_elem(rng.randint(0, Lmax), rng.randint(0, Lmax))


def S2(npairs=150000, Lmax=61, seed=20260714):
    print(f"[S2] adversarial long-length E(K1,K2) == 0 "
          f"(sides <= {Lmax}, {npairs} random pairs + targeted extremes)")
    rng = random.Random(seed)
    t0 = time.time()
    bad = 0
    # targeted extreme elements: tiny words vs long ones, equal lengths,
    # diffs 0..6, empty/(0)/(1) sides, dressings with one empty side
    extremes = []
    for aw in [(), (0,), (1,), (0, 1), (1, 0), (1, 0, 1)]:
        for bw in [(0,), (1, 0), alt_from(0, 39), alt_from(1, 40),
                   alt_from(0, 58)]:
            if ends0(aw) or ends0(bw):
                extremes.append(T_elem(aw, bw))
    for la in [0, 1, 2, 3, 25, 40, 59]:
        for lb in [0, 1, 2, 26, 41, 60]:
            extremes.append(D_elem(la, lb))
    for Lbig in [37, 38, 39, 40, 55, 56]:
        for d in range(0, 7):
            for f1 in (0, 1):
                for f2 in (0, 1):
                    aw, bw = alt_from(f1, Lbig), alt_from(f2, Lbig + d)
                    if ends0(aw) or ends0(bw):
                        extremes.append(T_elem(aw, bw))
    ne = len(extremes)
    for a in range(ne):
        for b in range(a, ne):
            if E_pair(extremes[a], extremes[b]) != 0:
                bad += 1
    npe = ne * (ne + 1) // 2
    chk(f"targeted extremes: all {npe} pairs pass", bad == 0)
    bad = 0
    for _ in range(npairs):
        K1 = rand_T(rng, Lmax) if rng.random() < 0.6 else rand_D(rng, Lmax)
        K2 = rand_T(rng, Lmax) if rng.random() < 0.6 else rand_D(rng, Lmax)
        if E_pair(K1, K2) != 0:
            bad += 1
    chk(f"random sampling: {npairs} pairs, sides <= {Lmax}: all pass "
        f"[{time.time()-t0:.0f}s]", bad == 0)


# ===================== S3 + S4: glue and full-matrix arc ====================
def S34(ks=(2, 3, 4), full_pd_ks=(2, 3)):
    import sympy as sp
    from verify_small_lambda import (LevelData, load_level, face_value,
                                     is_pd, pinv_exact)
    import numpy as np
    R = sp.Rational
    print("[S3] glue: y* as the Theorem-5.4 witness in the AUTHORS' cert "
          "block bases")
    for k in ks:
        L = LevelData(k)
        NC = L.NC
        Ns, Na, u, Ks, Bv, Z = load_level(k)
        # y* as a class vector in THEIR indexing (my formula, their classes)
        ystar_c = [R(ystar4(*L.crep[c]), 4) for c in range(NC)]
        Ms1 = sp.zeros(L.Ps[0].rows, L.Ps[0].rows)
        Ma1 = sp.zeros(L.Pa[0].rows, L.Pa[0].rows)
        for c in range(NC):
            if ystar_c[c] != 0:
                Ms1 += ystar_c[c] * L.Ps[c]
                Ma1 += ystar_c[c] * L.Pa[c]
        Rm = Ns.T * Ms1 * Ns
        # find scalar lam with Rm == lam * u u^T
        lam = None
        for i in range(u.rows):
            if u[i] != 0:
                lam = Rm[i, i] / (u[i] * u[i])
                break
        okr = (lam is not None and lam > 0 and Rm == lam * (u * u.T))
        chk(f"k={k}: Ns^T Ms(y*) Ns == lam u_cert u_cert^T with lam = {lam} > 0"
            " (basis match, sym)", okr)
        chk(f"k={k}: Na^T Ma(y*) Na == 0 (antisym block vanishes)",
            Na.T * Ma1 * Na == sp.zeros(Na.cols, Na.cols))
        # pins + objective identities for w := y*/lam
        ca = L.W((0,), ())
        fc = [L.W((1,), ())] + [L.W((), tuple((i % 2) for i in range(1, 2*j)))
                                for j in range(2, k + 1)]
        b0c = {L.W((0,), ()): sp.Integer(2), L.W((0,), (0,)): sp.Integer(1),
               L.W((0,), (1,)): sp.Integer(2), L.W((1,), (1,)): sp.Integer(-1)}
        y1w = [v / lam for v in ystar_c]
        chk(f"k={k}: w = y*/lam has pins 0, y1_<A0> == -1/8, B0.y1 == 0",
            all(y1w[c] == 0 for c in fc) and y1w[0] == 0
            and y1w[ca] == -R(1, 8)
            and sum(co * y1w[c] for c, co in b0c.items()) == 0)
        # rebuild Thm 5.4 (7.3) with the y*-witness: base, Slater, t0,
        # exact margins
        delta_base = [None, R(1, 2), R(1, 4)] + [R(1, 2)] * (k - 2)
        y0 = [face_value(*L.crep[c], delta_base) for c in range(NC)]
        Ms0 = sp.zeros(Ms1.rows, Ms1.rows)
        Ma0 = sp.zeros(Ma1.rows, Ma1.rows)
        for c in range(NC):
            if y0[c] != 0:
                Ms0 += y0[c] * L.Ps[c]
                Ma0 += y0[c] * L.Pa[c]
        _, Ms0p = pinv_exact(Ms0)
        _, Ma0p = pinv_exact(Ma0)
        dz = [None] + [sp.Integer(0)] * k
        y2s = [(sp.Integer(1) if c == 0 else sp.Integer(0))
               - face_value(*L.crep[c], dz) for c in range(NC)]
        Ms2 = sp.zeros(Ms1.rows, Ms1.rows)
        Ma2 = sp.zeros(Ma1.rows, Ma1.rows)
        for c in range(NC):
            if y2s[c] != 0:
                Ms2 += y2s[c] * L.Ps[c]
                Ma2 += y2s[c] * L.Pa[c]
        NsKs = Ns * Ks
        Ms1w = Ms1 / lam
        Ma1w = Ma1 / lam
        F = Ms1w * NsKs
        G = Ma1w * Na
        Qs = F.T * Ms0p * F
        Qa = G.T * Ma0p * G
        Cs = NsKs.T * Ms2 * NsKs
        Ca = Na.T * Ma2 * Na
        chk(f"k={k}: range identity with y*-witness: Ns^T Ms(w) NsKs == 0",
            Ns.T * Ms1w * NsKs == sp.zeros(Ns.cols, Ks.cols))
        rho = 0.0
        for Q, C in ((Qs, Cs), (Qa, Ca)):
            Qn = np.array([[float(Q[i, j]) for j in range(Q.cols)]
                           for i in range(Q.rows)])
            Cn = np.array([[float(C[i, j]) for j in range(C.cols)]
                           for i in range(C.rows)])
            Lc = np.linalg.cholesky(Cn)
            Mn = np.linalg.solve(Lc, np.linalg.solve(Lc, Qn).T)
            rho = max(rho, float(np.max(np.linalg.eigvalsh((Mn+Mn.T)/2))))
        t0r = None
        m = 1.02
        for _ in range(8):
            cand = R(int(rho * m * 10 ** 6) + 1, 10 ** 6)
            if is_pd(cand * Cs - Qs) and is_pd(cand * Ca - Qa):
                t0r = cand
                break
            m *= 1.5
        chk(f"k={k}: exact t0 = {t0r} for the y*-witness (Schur margins PD)"
            f" => v_{k} >= 1/(1024 t0) = {float(1/(1024*t0r)):.2e} > 0",
            t0r is not None)
        if k in full_pd_ks:
            # [S4] full unsymmetrized Gamma(y(s)) exact PD + raw objective
            lam0 = 1 / (64 * t0r)
            tt = t0r * lam0 ** 2
            svals = (R(1, 1000),)
            yfun = [y0[c] + 0 for c in range(NC)]
            for s in svals:
                ys = [y0[c] + s * lam0 * y1w[c] + s * s * tt * y2s[c]
                      for c in range(NC)]
                n = L.N
                Gm = sp.zeros(n, n)
                for p in range(n):
                    for q in range(n):
                        Gm[p, q] = ys[L.cls[L.s.idx[p, q]]]
                chk(f"k={k}: FULL Gamma(y(s)) exactly PD at s = {s} "
                    "(no swap-splitting used)", is_pd(Gm))
                # raw objective from unsymmetrized words
                def yw(aw, bw):
                    return ys[L.cls[L.s.words[
                        __import__('npa_general').canon(aw, bw)]]]
                Bs = (yw((0,), ()) + yw((), (0,)) + yw((0,), (0,))
                      + yw((0,), (1,)) + yw((1,), (0,)) - yw((1,), (1,))
                      - s / 2 * (yw((0,), ()) + yw((), (0,))))
                target = 4 - s + (1 / (1024 * t0r)) * s ** 2 + tt * s ** 3
                chk(f"k={k}: raw objective B_s.y(s) == 4 - s + g s^2 + "
                    f"t0 lam^2 s^3 exactly at s = {s}", Bs == target)
            _ = yfun


# ===================== S5: numeric SDP end-to-end ==========================
def S5():
    import numpy as np
    import cvxpy as cp
    from scenario_ext import ScenarioL
    from npa_general import Scenario
    print("[S5] numeric SDP c_k(s) end-to-end")

    def cmax(scen, s):
        y = cp.Variable(scen.NV)
        cons = [scen._moment(y) >> 0, scen._w(y, (), ()) == 1]
        obj = (scen._w(y, (0,), ()) + scen._w(y, (), (0,))
               + scen._w(y, (0,), (0,)) + scen._w(y, (0,), (1,))
               + scen._w(y, (1,), (0,)) - scen._w(y, (1,), (1,))
               - s / 2 * (scen._w(y, (0,), ()) + scen._w(y, (), (0,))))
        pr = cp.Problem(cp.Maximize(obj), cons)
        pr.solve(solver=cp.CLARABEL, max_iter=200000)
        return pr.value
    s1ab = Scenario(2, 2)
    s2 = ScenarioL(2)
    s3 = ScenarioL(3)
    g2, g3 = 3.092e-3, 3.291e-5   # their small-lambda floors
    svals = [0.05, 0.03, 0.02, 0.01, 0.005]
    rows = []
    okmono = True
    okfloor = True
    for s in svals:
        c1 = cmax(s1ab, s)
        c2 = cmax(s2, s)
        c3 = cmax(s3, s)
        rows.append((s, c1, c2, c3))
        okmono &= (c3 <= c2 + 1e-8 <= c1 + 2e-8)
        okfloor &= (c2 >= 4 - s + g2 * s * s - 1e-8)
        okfloor &= (c3 >= 4 - s + g3 * s * s - 1e-8)
        print(f"    s={s:.3f}: c_1AB-4+s={c1-4+s:+.3e} c_2-4+s={c2-4+s:+.3e} "
              f"c_3-4+s={c3-4+s:+.3e}")
    chk("nesting c_3 <= c_2 <= c_1AB at all tested s", okmono)
    chk(f"claimed floors hold: c_2 >= 4-s+{g2}s^2, c_3 >= 4-s+{g3}s^2",
        okfloor)
    # a_2 estimate: fit (c_2-4+s)/s^2 = a + b s + c s^2
    xs = np.array([r[0] for r in rows])
    q2 = np.array([(r[2] - 4 + r[0]) / r[0] ** 2 for r in rows])
    A = np.vstack([np.ones_like(xs), xs, xs * xs]).T
    coef = np.linalg.lstsq(A, q2, rcond=None)[0]
    a2est = coef[0]
    print(f"    a_2 fit estimate: {a2est:.5f} (certified exact: 0.0280587; "
          f"ladder bound 1/39 = {1/39:.5f})")
    chk("a_2 fit within 15% of certified 0.0280587 and above 1/39",
        abs(a2est - 0.0280587) < 0.15 * 0.0280587 and a2est > 1 / 39)
    # 1+AB basis containment in level-2 basis (sandwich direction)
    cont = all(any(bb == b2 for b2 in s2.basis) for bb in s1ab.basis)
    chk("1+AB monomial basis is a SUBSET of the level-2 basis "
        "(=> Gamma_1AB principal submatrix of Gamma_2 => c_2 <= c_1AB)", cont)


if __name__ == "__main__":
    T0 = time.time()
    which = set(sys.argv[1:]) or {"S1", "S2", "S3", "S5"}
    if "S1" in which:
        S1()
    if "S2" in which:
        S2()
    if "S3" in which:
        S34()
    if "S5" in which:
        S5()
    print(f"\n[{time.time()-T0:.0f}s] "
          + ("ALL AUDIT SPOT CHECKS PASS" if OK else "AUDIT FAILURES FOUND"))
    raise SystemExit(0 if OK else 1)
