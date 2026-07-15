"""Machine check of the no-repair positivity lemma and the exact arc:
paper npa_nonexactness_proof.tex, Lemma 6.1 and Theorem 6.3
(development record: the uniform framework lemma, pipeline_laws_note.md
sec 9).

LEMMA (no-repair positivity). Let M(s) = M0 + s M1 + s^2 M2 (symmetric, any
fixed size) and decompose R^n = Rg (+) U (+) W orthogonally, Rg = range(M0),
ker M0 = U (+) W, dim U <= 1. Assume
  (a) M0 >= 0, ker M0 = U (+) W exactly; c := min spec(M0|Rg) > 0;
  (b) the kernel-compressed first-order form is rank-one ALIGNED WITH U:
      P_ker M1 P_ker = lamhat P_U, lamhat > 0 (if U = 0: identically 0).
      [This is exactly N^T M1 N = lam u u^T: W = N.(u-perp), U = ker (-) W;
      in particular the U-W and W-W first-order blocks VANISH EXACTLY.]
  (c) nothing on the Rg-kernel coupling X1 = P_Rg M1 P_ker (arbitrary);
  (d) strict second-order Schur margin on W:
      H_WW := P_W M2 P_W - X1W^T M0^+ X1W >= mu I_W,  mu > 0.
With b1 = |M1|, b2 = |M2| (spectral), h, kappa <= b2 + b1^2/c the U-U / U-W
entries of H := P_ker M2 P_ker - X1^T M0^+ X1, and
      eps := 2 b1^2 (b1+b2)/c^2 + 4 b1 b2 / c + 2 b2^2 / c,
      s* := min{1, c/(2(b1+b2)), lamhat/(2(h+eps)), mu/(2 eps),
                lamhat mu / (4 (kappa+eps)^2)}          (drop lamhat terms if U=0),
then  M(s) >= 0 for ALL s in (0, s*]  (strictly PD on (0, s*)).
NO o(s^2) correction z(s) is needed: z == 0.

Proof (paper, Lemma 6.1; development record: note sec 9): Schur against
A(s) = P_Rg M(s) P_Rg >= (c/2) I;
the kernel Schur complement is
   S(s) = s lamhat P_U + s^2 H - Rem(s),   |Rem(s)| <= s^3 eps,
because (b) kills the kernel-kernel first-order block off U and M(s) is an
EXACT quadratic; then 2x2 block dominance: U-entry >= s lamhat/2,
W-block >= s^2 mu/2, cross <= s^2 (kappa+eps), and PSD follows from
(s lamhat/2)(s^2 mu/2) >= s^4 (kappa+eps)^2, i.e. s <= lamhat mu/(4(kappa+eps)^2).

CONSEQUENCE (paper Theorem 6.3, with the Theorem-5.4 jet y1 = lam*y1w,
y2 = t0 lam^2 y2* (development record: Theorem 7.3),
lam = 1/(64 t0)): y(s) = y0 + s y1 + s^2 y2 has y(s)_1 = 1 and
Gamma(y(s)) >= 0 for s in (0, s*] (both swap blocks), and the objective is
the EXACT CUBIC  B_s.y(s) = 4 - s + g s^2 + t0 lam^2 s^3,  g = 1/(1024 t0),
so   c_k(s) >= 4 - s + g s^2   on (0, s*]:   a_k >= g > 0.
This DISCHARGES the o(s^2)-framework flag for strict-margin jets:
R1(k) => a_k > 0 (not just v_k > 0).

CHECKS per level k = 2, 3, 4 (gates; all EXACT rational unless marked):
  (1) jet rebuild = verify_small_lambda ingredients (kernel exact, rank-1
      u-alignment, t0-margins PD, pins, B0.y2* = -4);
  (2) EXACT: Ms(s) and Ma(s) are PD at s = 1/1000 and 1/10000
      (swap congruence [Tm|Um] block-diagonalizes Gamma: full-matrix PSD
      follows; (1+S)Gamma(1-S) = 0 since [Gamma(y),S]=0, S^2=1);
  (3) EXACT: objective cubic coefficients (4, -1, 1/(1024 t0), t0 lam^2);
  (4) numeric (report): lemma constants c, b1, b2, lamhat, mu, kappa, eps
      and s* -- confirm s* >= 1e-3 so the tested s are inside the interval
      (CERTIFIED rational bounds for these constants and a certified
      rational s*_k are produced by verify_sstar_rational.py);
      eigenvalue margins of M(s) at s = 1e-3, 1e-4 vs predicted floors
      (at k = 4 the W-floor ~ s^2 lam^2 is below float resolution: the
      exact check (2) is the arbiter, stated honestly).
Exit 0 iff all gates pass.
"""
import sys
import time
import numpy as np
import sympy as sp

sys.path.insert(0, '/Users/antonp/research/why-tsirelson')
from verify_small_lambda import (LevelData, load_level, face_value,  # noqa: E402
                                 is_pd, pinv_exact)

R = sp.Rational
OK = True


def chk(name, cond):
    global OK
    OK = OK and bool(cond)
    print(f"  [{'OK' if cond else 'FAIL'}] {name}")


def npm(M):
    return np.array([[float(M[i, j]) for j in range(M.cols)]
                     for i in range(M.rows)])


def snorm(A):
    return float(np.linalg.norm(A, 2))


def build_jet(k):
    """Rebuild the Theorem-5.4 small-lambda jet exactly (as
    verify_small_lambda; development record: Theorem 7.3)."""
    L = LevelData(k)
    NC = L.NC
    Ns, Na, u, Ks, Bv, Z = load_level(k)
    nT = Bv.rows
    lamcol = NC - 1
    b0c = {L.W((0,), ()): sp.Integer(2), L.W((0,), (0,)): sp.Integer(1),
           L.W((0,), (1,)): sp.Integer(2), L.W((1,), (1,)): sp.Integer(-1)}
    ca = L.W((0,), ())
    fc = [L.W((1,), ())] + [L.W((), tuple((i % 2) for i in range(1, 2 * j)))
                            for j in range(2, k + 1)]
    delta_base = [None, R(1, 2), R(1, 4)] + [R(1, 2)] * (k - 2)
    y0 = [face_value(*L.crep[c], delta_base) for c in range(NC)]

    def blocks(y):
        Ms = sp.zeros(L.Ps[0].rows, L.Ps[0].rows)
        Ma = sp.zeros(L.Pa[0].rows, L.Pa[0].rows)
        for c in range(NC):
            if y[c] != 0:
                Ms += y[c] * L.Ps[c]
                Ma += y[c] * L.Pa[c]
        return Ms, Ma

    Ms0, Ma0 = blocks(y0)
    ns, na = Ms0.rows, Ma0.rows
    # witness (unit lambda, pinned)
    lamZ = [sum(Z[i, j] * Bv[i, lamcol] for i in range(nT))
            for j in range(Z.cols)]
    jz = next(j for j in range(Z.cols) if lamZ[j] != 0)
    wt = [Z[i, jz] / lamZ[jz] for i in range(nT)]
    y1w = [sp.Integer(0)] * NC
    for i in range(nT):
        if wt[i] != 0:
            for c in range(1, NC):
                y1w[c] += wt[i] * Bv[i, c - 1]
    Ms1w, Ma1w = blocks(y1w)
    # Slater point
    dzero = [None] + [sp.Integer(0)] * k
    y2s = [(sp.Integer(1) if c == 0 else sp.Integer(0))
           - face_value(*L.crep[c], dzero) for c in range(NC)]
    Ms2s, Ma2s = blocks(y2s)
    # exact structural gates (jet ingredients)
    chk("base kernels exact: Ms0 Ns = 0, Ma0 Na = 0, ranks complementary",
        Ms0 * Ns == sp.zeros(ns, Ns.cols) and Ma0 * Na == sp.zeros(na, Na.cols)
        and Ms0.rank() == ns - Ns.cols and Ma0.rank() == na - Na.cols)
    chk("rank-1 u-alignment: Ns^T Ms(y1w) Ns == u u^T, Na^T Ma(y1w) Na == 0, "
        "u^T Ks == 0 (kernel U-W and W-W first-order blocks vanish EXACTLY)",
        Ns.T * Ms1w * Ns == u * u.T
        and Na.T * Ma1w * Na == sp.zeros(Na.cols, Na.cols)
        and u.T * Ks == sp.zeros(1, Ks.cols))
    chk("witness pins + objective identities: lam(w)=1, y1w_<A0>=-1/8, "
        "B0.y1w=0, y2* pins=0, B0.y2*=-4",
        sum(wt[i] * Bv[i, lamcol] for i in range(nT)) == 1
        and y1w[ca] == -R(1, 8)
        and sum(co * y1w[c] for c, co in b0c.items()) == 0
        and all(y1w[c] == 0 for c in fc) and y2s[0] == 0
        and all(y2s[c] == 0 for c in fc)
        and sum(co * y2s[c] for c, co in b0c.items()) == -4)
    # transfer ratio t0 (exact certificate, numeric seed)
    NsKs = Ns * Ks
    RcS, Ms0p = pinv_exact(Ms0)
    RcA, Ma0p = pinv_exact(Ma0)
    Cs = NsKs.T * Ms2s * NsKs
    Ca = Na.T * Ma2s * Na
    F = Ms1w * NsKs
    G = Ma1w * Na
    Qs = F.T * Ms0p * F
    Qa = G.T * Ma0p * G
    rho = 0.0
    for Q, C in ((Qs, Cs), (Qa, Ca)):
        Lc = np.linalg.cholesky(npm(C))
        Mn = np.linalg.solve(Lc, np.linalg.solve(Lc, npm(Q)).T)
        rho = max(rho, float(np.max(np.linalg.eigvalsh((Mn + Mn.T) / 2))))
    t0 = None
    m = 1.02
    for _ in range(6):
        cand = R(int(rho * m * 10**6) + 1, 10**6)
        if is_pd(cand * Cs - Qs) and is_pd(cand * Ca - Qa):
            t0 = cand
            break
        m *= 1.5
    chk(f"strict W-margins: t0 = {t0} with t0*Cs - Qs > 0, t0*Ca - Qa > 0 "
        "(EXACT PD -> mu > 0 in the lemma)", t0 is not None)
    return dict(L=L, NC=NC, Ns=Ns, Na=Na, u=u, Ks=Ks, y0=y0, y1w=y1w,
                y2s=y2s, Ms0=Ms0, Ma0=Ma0, Ms1w=Ms1w, Ma1w=Ma1w,
                Ms2s=Ms2s, Ma2s=Ma2s, t0=t0, b0c=b0c, ca=ca, Qs=Qs, Qa=Qa,
                Cs=Cs, Ca=Ca)


def lemma_constants(M0, M1, M2, Nker, Wcols):
    """Numeric lemma constants for one block. Wcols: matrix whose columns
    span W inside ker M0 (None => W = whole kernel, U empty)."""
    A = npm(M0)
    B = npm(M1)
    C2 = npm(M2)
    n = A.shape[0]
    ev, V = np.linalg.eigh(A)
    tol = 1e-9 * max(1.0, float(np.max(np.abs(ev))))
    c = float(min(e for e in ev if e > tol))
    Qk, _ = np.linalg.qr(npm(Nker))
    dker = Nker.cols
    if Wcols is not None:
        Qw, _ = np.linalg.qr(npm(Wcols))
        Pw = Qw @ Qw.T
        Ur = Qk - Pw @ Qk
        uu, ss, _ = np.linalg.svd(Ur)
        e_u = uu[:, 0:1]                      # 1-dim U
        frames = np.hstack([e_u, Qw])
    else:
        Qw = Qk
        e_u = None
        frames = Qk
    b1, b2 = snorm(B), snorm(C2)
    M0p = np.linalg.pinv(A, rcond=1e-9)
    H = frames.T @ (C2 - B @ M0p @ B) @ frames
    if e_u is not None:
        lamhat = float((e_u.T @ B @ e_u)[0, 0])
        h = abs(float(H[0, 0]))
        kappa = float(np.linalg.norm(H[0, 1:]))
        muW = float(np.min(np.linalg.eigvalsh(H[1:, 1:])))
        # exact-zero first-order kernel off-U check (numeric echo of gate)
        K1 = frames.T @ B @ frames
        offU = float(np.max(np.abs(K1 - np.diag([K1[0, 0]] + [0] * (dker - 1)))))
    else:
        lamhat, h, kappa = None, 0.0, 0.0
        muW = float(np.min(np.linalg.eigvalsh(H)))
        K1 = Qk.T @ B @ Qk
        offU = float(np.max(np.abs(K1)))
    eps = 2 * b1**2 * (b1 + b2) / c**2 + 4 * b1 * b2 / c + 2 * b2**2 / c
    cands = [1.0, c / (2 * (b1 + b2)), muW / (2 * eps)]
    if lamhat is not None:
        cands += [lamhat / (2 * (h + eps)),
                  lamhat * muW / (4 * (kappa + eps)**2)]
    return dict(c=c, b1=b1, b2=b2, lamhat=lamhat, h=h, kappa=kappa,
                mu=muW, eps=eps, sstar=min(cands), offU=offU)


def run_level(k):
    print(f"\n===== level k = {k} =====")
    tstart = time.time()
    J = build_jet(k)
    t0 = J['t0']
    lam = 1 / (64 * t0)
    tt = t0 * lam**2
    Ms1 = lam * J['Ms1w']
    Ma1 = lam * J['Ma1w']
    Ms2 = tt * J['Ms2s']
    Ma2 = tt * J['Ma2s']

    # ---- lemma constants + s* (numeric report; sym has U, antisym U=0) ----
    Ls = lemma_constants(J['Ms0'], Ms1, Ms2, J['Ns'], J['Ns'] * J['Ks'])
    La = lemma_constants(J['Ma0'], Ma1, Ma2, J['Na'], None)
    print(f"    sym : c={Ls['c']:.3e} b1={Ls['b1']:.3e} b2={Ls['b2']:.3e} "
          f"lamhat={Ls['lamhat']:.3e} mu={Ls['mu']:.3e} kappa={Ls['kappa']:.3e}")
    print(f"          eps={Ls['eps']:.3e} offU(first-order kernel off-U)="
          f"{Ls['offU']:.1e}  s*_sym={Ls['sstar']:.3e}")
    print(f"    anti: c={La['c']:.3e} b1={La['b1']:.3e} b2={La['b2']:.3e} "
          f"mu={La['mu']:.3e} eps={La['eps']:.3e}  s*_anti={La['sstar']:.3e}")
    sstar = min(Ls['sstar'], La['sstar'])
    chk(f"lemma margins strictly positive: lamhat={Ls['lamhat']:.3e} > 0, "
        f"mu_s={Ls['mu']:.3e} > 0, mu_a={La['mu']:.3e} > 0 (numeric echo of "
        "exact t0-gates)", Ls['lamhat'] > 0 and Ls['mu'] > 0 and La['mu'] > 0)
    chk(f"s* = {sstar:.3e} >= 1e-3 (both tested s inside the lemma interval)",
        sstar >= 1e-3)

    # ---- EXACT PD of both blocks at s = 1e-3, 1e-4 (the lemma conclusion) --
    for s in (R(1, 1000), R(1, 10000)):
        Ms_s = J['Ms0'] + s * Ms1 + s**2 * Ms2
        Ma_s = J['Ma0'] + s * Ma1 + s**2 * Ma2
        ps, pa = is_pd(Ms_s), is_pd(Ma_s)
        chk(f"EXACT PD at s = {s}: Ms(s) > 0 and Ma(s) > 0 (=> full "
            f"Gamma(y(s)) > 0 by swap congruence) -- NO o(s^2) repair",
            ps and pa)
        # numeric eigen margins (report; k=4 W-floor is below float resolution)
        es = np.linalg.eigvalsh(npm(Ms_s))
        ea = np.linalg.eigvalsh(npm(Ma_s))
        fl = float(s)**2 * min(Ls['mu'], La['mu']) / 2
        print(f"      numeric min eig: sym {es[0]:+.2e}, anti {ea[0]:+.2e}; "
              f"predicted W-floor ~ s^2 mu/2 = {fl:.2e}"
              + ("  [below float resolution: exact check is the arbiter]"
                 if fl < 1e-12 else ""))

    # ---- EXACT objective cubic ----
    y0, y1w, y2s, b0c, ca = (J['y0'], J['y1w'], J['y2s'], J['b0c'], J['ca'])
    NC = J['NC']
    y1 = [lam * v for v in y1w]
    y2 = [tt * v for v in y2s]
    c0 = sum(co * y0[c] for c, co in b0c.items())
    c1 = sum(co * y1[c] for c, co in b0c.items()) - y0[ca]
    g = sum(co * y2[c] for c, co in b0c.items()) - y1[ca]
    c3 = -y2[ca]
    chk(f"EXACT objective cubic: B_s.y(s) = 4 - s + g s^2 + (t0 lam^2) s^3 "
        f"with g = 1/(1024 t0) = {g} = {float(g):.3e}",
        c0 == 4 and c1 == -1 and g == 1 / (1024 * t0) and c3 == tt and g > 0)
    chk(f"=> a_{k} >= {g} > 0  (c_k(s) >= 4 - s + g s^2 + t0 lam^2 s^3 on "
        f"(0, s*]; liminf (c_k(s)-4+s)/s^2 >= g)", g > 0)
    print(f"    [t0={t0} lam={float(lam):.3e}; {time.time()-tstart:.1f}s]")
    return float(g), sstar


if __name__ == "__main__":
    T0 = time.time()
    res = {}
    ks = [int(a) for a in sys.argv[1:]] or [2, 3, 4]
    for k in ks:
        res[k] = run_level(k)
    print("\nframework-lemma bounds: " + ", ".join(
        f"a_{k} >= {g:.3e} (s* = {ss:.1e})" for k, (g, ss) in res.items()))
    print(f"\n[{time.time()-T0:.1f}s] "
          + ("ALL CHECKS PASS -- strict-margin jets are EXACTLY feasible for "
             "small s > 0: the o(s^2)-framework flag is DISCHARGED; "
             "R1(k) => a_k > 0." if OK else "CHECK FAILURES"))
    raise SystemExit(0 if OK else 1)
