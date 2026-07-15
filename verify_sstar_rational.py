"""RATIONAL CERTIFICATION of the no-repair constants and of s*_k.

Paper: npa_nonexactness_proof.tex, Lemma 6.1 (no-repair positivity),
Theorem 6.3 (exact arc), Remark 7.1 (explicit constants).
(Development record: pipeline_laws_note.md sec 9.)

verify_framework_lemma.py proves the STRUCTURAL hypotheses of Lemma 6.1
exactly (kernels, rank-one u-alignment, strict t0-margins, PD at two sample
points) but reports the lemma CONSTANTS (c, b1, b2, lamhat, h, kappa, mu,
eps) and the resulting interval endpoint s* in floating point.  This
verifier certifies every constant by a rational bound in EXACT arithmetic,
in the direction that is safe for s*, and outputs a certified rational
    s*_cert(k)  <=  s*(true constants),
so that Lemma 6.1's conclusion  M(s) >= 0 for all s in (0, s*_cert]  holds
with an explicit RATIONAL interval endpoint.  Bounds and their certificates
(per level k = 2,3,4 and per swap block; sym has dim U = 1, antisym U = 0):

  c   >= c_lo    exact rational PD of  Rb^T M0 Rb - c_lo * Rb^T Rb, where
                 Rb = rational basis of range(M0) = null(Nker^T)
                 (generalized-eigenvalue lower bound; c enters s* only
                 through lower bounds).
  b1  <= b1_up   exact rational PD of  b_up^2 I - M^2  (spectral norm of a
  b2  <= b2_up   symmetric M is <= b_up iff M^2 <= b_up^2 I).
  lamhat exact   U = ker(M0) cap W^perp has a RATIONAL direction v
                 (1-dim nullspace of (Ks^T GK), GK = Nker^T Nker);
                 lamhat = (v^T K1 v)/(v^T GK v), an exact rational > 0.
  h   exact      |H_UU| = |v^T HK v|/(v^T GK v), exact rational, where
                 HK = Nker^T (M2 - M1 M0^+ M1) Nker (exact M0^+).
  kappa <= k_up  kappa^2 = (r GW^{-1} r^T)/(v^T GK v) is an exact rational
                 (r = v^T HK Ks, GW = Ks^T GK Ks = exact Gram of the W
                 frame); k_up is a rational with k_up^2 >= kappa^2.
  mu  >= mu_lo   exact rational PD of  HW - mu_lo * GW  (HW = Ks^T HK Ks):
                 lower bound on the smallest eigenvalue of the W-block of
                 H in an orthonormal frame of W.
  eps <= eps_up  the Lemma-6.1 formula evaluated at (b1_up, b2_up, c_lo)
                 (monotone: increasing in b's, decreasing in c).
  s*  >= s*_cert the Lemma-6.1 minimum evaluated at the certified bounds
                 (each entry is monotone in the constants in the certified
                 direction), taken over both swap blocks.

Also re-asserts exactly: K1 = lam * u u^T on the symmetric kernel and
K1 = 0 on the antisymmetric kernel (hypothesis (b)), u^T Ks = 0.

GATES (exit 0 iff all pass): every PD certificate verifies exactly;
lamhat, mu_lo, c_lo > 0; s*_cert > 0 rational; s*_cert >= 1e-3 (so the
exact-PD sample points s = 1e-3, 1e-4 of verify_framework_lemma lie inside
the certified interval).  (The certified endpoints sit 5-12% below the
float reports of verify_framework_lemma, as they must.)
"""
import sys
import time
import numpy as np
import sympy as sp

sys.path.insert(0, '/Users/antonp/research/why-tsirelson')
from verify_small_lambda import is_pd, pinv_exact         # noqa: E402
from verify_framework_lemma import build_jet, npm, snorm  # noqa: E402

R = sp.Rational
OK = True


def chk(name, cond):
    global OK
    OK = OK and bool(cond)
    print(f"  [{'OK' if cond else 'FAIL'}] {name}")


def rat_below(x, rel=0.02, den=10**12):
    """rational strictly below the float x (x > 0)."""
    return R(int(x * (1 - rel) * den), den)


def rat_above(x, rel=0.02, den=10**12):
    """rational strictly above the float x (x >= 0)."""
    return R(int(x * (1 + rel) * den) + 1, den)


def cert_c_lo(M0, Nker):
    """certified rational lower bound on min spec(M0 | range M0)."""
    Rb = sp.Matrix.hstack(*Nker.T.nullspace())
    A = Rb.T * M0 * Rb
    B = Rb.T * Rb
    L = np.linalg.cholesky(npm(B))
    Mn = np.linalg.solve(L, np.linalg.solve(L, npm(A)).T)
    est = float(np.min(np.linalg.eigvalsh((Mn + Mn.T) / 2)))
    rel = 0.02
    for _ in range(8):
        c_lo = rat_below(est, rel)
        if c_lo > 0 and is_pd(A - c_lo * B):
            return c_lo
        rel = 1 - (1 - rel) / 2
    return None


def cert_bnorm_up(M):
    """certified rational upper bound on the spectral norm of symmetric M."""
    est = snorm(npm(M))
    n = M.rows
    rel = 0.02
    for _ in range(8):
        b_up = rat_above(est, rel)
        if is_pd(b_up**2 * sp.eye(n) - M * M):
            return b_up
        rel = 2 * rel + 0.02
    return None


def cert_mu_lo(HW, GW):
    """certified rational lower bound on the smallest generalized
    eigenvalue of (HW, GW), GW > 0."""
    L = np.linalg.cholesky(npm(GW))
    Mn = np.linalg.solve(L, np.linalg.solve(L, npm(HW)).T)
    est = float(np.min(np.linalg.eigvalsh((Mn + Mn.T) / 2)))
    rel = 0.02
    for _ in range(8):
        mu_lo = rat_below(est, rel)
        if mu_lo > 0 and is_pd(HW - mu_lo * GW):
            return mu_lo
        rel = 1 - (1 - rel) / 2
    return None


def rat_sqrt_up(q):
    """rational k_up with k_up^2 >= q (q rational >= 0)."""
    est = float(q) ** 0.5
    rel = 0.01
    for _ in range(12):
        k_up = rat_above(est, rel, den=10**15)
        if k_up**2 >= q:
            return k_up
        rel = 2 * rel + 0.01
    return None


def certify_block(name, M0, M1, M2, Nker, Ks, lam, u):
    """Certified rational s* for one swap block.  Ks = W-frame coordinates
    inside the kernel (None => W = whole kernel, U = 0)."""
    c_lo = cert_c_lo(M0, Nker)
    b1_up = cert_bnorm_up(M1)
    b2_up = cert_bnorm_up(M2)
    chk(f"{name}: c >= {sp.nsimplify(c_lo)} > 0, "
        f"b1 <= {float(b1_up):.3e}, b2 <= {float(b2_up):.3e} "
        "(exact PD certificates)",
        c_lo is not None and b1_up is not None and b2_up is not None)
    if not (c_lo and b1_up and b2_up):
        return None
    GK = Nker.T * Nker
    K1 = Nker.T * M1 * Nker
    _, M0p = pinv_exact(M0)
    HK = Nker.T * (M2 - M1 * M0p * M1) * Nker
    eps_up = (2 * b1_up**2 * (b1_up + b2_up) / c_lo**2
              + 4 * b1_up * b2_up / c_lo + 2 * b2_up**2 / c_lo)
    if Ks is None:
        # U = 0: hypothesis (b) is K1 == 0 exactly; W = whole kernel
        chk(f"{name}: kernel first-order form vanishes exactly (U = 0)",
            K1 == sp.zeros(K1.rows, K1.cols))
        mu_lo = cert_mu_lo(HK, GK)
        chk(f"{name}: mu >= {float(mu_lo):.3e} > 0 (exact PD of "
            "HW - mu_lo GW)", mu_lo is not None and mu_lo > 0)
        if mu_lo is None:
            return None
        sstar = min(R(1), c_lo / (2 * (b1_up + b2_up)),
                    mu_lo / (2 * eps_up))
        print(f"    {name}: eps <= {float(eps_up):.3e}  "
              f"s*_cert = {sstar} = {float(sstar):.4e}")
        return sstar
    # symmetric block: dim U = 1
    chk(f"{name}: K1 == lam u u^T exactly and u^T Ks == 0 "
        "(hypothesis (b), rank-one aligned)",
        K1 == lam * (u * u.T) and u.T * Ks == sp.zeros(1, Ks.cols))
    ns = sp.Matrix(Ks.T * GK).nullspace()
    chk(f"{name}: U-direction rational and unique "
        f"(dim null((Ks^T GK)) = {len(ns)})", len(ns) == 1)
    if len(ns) != 1:
        return None
    v = ns[0]
    nv2 = (v.T * GK * v)[0, 0]
    uv = (u.T * v)[0, 0]
    lamhat = lam * uv**2 / nv2
    hUU = sp.Abs((v.T * HK * v)[0, 0]) / nv2
    rW = v.T * HK * Ks
    GW = Ks.T * GK * Ks
    HW = Ks.T * HK * Ks
    kappa2 = (rW * GW.inv() * rW.T)[0, 0] / nv2
    k_up = rat_sqrt_up(kappa2)
    mu_lo = cert_mu_lo(HW, GW)
    chk(f"{name}: lamhat = {float(lamhat):.4e} exact > 0, h = "
        f"{float(hUU):.3e} exact, kappa <= {float(k_up):.3e} "
        f"(kappa^2 exact), mu >= {float(mu_lo):.3e} > 0 (exact PD)",
        lamhat > 0 and k_up is not None and mu_lo is not None and mu_lo > 0)
    if k_up is None or mu_lo is None:
        return None
    sstar = min(R(1), c_lo / (2 * (b1_up + b2_up)),
                lamhat / (2 * (hUU + eps_up)),
                mu_lo / (2 * eps_up),
                lamhat * mu_lo / (4 * (k_up + eps_up)**2))
    print(f"    {name}: eps <= {float(eps_up):.3e}  "
          f"s*_cert = {sstar} = {float(sstar):.4e}")
    return sstar


def run_level(k):
    print(f"\n===== level k = {k} =====")
    t = time.time()
    J = build_jet(k)   # exact structural gates re-run here
    t0 = J['t0']
    lam = 1 / (64 * t0)
    tt = t0 * lam**2
    ss = certify_block("sym ", J['Ms0'], lam * J['Ms1w'], tt * J['Ms2s'],
                       J['Ns'], J['Ks'], lam, J['u'])
    sa = certify_block("anti", J['Ma0'], lam * J['Ma1w'], tt * J['Ma2s'],
                       J['Na'], None, lam, None)
    if ss is None or sa is None:
        chk(f"k={k}: certification complete", False)
        return None
    sstar = min(ss, sa)
    chk(f"k={k}: certified rational s*_{k} = {sstar} = {float(sstar):.4e}"
        " > 0", sstar > 0)
    chk(f"k={k}: s*_cert >= 1e-3 (the exact-PD sample points s = 1e-3, "
        "1e-4 lie inside the certified interval)", sstar >= R(1, 1000))
    print(f"    [t0 = {t0}; {time.time()-t:.1f}s]")
    return sstar


if __name__ == "__main__":
    T0 = time.time()
    ks = [int(a) for a in sys.argv[1:]] or [2, 3, 4]
    res = {}
    for k in ks:
        res[k] = run_level(k)
    print("\ncertified rational interval endpoints: " + ", ".join(
        f"s*_{k} >= {v} ({float(v):.3e})" for k, v in res.items() if v))
    print(f"\n[{time.time()-T0:.1f}s] "
          + ("ALL RATIONAL s* CERTIFICATES PASS: Lemma 6.1 holds on "
             "(0, s*_cert] with explicit rational endpoints."
             if OK else "CERTIFICATION FAILURES"))
    raise SystemExit(0 if OK else 1)
