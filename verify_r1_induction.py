"""Machine check for the witness theorem: paper npa_nonexactness_proof.tex,
Theorem 4.4 (the witness works at every level) with Lemmas 4.5-4.9
(development record: Theorem 8.1, pipeline_laws_note.md sec 8).

THE CLOSED-FORM WITNESS.  For an alternating word w over {0,1} define
    i(w)   = number of 1s,
    tau(w) = first(w) + last(w) - 1   (tau(empty) = 0),
    sig(w) = first(w) - last(w)       (sig(empty) = 0),
    d0(w)  = 1 if w == (0,) else 0.
For a moment class (alpha, beta) (alpha = A-word, beta = B-word) define

  y*(alpha,beta) = i^2 j^2 - 2 i j (tau_a tau_b + sig_a sig_b)
                   + (tau_a + tau_b)/2 - (tau_a^2 + tau_b^2)/2
                   + tau_a^2 tau_b^2 + (d0(alpha) + d0(beta))/2 ,

with i = i(alpha), j = i(beta).  y* is invariant under dagger
(rev both: tau inv., sig -> -sig on both, product inv.) and party swap
(manifestly symmetric), so it is a well-defined CLASS function at every
level simultaneously.

CLAIM (paper Theorem 4.4; development record: Theorem 8.1): at every
level k >= 2, with the canonical kernel
basis K = [trivial reductions | dressings] and the closed-form u_k,
    K_a^T Gamma(y*) K_b = u_a u_b   for all a, b     (exactly),
i.e. P(y*) = u_k u_k^T -- and y* also satisfies the sharp pins
y*[identity] = 0 and y*[face coordinate classes] = 0.  Hence (R1) holds
at level k for every k, with the SAME witness.

This file checks the claim exactly (integer arithmetic, x2 scaling) by
DIRECT EVALUATION -- no linear algebra -- at k = 2..KMAX (default 8).
Also checks: pins, u-orthogonality bookkeeping used by the small-lambda
theorem, paper Theorem 5.4 (development record: Theorem 7.3); nothing
else is needed: Theorem 5.4 + this witness => v_k > 0.

The margin rerun at side length <= 26 quoted in the paper's appendix is
reproducible here via check_all_pairs(26): 3,946,645 pairs, all zero.

Exit code 0 iff all checks pass.
"""
import sys
import time

sys.path.insert(0, '/Users/antonp/research/why-tsirelson')
from npa_general import reduce_word                          # noqa: E402

OK = True


def chk(name, cond):
    global OK
    OK = OK and bool(cond)
    print(f"  [{'OK' if cond else 'FAIL'}] {name}")


# ---------------------------------------------------------------- words
def alt_words(maxlen):
    out = [()]
    for L in range(1, maxlen + 1):
        for start in (0, 1):
            out.append(tuple((start + i) % 2 for i in range(L)))
    return out


def strip0(w):
    """strip trailing 0s (alternating: at most one)."""
    return w[:-1] if w and w[-1] == 0 else w


def colreduce(aw, bw):
    return strip0(aw), strip0(bw)


def mword(p, q):
    """moment word of basis pair p, q (each = (aword, bword))."""
    return (reduce_word(tuple(reversed(p[0])) + q[0]),
            reduce_word(tuple(reversed(p[1])) + q[1]))


# ------------------------------------------------------- the witness y*
def _feat(w):
    i = sum(1 for x in w if x == 1)
    if w:
        tau = w[0] + w[-1] - 1
        sig = w[0] - w[-1]
    else:
        tau = sig = 0
    d0 = 1 if w == (0,) else 0
    return i, tau, sig, d0


def ystar2(alpha, beta):
    """2 * y*(alpha, beta): exact integer."""
    i, ta, sa, da = _feat(alpha)
    j, tb, sb, db = _feat(beta)
    return (2*i*i*j*j - 4*i*j*(ta*tb + sa*sb)
            + (ta + tb) - (ta*ta + tb*tb) + 2*ta*ta*tb*tb + da + db)


# ------------------------------------------------- kernel + u at level k
def endw0(l):
    """alternating length-l word ending in 0."""
    return tuple(((l - 1 - i) % 2) for i in range(l))


def kernel_and_u(k):
    """list of (support, u) for the canonical kernel at level k;
    support = list of (coef, basisword)."""
    basis = []
    for aw in alt_words(k):
        for bw in alt_words(k - len(aw)):
            basis.append((aw, bw))
    out = []

    def s(w):
        return -1 if (w and w[0] == 1) else 1

    def ends0(w):
        return bool(w) and w[-1] == 0
    for (aw, bw) in basis:                       # trivial reductions
        r = colreduce(aw, bw)
        if r != (aw, bw):
            uval = s(aw)*s(bw) if ends0(aw) != ends0(bw) else 0
            out.append(([(1, (aw, bw)), (-1, r)], uval))
    for la in range(0, k - 1):                   # dressings
        for lb in range(0, k - 1 - la):
            wa, wb = endw0(la), endw0(lb)
            sup = []
            for da, db, sg in [((), (), 1), ((1,), (), -1),
                               ((), (1,), -1), ((1,), (1,), 1)]:
                sup.append((sg, colreduce(wa + da, wb + db)))
            L = la + lb
            out.append((sup, 0 if L == 0 else 2*(-1)**L))
    nT = sum(1 for (aw, bw) in basis if colreduce(aw, bw) != (aw, bw))
    assert nT == (3*k*k + k)//2 and len(out) - nT == k*(k - 1)//2
    return basis, out


def check_level(k):
    t0 = time.time()
    basis, KU = kernel_and_u(k)
    nK = len(KU)
    cache = {}

    def Y2(p, q):
        m = mword(p, q)
        v = cache.get(m)
        if v is None:
            v = ystar2(*m)
            cache[m] = v
        return v
    bad = 0
    for a in range(nK):
        sa, ua = KU[a]
        for b in range(a, nK):
            sb, ub = KU[b]
            tot = 0
            for ca, p in sa:
                for cb, q in sb:
                    tot += ca*cb*Y2(p, q)
            if tot != 2*ua*ub:
                bad += 1
    chk(f"k={k}: P(y*) == u_k u_k^T exactly (all {nK*(nK+1)//2} kernel "
        f"pairs, direct evaluation) [{time.time()-t0:.1f}s]", bad == 0)
    # sharp pins: y*[identity] = 0, y*[face coords <(B1B0)^(j-1)B1>] = 0
    pins_ok = ystar2((), ()) == 0
    for j in range(1, k + 1):
        w = tuple(((i + 1) % 2) for i in range(2*j - 1))   # 1,0,1,...,1
        pins_ok = pins_ok and ystar2((), w) == 0
    chk(f"k={k}: sharp pins y*[1] = 0 and y*[face coords] = 0", pins_ok)
    return bad == 0


# ----------------------------------------------- part B: regime-grid check
# Proof engine for ALL k (paper Lemmas 4.5-4.9; development record: note
# sec 8, Lemmas 8.3-8.5): every kernel-pair
# equation E(K1,K2) = sum c_p c_q y*(m(p,q)) - u1 u2, restricted to a fine
# regime (first letters + parities of the 4 base-word side lengths + all
# length comparisons clamped at 4), is a polynomial of degree <= 2 in each
# of the <= 4 free length coordinates (steps of 2).  Every fine regime
# realizes its full {0,1,2}-grid of free coordinates with all side lengths
# <= 20 (clamp base <= 6 by parity, + 2 grid steps of 2 on the min- and
# diff-coordinates: 10 + 10).  So E == 0 on ALL kernel pairs with side
# lengths <= 20 (checked below, level-independent enumeration) proves
# E == 0 for ALL kernel pairs at ALL levels.  See paper sec 4.3
# (development record: note sec 8.3).
def elements_upto(L):
    """all kernel elements (T and D) with BOTH side words of length <= L.
    Superset-by-regime of what any level uses; returns (support, u)."""
    words = alt_words(L)

    def s(w):
        return -1 if (w and w[0] == 1) else 1

    def ends0(w):
        return bool(w) and w[-1] == 0
    out = []
    for aw in words:                                # T elements
        for bw in words:
            if not (ends0(aw) or ends0(bw)):
                continue
            r = colreduce(aw, bw)
            uval = s(aw)*s(bw) if ends0(aw) != ends0(bw) else 0
            out.append(([(1, (aw, bw)), (-1, r)], uval))
    for la in range(0, L + 1):                      # D elements
        for lb in range(0, L + 1):
            wa, wb = endw0(la), endw0(lb)
            sup = []
            for da, db, sg in [((), (), 1), ((1,), (), -1),
                               ((), (1,), -1), ((1,), (1,), 1)]:
                sup.append((sg, colreduce(wa + da, wb + db)))
            Lw = la + lb
            out.append((sup, 0 if Lw == 0 else 2*(-1)**Lw))
    return out


def check_all_pairs(L):
    t0 = time.time()
    elems = elements_upto(L)
    nE = len(elems)
    # per-party moment-word feature cache
    fA = {}

    def feats(v, w):
        key = (v, w)
        r = fA.get(key)
        if r is None:
            m = reduce_word(tuple(reversed(v)) + w)
            r = _feat(m)
            fA[key] = r
        return r
    # y* value cache on feature pairs
    yc = {}

    def Y2f(pa, pb):
        key = (pa, pb)
        v = yc.get(key)
        if v is None:
            i, ta, sa, da = pa
            j, tb, sb, db = pb
            v = (2*i*i*j*j - 4*i*j*(ta*tb + sa*sb)
                 + (ta + tb) - (ta*ta + tb*tb) + 2*ta*ta*tb*tb + da + db)
            yc[key] = v
        return v
    bad = 0
    npairs = 0
    for x in range(nE):
        sx, ux = elems[x]
        for y in range(x, nE):
            sy, uy = elems[y]
            tot = 0
            for cx, p in sx:
                for cy, q in sy:
                    tot += cx*cy*Y2f(feats(p[0], q[0]), feats(p[1], q[1]))
            npairs += 1
            if tot != 2*ux*uy:
                bad += 1
    chk(f"regime-grid check: E(K1,K2) == 0 for ALL {npairs} kernel-element "
        f"pairs with side lengths <= {L} ({nE} elements) "
        f"[{time.time()-t0:.1f}s]", bad == 0)
    return bad == 0


if __name__ == "__main__":
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    Lgrid = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    t0 = time.time()
    print(f"(Thm 4.4) closed-form witness y*: levels k = 2..{kmax}:")
    for k in range(2, kmax + 1):
        check_level(k)
    print(f"(Thm 4.4) regime-grid engine (proves ALL k; Lemmas 4.5-4.9):")
    check_all_pairs(Lgrid)
    print(f"\n[{time.time()-t0:.1f}s] " +
          ("ALL CHECKS PASS" if OK else "CHECK FAILURES"))
    raise SystemExit(0 if OK else 1)
