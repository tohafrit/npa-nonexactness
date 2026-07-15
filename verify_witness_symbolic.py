#!/usr/bin/env python3
"""FULLY SYMBOLIC verification of Theorem 8.1 (pipeline_laws_note.md sec 8,
paper Theorem "the witness works at every level"): the regime step done by
machine, with NO interpolation argument.

WHAT IS PROVED HERE.  For every pair of kernel elements (T or D, any base
side lengths, at any level) the equation
    E(K1,K2) = sum_{p,q} c_p c_q y*(m(p,q)) - u(K1) u(K2) = 0
is verified as follows:

  (1) REGIME ENUMERATION.  An alternating word is determined by its first
      letter F and length L.  Per party, the ordered pair of base side
      lengths (l1, l2) of (K1, K2) is classified by
          mu = min(l1,l2):  exact value in [0..3], or "large" = M0 + 2n,
                            M0 in {4,5} (parity), n >= 0 a free symbol;
          d  = l1 - l2:     exact value in [-3..3], or "large" =
                            +-(D0 + 2m), D0 in {4,5}, m >= 0 free.
      A REGIME = (kind of K1, kind of K2) in {T,D}^2, one such joint
      length class per party, and the first letters of the four base side
      words (free in {0,1} for T sides of nonzero length; empty sides have
      no letter; D-side letters are parity-forced).  The classifier
      classify_lengths() is a TOTAL function: every ordered pair of valid
      kernel elements, with arbitrary side lengths, lands in exactly one
      enumerated regime at some nonnegative integer symbol values.

  (2) SYMBOLIC ENGINE.  Within a regime the lengths of all support words
      (base / strip-trailing-0 / append-1) are affine in the free symbols
      (n_A, m_A, n_B, m_B) with even symbol coefficients, so every length
      has a definite parity.  The moment word m(v,w) = reduce(rev(v).w) is
      computed by the locality lemma (Lemma 8.2): its first letter and
      length are produced by a case analysis whose every branch decision
      (emptiness, first-letter equality, length trichotomy, L == 1 for the
      d0 flag) is made by a GUARDED sign oracle that is sound for all
      nonnegative symbol values and RAISES if a decision is not uniform
      over the regime.  No decision ever raises (checked), so the word
      features (i, tau, sigma, d0) are exact affine/constant functions of
      the symbols, y*(m(p,q)) is an exact polynomial, and
      E(K1,K2) is assembled as an exact sympy polynomial in <= 4 symbols.

  (3) IDENTITY CHECK.  sympy expand(E) == 0 in EVERY regime.  Since the
      classifier is total and the engine is sound pointwise, E == 0 for
      every kernel pair at every level: Theorem 8.1, with no interpolation
      / degree-bound / grid-coverage reasoning in the trust base.

  (4) CROSS-CHECKS against exact integer ground truth (self-contained
      re-implementation, no imports from the repo):
        - every concrete kernel-element pair with side lengths <= LCC
          (default 20; 1,413,721 unordered pairs as in the grid check)
          has E == 0 exactly (integers) AND classifies into an enumerated
          regime; the set of realized regimes EQUALS the enumerated set;
        - on a random sample of pairs, the symbolic support words, moment
          words, features and u's, instantiated at the classifier's symbol
          values, coincide with the concrete ones.

Exit code 0 iff every regime's E is identically zero and all cross-checks
pass.  Usage:  python3 verify_witness_symbolic.py [symbolic|concrete|all]
                                                  [LCC] [sample]
"""
import random
import sys
import time

import sympy as sp

nA, mA, nB, mB = sp.symbols('n_A m_A n_B m_B', integer=True,
                            nonnegative=True)
SYMS = (nA, mA, nB, mB)
CL = 3                      # exact-window clamp: mu in [0..CL], |d| <= CL
BASES = (CL + 1, CL + 2)    # "large" bases by parity: 4, 5

OK = True


def chk(name, cond):
    global OK
    OK = OK and bool(cond)
    print(f"  [{'OK' if cond else 'FAIL'}] {name}")


class Undecidable(Exception):
    """A regime branch decision would depend on the free symbols."""


# ------------------------------------------------ guarded affine oracle
_AFF = {}


def affine_parts(e):
    r = _AFF.get(e)
    if r is None:
        ee = sp.expand(e)
        c0 = ee.subs({s: 0 for s in SYMS})
        cs = tuple(ee.coeff(s) for s in SYMS)
        r = (int(c0), tuple(int(c) for c in cs))
        assert sp.expand(ee - c0 - sum(c*s for c, s in zip(r[1], SYMS))) == 0
        _AFF[e] = r
    return r


def sign_of(e):
    """-1 / 0 / +1, valid for ALL nonneg integer symbol values, else raise."""
    c0, cs = affine_parts(sp.sympify(e))
    if c0 == 0 and all(c == 0 for c in cs):
        return 0
    if c0 > 0 and all(c >= 0 for c in cs):
        return 1
    if c0 < 0 and all(c <= 0 for c in cs):
        return -1
    raise Undecidable(f"sign of {e}")


def parity_of(e):
    c0, cs = affine_parts(sp.sympify(e))
    if not all(c % 2 == 0 for c in cs):
        raise Undecidable(f"parity of {e}")
    return c0 % 2


# ------------------------------------- symbolic words:  (F, L),  F=None
# iff the word is empty.  Alternating word = (first letter, length).
def W(F, L):
    L = sp.expand(sp.sympify(L))
    s = sign_of(L)
    assert s >= 0, "negative word length"
    if s == 0:
        return (None, sp.Integer(0))
    assert F in (0, 1)
    return (F, L)


def lastletter(F, L):
    return (F + parity_of(L) + 1) % 2


def ends0(F, L):
    return F is not None and lastletter(F, L) == 0


def strip_word(F, L):
    """strip one trailing 0 if present (alternating: at most one)."""
    if not ends0(F, L):
        return (F, L)
    return W(F, L - 1) if sign_of(L - 1) > 0 else (None, sp.Integer(0))


def s_of(F):
    return -1 if F == 1 else 1


def mword(w1, w2):
    """(F, L) of reduce(rev(v) . w)  -- Lemma 8.2 case analysis."""
    F1, L1 = w1
    F2, L2 = w2
    if F1 is None and F2 is None:
        return (None, sp.Integer(0))
    if F1 is None:
        return w2
    revF1 = lastletter(F1, L1)
    if F2 is None:
        return (revF1, L1)
    if F1 != F2:                       # no cancellation
        return (revF1, sp.expand(L1 + L2))
    t = sign_of(L1 - L2)               # full common-prefix cancellation
    if t == 0:
        return (None, sp.Integer(0))
    if t < 0:
        return ((F2 + parity_of(L1)) % 2, sp.expand(L2 - L1))
    return (revF1, sp.expand(L1 - L2))


def feats(w):
    """(i, tau, sig, d0) of an alternating word; i is a sympy expr."""
    F, L = w
    if F is None:
        return (sp.Integer(0), 0, 0, 0)
    last = lastletter(F, L)
    tau = F + last - 1
    sig = F - last
    i = sp.expand((L + tau) / 2)
    if F == 0:
        d0 = 1 if sign_of(L - 1) == 0 else 0   # raises if not uniform
    else:
        d0 = 0
    return (i, tau, sig, d0)


# ------------------------------------------------------ kernel elements
def t_element(wa, wb):
    """T element from base sides; None if invalid (nothing to reduce)."""
    ea, eb = ends0(*wa), ends0(*wb)
    if not (ea or eb):
        return None
    sup = [(1, (wa, wb)), (-1, (strip_word(*wa), strip_word(*wb)))]
    u = s_of(wa[0]) * s_of(wb[0]) if (ea != eb) else 0
    return (sup, u)


def dd_side(L):
    """[strip, append] support words of a dressing side endw0(L)."""
    if sign_of(L) == 0:
        return [(None, sp.Integer(0)), (1, sp.Integer(1))]
    F = (parity_of(L) + 1) % 2         # first letter of endw0(L)
    stripped = W(F, L - 1) if sign_of(L - 1) > 0 else (None, sp.Integer(0))
    return [stripped, (F, sp.expand(L + 1))]


def d_element(La, Lb):
    A, B = dd_side(La), dd_side(Lb)
    sup = [(1, (A[0], B[0])), (-1, (A[1], B[0])),
           (-1, (A[0], B[1])), (1, (A[1], B[1]))]
    if sign_of(La) == 0 and sign_of(Lb) == 0:
        u = 0
    else:
        u = 2 * (-1) ** ((parity_of(La) + parity_of(Lb)) % 2)
    return (sup, u)


def make_element(kind, fa, la, fb, lb):
    if kind == 'T':
        return t_element(W(fa, la), W(fb, lb))
    return d_element(la, lb)


# ------------------------------------------------------- E of a regime
_YC = {}


def Y2(fA, fB):
    """2 * y* on a feature pair (exact; matches verify_r1_induction.py)."""
    key = (fA, fB)
    v = _YC.get(key)
    if v is None:
        i, ta, sa, da = fA
        j, tb, sb, db = fB
        v = (2*i*i*j*j - 4*i*j*(ta*tb + sa*sb)
             + (ta + tb) - (ta*ta + tb*tb) + 2*ta*ta*tb*tb + da + db)
        _YC[key] = v
    return v


_FC = {}


def pair_feats(w1, w2):
    key = (w1, w2)
    r = _FC.get(key)
    if r is None:
        r = feats(mword(w1, w2))
        _FC[key] = r
    return r


def E_of(el1, el2):
    sup1, u1 = el1
    sup2, u2 = el2
    tot = sp.Integer(0)
    for c1, (a1, b1) in sup1:
        for c2, (a2, b2) in sup2:
            tot += c1 * c2 * Y2(pair_feats(a1, a2), pair_feats(b1, b2))
    return sp.expand(tot - 2 * u1 * u2)


# --------------------------------------------------- regime enumeration
def party_classes(n, m):
    """all joint length classes (id, l1, l2) for one party; n, m = its
    free symbols (min- and difference-coordinate, stepping by 2)."""
    mus = [('e', v) for v in range(CL + 1)] + [('L', b) for b in BASES]
    dts = ([('e', v) for v in range(-CL, CL + 1)]
           + [('L', sg * b) for b in BASES for sg in (1, -1)])
    out = []
    for mu in mus:
        lmin = sp.Integer(mu[1]) if mu[0] == 'e' else mu[1] + 2 * n
        for dt in dts:
            if dt[0] == 'e':
                l1 = sp.expand(lmin + max(dt[1], 0))
                l2 = sp.expand(lmin + max(-dt[1], 0))
            elif dt[1] > 0:
                l1 = sp.expand(lmin + dt[1] + 2 * m)
                l2 = lmin
            else:
                l1 = lmin
                l2 = sp.expand(lmin - dt[1] + 2 * m)
            out.append(((mu, dt), l1, l2))
    return out


def letter_opts(kind, L):
    if kind == 'D':
        return [None]                  # parity-forced inside d_element
    return [None] if sign_of(L) == 0 else [0, 1]


def run_symbolic(verbose=True):
    t0 = time.time()
    clsA = party_classes(nA, mA)
    clsB = party_classes(nB, mB)
    enumerated = {}                    # key -> nfree
    nonzero, undecid = [], []
    for k1 in ('T', 'D'):
        for k2 in ('T', 'D'):
            for idA, l1a, l2a in clsA:
                for idB, l1b, l2b in clsB:
                    for f1a in letter_opts(k1, l1a):
                        for f1b in letter_opts(k1, l1b):
                            el1 = make_element(k1, f1a, l1a, f1b, l1b)
                            if el1 is None:
                                continue
                            for f2a in letter_opts(k2, l2a):
                                for f2b in letter_opts(k2, l2b):
                                    el2 = make_element(k2, f2a, l2a,
                                                       f2b, l2b)
                                    if el2 is None:
                                        continue
                                    key = (k1, k2, idA, idB,
                                           f1a, f1b, f2a, f2b)
                                    try:
                                        e = E_of(el1, el2)
                                    except Undecidable as ex:
                                        undecid.append((key, str(ex)))
                                        continue
                                    nf = len((l1a + l2a + l1b + l2b)
                                             .free_symbols)
                                    enumerated[key] = nf
                                    if e != 0:
                                        nonzero.append((key, e))
    nreg = len(enumerated)
    byfree = {}
    for v in enumerated.values():
        byfree[v] = byfree.get(v, 0) + 1
    if verbose:
        print(f"  regimes enumerated: {nreg}  "
              f"(by #free symbols: {dict(sorted(byfree.items()))})")
        for key, e in nonzero[:10]:
            print(f"    NONZERO regime {key}: E = {e}")
        for key, ex in undecid[:10]:
            print(f"    UNDECIDABLE regime {key}: {ex}")
    chk(f"symbolic sweep: E(K1,K2) identically 0 (sympy expand) in ALL "
        f"{nreg} regimes, 0 undecidable [{time.time()-t0:.1f}s]",
        not nonzero and not undecid and nreg > 0)
    return enumerated


# ===================================================== concrete ground truth
# (self-contained: no imports from the repo)
def reduce_c(w):
    out = []
    for x in w:
        if out and out[-1] == x:
            out.pop()
        else:
            out.append(x)
    return tuple(out)


def altword(F, L):
    return tuple((F + i) % 2 for i in range(L))


def alt_words_c(maxlen):
    return [()] + [altword(F, L) for L in range(1, maxlen + 1)
                   for F in (0, 1)]


def strip0_c(w):
    return w[:-1] if w and w[-1] == 0 else w


def feats_c(w):
    i = sum(1 for x in w if x == 1)
    tau = (w[0] + w[-1] - 1) if w else 0
    sig = (w[0] - w[-1]) if w else 0
    return (i, tau, sig, 1 if w == (0,) else 0)


def elements_c(L):
    """(support, u, meta); meta = (kind, la, fa, lb, fb)."""
    out = []
    for aw in alt_words_c(L):                       # T elements
        for bw in alt_words_c(L):
            ea = bool(aw) and aw[-1] == 0
            eb = bool(bw) and bw[-1] == 0
            if not (ea or eb):
                continue
            sup = [(1, (aw, bw)), (-1, (strip0_c(aw), strip0_c(bw)))]
            sa = -1 if (aw and aw[0] == 1) else 1
            sb = -1 if (bw and bw[0] == 1) else 1
            u = sa * sb if (ea != eb) else 0
            out.append((sup, u, ('T', len(aw), aw[0] if aw else None,
                                 len(bw), bw[0] if bw else None)))
    for la in range(0, L + 1):                      # D elements
        for lb in range(0, L + 1):
            wa = altword((la - 1) % 2, la) if la else ()
            wb = altword((lb - 1) % 2, lb) if lb else ()
            sup = []
            for da, db, sg in [((), (), 1), ((1,), (), -1),
                               ((), (1,), -1), ((1,), (1,), 1)]:
                sup.append((sg, (strip0_c(wa + da), strip0_c(wb + db))))
            u = 0 if la + lb == 0 else 2 * (-1) ** (la + lb)
            out.append((sup, u, ('D', la, None, lb, None)))
    return out


def classify_lengths(l1, l2):
    """total classifier: ordered length pair -> class id + (n, m) values."""
    mu, d = min(l1, l2), l1 - l2
    if mu <= CL:
        mu_c, nv = ('e', mu), 0
    else:
        b = BASES[0] + ((mu - BASES[0]) % 2)
        mu_c, nv = ('L', b), (mu - b) // 2
    if abs(d) <= CL:
        d_c, mv = ('e', d), 0
    else:
        b = BASES[0] + ((abs(d) - BASES[0]) % 2)
        d_c, mv = ('L', b if d > 0 else -b), (abs(d) - b) // 2
    return (mu_c, d_c), (nv, mv)


def classify_pair(m1, m2):
    k1, la1, fa1, lb1, fb1 = m1
    k2, la2, fa2, lb2, fb2 = m2
    idA, (nv, mv) = classify_lengths(la1, la2)
    idB, (nv2, mv2) = classify_lengths(lb1, lb2)
    return ((k1, k2, idA, idB, fa1, fb1, fa2, fb2), (nv, mv, nv2, mv2))


def Y2c(fA, fB):
    i, ta, sa, da = fA
    j, tb, sb, db = fB
    return (2*i*i*j*j - 4*i*j*(ta*tb + sa*sb)
            + (ta + tb) - (ta*ta + tb*tb) + 2*ta*ta*tb*tb + da + db)


def run_concrete(enumerated, LCC, nsample):
    t0 = time.time()
    elems = elements_c(LCC)
    nE = len(elems)
    fcache = {}

    def pf(v, w):
        key = (v, w)
        r = fcache.get(key)
        if r is None:
            r = feats_c(reduce_c(tuple(reversed(v)) + w))
            fcache[key] = r
        return r
    bad = badkey = 0
    realized = set()
    npairs = 0
    for x in range(nE):
        sx, ux, mx = elems[x]
        for y in range(x, nE):
            sy, uy, my = elems[y]
            tot = 0
            for cx, (pa, pb) in sx:
                for cy, (qa, qb) in sy:
                    fA, fB = pf(pa, qa), pf(pb, qb)
                    tot += cx * cy * Y2c(fA, fB)
            npairs += 1
            if tot != 2 * ux * uy:
                bad += 1
            for mm1, mm2 in ((mx, my), (my, mx)):
                key, _ = classify_pair(mm1, mm2)
                realized.add(key)
                if key not in enumerated:
                    badkey += 1
    chk(f"concrete ground truth: E == 0 (exact integers) for all {npairs} "
        f"unordered element pairs, side lengths <= {LCC} "
        f"[{time.time()-t0:.1f}s]", bad == 0)
    chk(f"classifier totality: all {nE*nE} ordered pairs land in "
        f"enumerated regimes", badkey == 0)
    missing = set(enumerated) - realized
    extra = realized - set(enumerated)
    chk(f"coverage accounting: realized regimes ({len(realized)}) == "
        f"enumerated regimes ({len(enumerated)})",
        not missing and not extra)
    if missing:
        print(f"    e.g. enumerated-but-unrealized: {sorted(missing)[:5]}")

    # sampled deep cross-check: symbolic engine vs concrete words
    t1 = time.time()
    rng = random.Random(20260715)
    clsA = {i: (l1, l2) for i, l1, l2 in party_classes(nA, mA)}
    clsB = {i: (l1, l2) for i, l1, l2 in party_classes(nB, mB)}
    mismatch = 0
    for _ in range(nsample):
        sx, ux, mx = elems[rng.randrange(nE)]
        sy, uy, my = elems[rng.randrange(nE)]
        key, vals = classify_pair(mx, my)
        subs = dict(zip(SYMS, vals))
        (k1, k2, idA, idB, f1a, f1b, f2a, f2b) = key
        l1a, l2a = clsA[idA]
        l1b, l2b = clsB[idB]
        el1 = make_element(k1, f1a, l1a.subs(subs), f1b, l1b.subs(subs))
        el2 = make_element(k2, f2a, l2a.subs(subs), f2b, l2b.subs(subs))
        # instantiated symbolic elements == concrete elements?
        for (el, (sc, uc, _)) in ((el1, (sx, ux, mx)), (el2, (sy, uy, my))):
            sup, u = el
            if u != uc:
                mismatch += 1
            for (c, (wA, wB)), (cc, (va, vb)) in zip(sup, sc):
                for (F, L), v in ((wA, va), (wB, vb)):
                    if c != cc or int(L) != len(v) or \
                            F != (v[0] if v else None):
                        mismatch += 1
        # per-pairing symbolic moment words == concrete reduce()?
        for (c1, (a1, b1)), (cc1, (pa, pb)) in zip(el1[0], sx):
            for (c2, (a2, b2)), (cc2, (qa, qb)) in zip(el2[0], sy):
                for (wm, vm) in ((mword(a1, a2),
                                  reduce_c(tuple(reversed(pa)) + qa)),
                                 (mword(b1, b2),
                                  reduce_c(tuple(reversed(pb)) + qb))):
                    F, L = wm
                    if int(L) != len(vm) or F != (vm[0] if vm else None):
                        mismatch += 1
    chk(f"sampled instantiation: symbolic supports, u's and moment words "
        f"reproduce concrete ground truth on {nsample} random pairs "
        f"[{time.time()-t1:.1f}s]", mismatch == 0)


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else 'all'
    LCC = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    nsample = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
    t0 = time.time()
    print("(sec 8.5) symbolic regime verification of Theorem 8.1:")
    enumerated = run_symbolic() if stage in ('all', 'symbolic') else None
    if stage in ('all', 'concrete'):
        if enumerated is None:
            enumerated = run_symbolic(verbose=False)
        run_concrete(enumerated, LCC, nsample)
    print(f"\n[{time.time()-t0:.1f}s] "
          + ("ALL CHECKS PASS" if OK else "CHECK FAILURES"))
    raise SystemExit(0 if OK else 1)
