#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CLEAN-ROOM verification of the combinatorial core of
    npa_nonexactness_proof.tex   (Sections 2 and 4)
implemented from the PAPER TEXT ALONE.  No project code was consulted.
Pure Python stdlib; exact integer arithmetic throughout (y* is stored as 2*y*,
which is an integer by eq. (ystar)).

Conventions taken from the paper:
  * A word is an alternating tuple over {0,1} (letters X_0, X_1 of its party).
  * A word pair is (a,b), a Alice-word, b Bob-word;  W_k = {(a,b): |a|+|b|<=k}.
  * class(a,b): orbit of (a,b) under dagger (reverse BOTH parts) and party
    swap; canonical representative = lexicographic minimum of the 4 images.
  * Gamma_k(y)_{pq} = y[class(reduce(rev(a_p) a_q), reduce(rev(b_p) b_q))].
  * Trivial reductions T_q = e_q - e_{red q}, one per basis word q with a
    trailing 0 in either part; red q strips one trailing 0 from EACH part
    that ends in 0 (resolution of both-parts case: Lemma "tensor form",
    T = e_{(a,b)} - e_{(a',b')} when both end in 0).
  * Dressings D_W, W=(w_a,w_b), neither part ends in 1, |W|<=k-2:
    four corners W.(1, A_1, B_1, A_1B_1) with signs +,-,-,+ "written
    post-reduction".  Resolution (forced by Lemma "tensor form" and by the
    independence proof's claim that dressings are supported on IRREDUCIBLE
    words): each corner part additionally has its trailing 0 stripped, i.e.
    D_W = (e_{w_a'} - e_{w_a.1}) (x) (e_{w_b'} - e_{w_b.1}).
  * u_k:  u[T_{(a,b)}] = s(a)s(b) if exactly one part ends in 0, else 0,
          s(w) = -1 if w starts with 1 else +1 (empty -> +1);
          u[D_W] = 2(-1)^{|W|} for W != empty,  u[D_empty] = 0.
  * y*(alpha,beta) = i^2 j^2 - 2ij(ta*tb + sa*sb) + (ta+tb)/2
                     - (ta^2+tb^2)/2 + ta^2 tb^2 + (d0(alpha)+d0(beta))/2,
    i=#1s, tau=first+last-1, sigma=first-last, d0=[w==(0)], tau=sigma=0 on
    the empty word.

Exit code 0 iff every PASS/FAIL check passes.  The trailing D1 line is a
DIAGNOSTIC (naive reading of "post-reduction"), not a pass/fail check.
"""

import sys
from fractions import Fraction

# ----------------------------------------------------------------- words ---

def alt_words(length):
    """All alternating 0/1 words of the given length."""
    if length == 0:
        return [()]
    return [tuple((s + t) % 2 for t in range(length)) for s in (0, 1)]

def reduce_word(w):
    """Reduce a letter sequence modulo X^2 = 1 (stack cancellation)."""
    st = []
    for x in w:
        if st and st[-1] == x:
            st.pop()
        else:
            st.append(x)
    return tuple(st)

def rev(w):
    return tuple(reversed(w))

def mword(p, q):
    """Moment word pair of basis pairs p,q: per-party reduce(rev(.)*.)."""
    return (reduce_word(rev(p[0]) + q[0]), reduce_word(rev(p[1]) + q[1]))

def cls(wp):
    """Canonical class representative: min over dagger and party swap."""
    a, b = wp
    ra, rb = rev(a), rev(b)
    return min((a, b), (ra, rb), (b, a), (rb, ra))

IDENTITY = ((), ())

# ----------------------------------------------------------------- basis ---

def basis(k):
    out = []
    for la in range(k + 1):
        for lb in range(k + 1 - la):
            for a in alt_words(la):
                for b in alt_words(lb):
                    out.append((a, b))
    return out

def all_classes(k):
    """Classes of word pairs with total length <= 2k, mod dagger+swap."""
    seen = set()
    for la in range(2 * k + 1):
        for lb in range(2 * k + 1 - la):
            for a in alt_words(la):
                for b in alt_words(lb):
                    seen.add(cls((a, b)))
    return seen

# ------------------------------------------------------------------- y* ----

def _feats(w):
    i = sum(w)
    if not w:
        return i, 0, 0, 0
    tau = w[0] + w[-1] - 1
    sig = w[0] - w[-1]
    d0 = 1 if w == (0,) else 0
    return i, tau, sig, d0

def ystar2(wp):
    """2 * y*(alpha,beta): an integer, from eq. (ystar) of the paper."""
    a, b = wp
    ia, ta, sa, da = _feats(a)
    ib, tb, sb, db = _feats(b)
    return (2 * ia * ia * ib * ib
            - 4 * ia * ib * (ta * tb + sa * sb)
            + (ta + tb)
            - (ta * ta + tb * tb)
            + 2 * ta * ta * tb * tb
            + (da + db))

# -------------------------------------------------------- relation lattice --

def _s(w):
    """Sign s(w) of Definition (gain vector): -1 iff w starts with 1."""
    return -1 if (w and w[0] == 1) else 1

def trivial_reductions(k):
    """[(vector-dict, u-value, label)] for all T_q, q in W_k reducible."""
    out = []
    for q in basis(k):
        a, b = q
        ea = bool(a) and a[-1] == 0
        eb = bool(b) and b[-1] == 0
        if not (ea or eb):
            continue
        ra = a[:-1] if ea else a
        rb = b[:-1] if eb else b
        vec = {q: 1}
        red = (ra, rb)
        vec[red] = vec.get(red, 0) - 1
        u = _s(a) * _s(b) if (ea != eb) else 0
        out.append((vec, u, ('T', q)))
    return out

def _end0_word(length):
    """The unique alternating word of given length NOT ending in 1."""
    if length == 0:
        return ()
    first = 0 if length % 2 == 1 else 1
    return tuple((first + t) % 2 for t in range(length))

def _dcorner(w, strip):
    """D^Delta(w) = e_{w'} - e_{w.1}  (strip=True) or naive e_w - e_{w.1}."""
    w1 = w + (1,)
    w0 = (w[:-1] if (strip and w) else w)
    return [(w0, 1), (w1, -1)]

def dressings(k, strip=True):
    """[(vector-dict, u-value, label)] for all D_W, |W|<=k-2, no part ends 1."""
    out = []
    for la in range(k - 1):
        for lb in range(k - 1 - la):
            wa, wb = _end0_word(la), _end0_word(lb)
            vec = {}
            for (x, cx) in _dcorner(wa, strip):
                for (y, cy) in _dcorner(wb, strip):
                    key = (x, y)
                    vec[key] = vec.get(key, 0) + cx * cy
            u = 0 if (la == 0 and lb == 0) else 2 * (-1) ** (la + lb)
            out.append((vec, u, ('D', (wa, wb))))
    return out

def kernel_basis(k, strip=True):
    """N_k = [T's in basis order | D's in (l_a,l_b) order]."""
    return trivial_reductions(k) + dressings(k, strip)

# -------------------------------------------------------------- lin. alg. --

def rank_of_columns(cols, word_index):
    """Exact rank of the integer matrix whose columns are the given dicts."""
    rows = []
    for vec, _, _ in cols:
        r = [Fraction(0)] * len(word_index)
        for w, c in vec.items():
            r[word_index[w]] = Fraction(c)
        rows.append(r)
    rank = 0
    ncols = len(word_index)
    pivot_col = 0
    r = 0
    while r < len(rows) and pivot_col < ncols:
        piv = None
        for rr in range(r, len(rows)):
            if rows[rr][pivot_col] != 0:
                piv = rr
                break
        if piv is None:
            pivot_col += 1
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        pv = rows[r][pivot_col]
        for rr in range(r + 1, len(rows)):
            f = rows[rr][pivot_col] / pv
            if f:
                rows[rr] = [x - f * y for x, y in zip(rows[rr], rows[r])]
        rank += 1
        r += 1
        pivot_col += 1
    return rank

# ------------------------------------------------------------- main check --

class YCache:
    def __init__(self):
        self.c = {}
    def val2(self, p, q):
        m = mword(p, q)
        key = cls(m)
        v = self.c.get(key)
        if v is None:
            v = ystar2(key)
            self.c[key] = v
        return v

def witness_identity(kern, cache):
    """Check  N^T Gamma(y*) N == u u^T  exactly.  Returns (ok, first_bad)."""
    n = len(kern)
    for i in range(n):
        vi, ui, li = kern[i]
        items_i = list(vi.items())
        for j in range(i, n):
            vj, uj, lj = kern[j]
            e2 = 0
            for p, cp in items_i:
                for q, cq in vj.items():
                    e2 += cp * cq * cache.val2(p, q)
            if e2 != 2 * ui * uj:
                return False, (li, lj, Fraction(e2, 2), ui * uj)
    return True, None

# --------------------------------------------------------------- checks ----

RESULTS = []

def check(name, ok, detail=""):
    RESULTS.append(ok)
    print("%-58s %s%s" % (name, "PASS" if ok else "FAIL",
                          ("   " + detail) if detail else ""))

def main():
    KS = range(2, 8)

    # --- counting laws, rank, index-map sanity, per level -------------------
    for k in KS:
        B = basis(k)
        check("C1 k=%d  |W_k| = 2k^2+2k+1 = %d" % (k, 2 * k * k + 2 * k + 1),
              len(B) == 2 * k * k + 2 * k + 1, "got %d" % len(B))

        classes = all_classes(k)
        want = (5 * k + 2) * (k + 1) // 2
        check("C2 k=%d  #classes = (5k+2)(k+1)/2 = %d" % (k, want),
              len(classes) == want, "got %d" % len(classes))

        Ts = trivial_reductions(k)
        want = (3 * k * k + k) // 2
        check("C3 k=%d  #trivial reductions = (3k^2+k)/2 = %d" % (k, want),
              len(Ts) == want, "got %d" % len(Ts))

        Ds = dressings(k)
        want = k * (k - 1) // 2
        check("C4 k=%d  #dressings = k(k-1)/2 = %d" % (k, want),
              len(Ds) == want, "got %d" % len(Ds))

        kern = Ts + Ds
        widx = {w: t for t, w in enumerate(B)}
        supported = all(w in widx for vec, _, _ in kern for w in vec)
        rk = rank_of_columns(kern, widx) if supported else -1
        check("C5 k=%d  N_k columns independent (rank = 2k^2 = %d)"
              % (k, 2 * k * k),
              supported and len(kern) == 2 * k * k and rk == 2 * k * k,
              "rank %d, cols %d, in-basis %s" % (rk, len(kern), supported))

        # C6: Gamma(e_identity) = I ; C7: entry classes = class set
        ok6 = True
        seen = set()
        for p in B:
            for q in B:
                key = cls(mword(p, q))
                seen.add(key)
                if (key == IDENTITY) != (p == q):
                    ok6 = False
        check("C6 k=%d  Gamma_k(e_identity) = I" % k, ok6)
        check("C7 k=%d  Gamma entries exhaust the <=2k class set" % k,
              seen == classes,
              "" if seen == classes else
              "entries %d vs classes %d" % (len(seen), len(classes)))

    # --- y* well-definedness on classes (largest class set, k=7) ------------
    ok = True
    for (a, b) in all_classes(7):
        vals = {ystar2((a, b)), ystar2((rev(a), rev(b))),
                ystar2((b, a)), ystar2((rev(b), rev(a)))}
        if len(vals) != 1:
            ok = False
    check("C8 y* invariant under dagger and party swap (all classes, k=7)", ok)

    # --- pins ---------------------------------------------------------------
    ok = ystar2(IDENTITY) == 0
    for j in range(1, 8):
        beta_j = tuple((1 + t) % 2 for t in range(2 * j - 1))  # 1(01)^{j-1}
        ok = ok and beta_j[0] == 1 and beta_j[-1] == 1 \
                and ystar2(((), beta_j)) == 0
    check("C9 pins: y*[identity]=0 and y*[(empty,beta_j)]=0, j=1..7", ok)

    # --- stated consistency values (Lemma 'identities' at lam = 1) ----------
    yA0   = Fraction(ystar2(((0,), ())), 2)
    yA0B0 = Fraction(ystar2(((0,), (0,))), 2)
    yA0B1 = Fraction(ystar2(((0,), (1,))), 2)
    yA1B1 = Fraction(ystar2(((1,), (1,))), 2)
    yB1   = Fraction(ystar2(((), (1,))), 2)
    B0dot = 2 * yA0 + yA0B0 + 2 * yA0B1 - yA1B1
    check("C10a  y*_<A0> = -1/2", yA0 == Fraction(-1, 2), "got %s" % yA0)
    check("C10b  B_0 . y* = 0", B0dot == 0, "got %s" % B0dot)
    check("C10c  y*_<A0B0> = 0", yA0B0 == 0, "got %s" % yA0B0)
    check("C10d  y*_<A1B1> = 2 y*_<B1>", yA1B1 == 2 * yB1,
          "got %s vs %s" % (yA1B1, yB1))

    # --- tower property ------------------------------------------------------
    ok = True
    for k in range(3, 8):
        prev = {lab: u for _, u, lab in kernel_basis(k - 1)}
        cur = {lab: u for _, u, lab in kernel_basis(k)}
        for lab, u in prev.items():
            if lab not in cur or cur[lab] != u:
                ok = False
    check("C11 u_k extends u_{k-1} (labels and values), k=3..7", ok)

    # --- THE MAIN CHECK -------------------------------------------------------
    cache = YCache()
    for k in KS:
        kern = kernel_basis(k)
        ok, bad = witness_identity(kern, cache)
        check("C12 k=%d  N^T Gamma(y*) N = u u^T  (%d^2 kernel pairs)"
              % (k, len(kern)), ok,
              "" if ok else "first bad pair %s: E=%s, uu=%s"
              % ((bad[0], bad[1]), bad[2], bad[3]))

    # --- diagnostic: naive reading of "post-reduction" -----------------------
    diag = []
    for k in (3, 4):
        kern = kernel_basis(k, strip=False)
        ok, _ = witness_identity(kern, YCache())
        diag.append("k=%d:%s" % (k, "holds" if ok else "fails"))
    print("D1 [diagnostic, not scored] naive dressing corners "
          "(no trailing-0 strip): identity %s" % ", ".join(diag))

    n_pass = sum(RESULTS)
    print("-" * 72)
    print("TOTAL: %d/%d checks passed" % (n_pass, len(RESULTS)))
    return 0 if n_pass == len(RESULTS) else 1

if __name__ == "__main__":
    sys.exit(main())
