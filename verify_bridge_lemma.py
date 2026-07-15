"""Machine check of the swap-split bridge lemma
(npa_nonexactness_proof.tex, Lemma 5.1 [swap-split bridge]), k = 2..5.
Uses only the fresh (audit) implementation of words/kernel/u/y*, not the
certificate bases.

Claims:
 (1) S (party swap on basis-word pairs) permutes the canonical T/D kernel
     basis: S.T_{(a,b)} = T_{(b,a)}, S.D_{(wa,wb)} = D_{(wb,wa)}  (exactly,
     as coefficient vectors on basis words).
 (2) u is invariant under the induced permutation pi: u(pi(K)) = u(K).
 (3) Bridge: with Lambda_s = span of {K + S.K} (+ fixed K), Lambda_a =
     span of {K - S.K}: the compression of Gamma(y*) to Lambda_a is 0, and
     to Lambda_s is  utilde utilde^T  with utilde_i = u . c^{(i)}
     (coefficients of the i-th symmetric combination in the kernel basis).
 (4) dim Lambda_s + dim Lambda_a = 2k^2 (split is exhaustive), and the
     symmetric/antisym combination sets are each linearly independent.
"""
import sys
sys.path.insert(0, '/Users/antonp/research/why-tsirelson')
from audit_spot import (red, altw, strip0, ends0, endw0, T_elem, D_elem,
                        level_kernel, mom, ystar4, sgn)
from fractions import Fraction

def swap_support(sup):
    return sorted((c, (q[1], q[0])) for c, q in sup)

def norm_sup(sup):
    return sorted((c, q) for c, q in sup)

ok = True
for k in range(2, 6):
    ker = level_kernel(k)
    n = len(ker)
    assert n == 2 * k * k
    # identify each kernel element by its normalized support
    key2idx = {tuple(norm_sup(s)): i for i, (s, u) in enumerate(ker)}
    pi = []
    for s, u in ker:
        j = key2idx.get(tuple(swap_support(s)))
        assert j is not None, "S does not permute the kernel basis!"
        pi.append(j)
    # (1) check the *labelled* form of the permutation: rebuild expected partner
    lab_ok = True
    # reconstruct labels: T elements first (aw,bw enumeration), then D (la,lb)
    labels = []
    for aw in altw(k):
        for bw in altw(k - len(aw)):
            if ends0(aw) or ends0(bw):
                labels.append(('T', aw, bw))
    for la in range(0, k - 1):
        for lb in range(0, k - 1 - la):
            labels.append(('D', la, lb))
    assert len(labels) == n
    for i, lab in enumerate(labels):
        if lab[0] == 'T':
            exp_sup, exp_u = T_elem(lab[2], lab[1])
        else:
            exp_sup, exp_u = D_elem(lab[2], lab[1])
        j = pi[i]
        if norm_sup(exp_sup) != norm_sup(ker[j][0]) or exp_u != ker[j][1]:
            lab_ok = False
    # (2) u invariance
    u_ok = all(ker[pi[i]][1] == ker[i][1] for i in range(n))
    # involution
    inv_ok = all(pi[pi[i]] == i for i in range(n))
    # (3) build symmetric / antisymmetric combinations as coefficient
    # vectors c over the kernel basis; Gram of Gamma(y*) on lattice =
    # c^T (u u^T + 0) c ... but compute *directly* from y* pairings:
    # <K_c, K_c'>_{y*} = sum_m,m' c_m c'_m' E_pair-ish full value
    # pair value P[m,m'] = sum over supports of y*(mom) (times 4)
    P = [[None] * n for _ in range(n)]
    def pval(i, j):
        if P[i][j] is None:
            s1, _ = ker[i]; s2, _ = ker[j]
            tot = 0
            for c1, p in s1:
                for c2, q in s2:
                    tot += c1 * c2 * ystar4(*mom(p, q))
            P[i][j] = P[j][i] = Fraction(tot, 4)
        return P[i][j]
    fixed = [i for i in range(n) if pi[i] == i]
    pairs = sorted({tuple(sorted((i, pi[i]))) for i in range(n) if pi[i] != i})
    sym = [[(i, 1)] for i in fixed] + [[(i, 1), (j, 1)] for i, j in pairs]
    asym = [[(i, 1), (j, -1)] for i, j in pairs]
    dim_ok = (len(sym) + len(asym) == n)
    # utilde for symmetric combos
    ut = [sum(c * ker[i][1] for i, c in comb) for comb in sym]
    # antisym compression == 0
    a_ok = True
    for x in asym:
        for y in asym:
            v = sum(cx * cy * pval(ix, iy) for ix, cx in x for iy, cy in y)
            if v != 0:
                a_ok = False
    # sym compression == ut ut^T
    s_ok = True
    for ai, x in enumerate(sym):
        for bi_, y in enumerate(sym):
            v = sum(cx * cy * pval(ix, iy) for ix, cx in x for iy, cy in y)
            if v != Fraction(ut[ai] * ut[bi_]):
                s_ok = False
    # cross compression sym-antisym == 0 (needed for the block-diag claim)
    c_ok = True
    for x in sym:
        for y in asym:
            v = sum(cx * cy * pval(ix, iy) for ix, cx in x for iy, cy in y)
            if v != 0:
                c_ok = False
    # u.c = 0 automatically for antisym combos
    uc_ok = all(sum(c * ker[i][1] for i, c in comb) == 0 for comb in asym)
    res = all([lab_ok, u_ok, inv_ok, dim_ok, a_ok, s_ok, c_ok, uc_ok])
    ok = ok and res
    print(f"k={k}: perm-labels={lab_ok} u-invar={u_ok} involution={inv_ok} "
          f"dims({len(sym)}+{len(asym)}={n})={dim_ok} antisym0={a_ok} "
          f"sym=utut^T={s_ok} cross0={c_ok} u.c_asym=0={uc_ok}")
print("ALL PASS" if ok else "FAILURES")
sys.exit(0 if ok else 1)
