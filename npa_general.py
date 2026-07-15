"""NPA level 1+AB (= almost-quantum set Q~) for general (mA, mB) binary-outcome
Bell scenarios, with marginals.

Monomial basis: {1} u {A_i} u {B_j} u {A_i B_j}   (1 + mA + mB + mA*mB).
Algebra: A_i^2 = B_j^2 = 1, [A_i, B_j] = 0; A's don't commute among themselves.
Real (symmetrized) picture: <w> identified with <w^dagger> — Re of any complex PSD
solution satisfies these constraints and any real solution is a complex one, so the
projection onto (marginals, correlators) is exactly Q~'s.

API:
  scen = Scenario(mA, mB)
  val, sol = scen.max_linear(cA, cB, cAB)     # max sum cA_i<A_i>+cB_j<B_j>+cAB_ij<A_iB_j>
  feas     = scen.membership(a, b, c)          # behavior in Q~?
Calibrations in __main__ (must pass): CHSH, I_2222, I_3322.
"""
import numpy as np
import cvxpy as cp

def reduce_word(w):
    out = []
    for x in w:
        if out and out[-1] == x:
            out.pop()
        else:
            out.append(x)
    return tuple(out)

def dagger(aw, bw):
    return (tuple(reversed(aw)), tuple(reversed(bw)))

def canon(aw, bw):
    return min((aw, bw), dagger(aw, bw))

class Scenario:
    def __init__(self, mA, mB):
        self.mA, self.mB = mA, mB
        basis = [((), ())]
        basis += [((i,), ()) for i in range(mA)]
        basis += [((), (j,)) for j in range(mB)]
        basis += [((i,), (j,)) for i in range(mA) for j in range(mB)]
        self.basis = basis
        N = len(basis)
        self.N = N
        words = {}
        idx = np.zeros((N, N), dtype=int)
        for p in range(N):
            ap, bp = basis[p]
            for q in range(N):
                aq, bq = basis[q]
                w = canon(reduce_word(tuple(reversed(ap)) + aq),
                          reduce_word(tuple(reversed(bp)) + bq))
                if w not in words:
                    words[w] = len(words)
                idx[p, q] = words[w]
        self.words, self.idx = words, idx
        self.NV = len(words)

    def _moment(self, y):
        return cp.vstack([cp.hstack([y[self.idx[p, q]] for q in range(self.N)])
                          for p in range(self.N)])

    def _w(self, y, aw, bw):
        return y[self.words[canon(aw, bw)]]

    def _base(self):
        y = cp.Variable(self.NV)
        cons = [self._moment(y) >> 0, self._w(y, (), ()) == 1]
        return y, cons

    def _lin(self, y, cA, cB, cAB):
        terms = []
        for i in range(self.mA):
            if cA is not None and cA[i]:
                terms.append(cA[i] * self._w(y, (i,), ()))
        for j in range(self.mB):
            if cB is not None and cB[j]:
                terms.append(cB[j] * self._w(y, (), (j,)))
        for i in range(self.mA):
            for j in range(self.mB):
                if cAB[i][j]:
                    terms.append(cAB[i][j] * self._w(y, (i,), (j,)))
        return sum(terms)

    def max_linear(self, cA, cB, cAB, solver=cp.CLARABEL, **kw):
        y, cons = self._base()
        prob = cp.Problem(cp.Maximize(self._lin(y, cA, cB, cAB)), cons)
        prob.solve(solver=solver, **kw)
        return prob.value, y.value

    def membership(self, a, b, c, solver=cp.CLARABEL):
        y, cons = self._base()
        for i in range(self.mA):
            cons.append(self._w(y, (i,), ()) == a[i])
        for j in range(self.mB):
            cons.append(self._w(y, (), (j,)) == b[j])
        for i in range(self.mA):
            for j in range(self.mB):
                cons.append(self._w(y, (i,), (j,)) == c[i][j])
        prob = cp.Problem(cp.Minimize(0), cons)
        prob.solve(solver=solver)
        return prob.status in ("optimal", "optimal_inaccurate")

def cg_Innn22_coeffs(n):
    """Collins-Gisin I_nn22 in probability form -> correlator form.
    Prob form: sum_i (-(n-i)) P_A(1|i) [i=1..n]  - P_B(1|1)
               + sum_ij M_ij P(11|ij),  M_ij = +1 (i+j<=n+1), -1 (i+j=n+2), 0 else.
    LHV bound 0; NS max (n-1)/2; n=2 -> (S_chsh-2)/4.
    Returns (const, cA, cB, cAB) with I = const + sum cA<A> + sum cB<B> + sum cAB<AB>,
    P(11|ij)=(1+<Ai>+<Bj>+<AiBj>)/4, P_A(1|i)=(1+<Ai>)/2."""
    M = np.zeros((n, n))
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            s = i + j
            M[i - 1, j - 1] = 1.0 if s <= n + 1 else (-1.0 if s == n + 2 else 0.0)
    alphaA = np.array([-(n - i) for i in range(1, n + 1)], dtype=float)
    alphaB = np.zeros(n); alphaB[0] = -1.0
    const = alphaA.sum() / 2 + alphaB.sum() / 2 + M.sum() / 4
    cA = alphaA / 2 + M.sum(axis=1) / 4
    cB = alphaB / 2 + M.sum(axis=0) / 4
    cAB = M / 4
    return const, cA, cB, cAB

def max_Innn22(n, **kw):
    scen = Scenario(n, n)
    const, cA, cB, cAB = cg_Innn22_coeffs(n)
    v, _ = scen.max_linear(cA, cB, cAB, **kw)
    return const + v

if __name__ == "__main__":
    ok = True
    # CHSH via general engine
    s = Scenario(2, 2)
    v, _ = s.max_linear(None, None, [[1, 1], [1, -1]])
    t = 2 * np.sqrt(2)
    print(f"CHSH  : {v:.10f}  target {t:.10f}  diff {v-t:+.2e}")
    ok &= abs(v - t) < 1e-7
    # I_2222 = (2sqrt2-2)/4
    v2 = max_Innn22(2)
    t2 = (2 * np.sqrt(2) - 2) / 4
    print(f"I_2222: {v2:.10f}  target {t2:.10f}  diff {v2-t2:+.2e}")
    ok &= abs(v2 - t2) < 1e-7
    # I_3322 at level 1+AB: NPA 2008 (NJP 10, 073013) Table 3: I^{1+AB} = 0.2515
    # (level 2 = 0.25091, level 3 = 0.25089). Anchor: must round to 0.2515.
    v3 = max_Innn22(3)
    print(f"I_3322: {v3:.10f}  published 1+AB anchor 0.2515 (NPA'08 Table 3)")
    ok &= 0.25145 <= v3 <= 0.25155
    print("CALIBRATION:", "PASS" if ok else "FAIL")
    raise SystemExit(0 if ok else 1)
