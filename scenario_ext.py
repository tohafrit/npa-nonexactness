"""NPA levels 2 and 3 for the (2,2,2) scenario (binary settings/outcomes,
marginals included), reusing npa_general's word algebra.

Basis = all monomials (a-word, b-word) with alternating words and total length
<= level. Level 1+AB = Scenario(2,2) from npa_general (= almost-quantum Q~).
Ray membership: max t s.t. t*u in set, u an 8-vector
(mA0, mA1, mB0, mB1, C00, C01, C10, C11). One SDP per (u, level).
CONTROLS:
  - pure-correlator u: t_{Q~} == t_{L2} == t_{L3} (sandwich lemma; gap must be 0)
  - CHSH direction: correlator body radius known analytically:
    u = (0,0,0,0, 1,1,1,-1)/norm ... t*4-correlator sum hits 2sqrt2 boundary.
"""
import numpy as np
import cvxpy as cp
from npa_general import Scenario, canon, reduce_word

def alt_words(maxlen):
    """alternating words over {0,1} up to length maxlen (incl. empty)."""
    out = [()]
    for L in range(1, maxlen + 1):
        for start in (0, 1):
            w = tuple((start + i) % 2 for i in range(L))
            out.append(w)
    return out

class ScenarioL(Scenario):
    """NPA level-L (L=2 or 3) for (2,2) binary settings."""
    def __init__(self, level):
        self.mA, self.mB = 2, 2
        aws = alt_words(level)
        basis = []
        for aw in aws:
            for bw in alt_words(level - len(aw)):
                basis.append((aw, bw))
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

def ray_max_t(scen, u, solver=cp.CLARABEL):
    """max t such that behavior t*u is in the set. u = 8-vector."""
    y = cp.Variable(scen.NV)
    t = cp.Variable()
    cons = [scen._moment(y) >> 0, scen._w(y, (), ()) == 1]
    mA = u[0:2]; mB = u[2:4]; C = np.array(u[4:8]).reshape(2, 2)
    for x in range(2):
        cons.append(scen._w(y, (x,), ()) == t * mA[x])
        cons.append(scen._w(y, (), (x,)) == t * mB[x])
    for x in range(2):
        for yy in range(2):
            cons.append(scen._w(y, (x,), (yy,)) == t * C[x, yy])
    prob = cp.Problem(cp.Maximize(t), cons)
    prob.solve(solver=solver)
    return t.value

if __name__ == "__main__":
    s1 = Scenario(2, 2)      # 1+AB = Q~
    s2 = ScenarioL(2)
    s3 = ScenarioL(3)
    print(f"basis sizes: 1+AB {s1.N}, L2 {s2.N}, L3 {s3.N}; words {s1.NV}/{s2.NV}/{s3.NV}")
    # control 1: CHSH correlator direction — all levels must give same t, and
    # t*(sum of correlators pattern) at boundary: correlators t*(1,1,1,-1);
    # CHSH = 4t = 2sqrt2 => t = 1/sqrt2.
    u = np.array([0, 0, 0, 0, 1, 1, 1, -1], float)
    ts = [ray_max_t(s, u) for s in (s1, s2, s3)]
    print("CHSH dir t:", [f"{t:.8f}" for t in ts], f" target {1/np.sqrt(2):.8f}")
    # control 2: random correlator-only direction — gaps must be 0
    rng = np.random.default_rng(3)
    for _ in range(3):
        v = rng.normal(size=4)
        u = np.concatenate([[0, 0, 0, 0], v / np.linalg.norm(v)])
        ts = [ray_max_t(s, u) for s in (s1, s2, s3)]
        print("corr-only t:", [f"{t:.8f}" for t in ts],
              f" gap12 {ts[0]-ts[1]:+.2e} gap13 {ts[0]-ts[2]:+.2e}")
    # control 3: slice-3-like direction (biased marginals)
    u = np.array([0.3, 0.3, 0.3, 0.3, 0.9, 0.9, 0.9, 0.1], float)
    u /= np.linalg.norm(u)
    ts = [ray_max_t(s, u) for s in (s1, s2, s3)]
    print("biased dir t:", [f"{t:.8f}" for t in ts],
          f" gap13 {ts[0]-ts[2]:+.3e}")
