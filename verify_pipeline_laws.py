"""Machine checks for the structural laws of the paper
npa_nonexactness_proof.tex, Section 3: Theorem 3.1 (L1 face law, via
Lemma 3.2 dual certificate, Lemma 3.4 linear closure, Lemma 3.6 product
law), Theorem 3.7 (L2 kernel law), Theorem 3.9 (L4 uniform Slater), plus
the u_k facts of Definition 4.1 / Lemma 4.2 and Proposition 2.1.
(Development record: the level-uniform laws (L1)-(L4) of
pipeline_laws_note.md; FINDINGS_wt.md N+49..N+51b;
verify_v{2,3,4}_rational_base.py.)

Checks are grouped by law; everything here is EXACT (sympy rationals /
DomainMatrix over QQ), no floats in any verdict.  Levels k = 2,3,4 throughout,
k = 5 where it is a NEW data point (kernel counts, dim V_k, u_5).

Sections:
  [B] counting laws (paper Prop 2.1):  N(k) = 2k^2+2k+1 basis words,
      NC(k) = (5k+2)(k+1)/2 classes, kernel = (3k^2+k)/2 trivial
      reductions + k(k-1)/2 dressings = 2k^2, all independent.
  [C] L1 face law (paper Thm 3.1, Lemmas 3.2/3.4/3.6): rational max-rank
      dual Z0* = W E W^T (E > 0 explicit),
      identity 4 - B0.y = <Z0*, Gamma_{1+AB}(y)>, linear closure rank 15,
      level-k solution-set dim k with m_ij = delta_i + delta_j - 1.
  [D] L2 kernel law (paper Thm 3.7, incl. the base-point integrals of its
      defect realization): core Gram G_k(delta) formula, defect-realization
      moments (sin^2 measure) = base family (1/2,1/4,1/2,...,1/2), G_k PD
      exact, hence ker Gamma_k(base) = relation lattice exactly,
      ranks (k+1,k).
  [E] L3 first-order law (development record; the dim-V_k law is NOT used
      by the paper, see its Sec 9 "dimension law"; u_k = paper Def 4.1,
      identities = paper Lemma 4.2): dim V_k = 2k+3, u_k integer extension
      tower (u_5 = new data point), identities y1_<A0> = -lam/8,
      B0.y1 = 0 on V_k.
  [F] L4 uniform Slater (paper Thm 3.9): y2* = e_{class(1)} - y0(delta=0)
      has compressed pair ((NsKs)^T Ps0 (NsKs), Na^T Pa0 Na) exactly PD at
      every level.

Exit code 0 iff all checks pass.  Runtime: ~2 s.
"""
import sys
import time
import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix

sys.path.insert(0, '/Users/antonp/research/why-tsirelson')
from scenario_ext import ScenarioL          # noqa: E402
from npa_general import canon               # noqa: E402

R = sp.Rational
OK = True


def chk(name, cond):
    global OK
    OK = OK and bool(cond)
    print(f"  [{'OK' if cond else 'FAIL'}] {name}")


def drank(rows):
    """exact rank over QQ of a list-of-lists (ints/Rationals)."""
    if not rows:
        return 0
    def q(x):
        x = sp.Rational(x)
        return QQ(int(x.p), int(x.q))
    M = DomainMatrix([[q(x) for x in r] for r in rows],
                     (len(rows), len(rows[0])), QQ)
    return M.rank()


# ---------------------------------------------------------------- level data
class Level:
    """Exact class tables + canonical integer kernel for ScenarioL(k)."""

    def __init__(self, k):
        self.k = k
        s = ScenarioL(k)
        self.s = s
        self.N = s.N
        self.bi = {b: i for i, b in enumerate(s.basis)}
        inv = {v: kk for kk, v in s.words.items()}
        swp = [s.words[canon(bw, aw)] for i in range(s.NV)
               for aw, bw in [inv[i]]]
        cls = [-1]*s.NV
        nc = 0
        for i in range(s.NV):
            if cls[i] < 0:
                cls[i] = cls[swp[i]] = nc
                nc += 1
        self.cls, self.NC = cls, nc
        self.inv = inv          # word-index -> (aw,bw)
        # class representative word (first occurrence)
        self.crep = {}
        for i in range(s.NV):
            self.crep.setdefault(cls[i], inv[i])

    def W(self, a, b):
        """class of the moment word (a,b)."""
        return self.cls[self.s.words[canon(tuple(a), tuple(b))]]

    # canonical kernel: trivial trailing-0 reductions + R1 dressings
    @staticmethod
    def colreduce(aw, bw):
        aw, bw = list(aw), list(bw)
        ch = True
        while ch:
            ch = False
            if bw and bw[-1] == 0:
                bw.pop()
                ch = True
            if aw and aw[-1] == 0:
                aw.pop()
                ch = True
        return (tuple(aw), tuple(bw))

    def kernel_vectors(self):
        """returns (Ttriv list, Dress list) of length-N integer vectors."""
        N, bi = self.N, self.bi
        Ts = []
        for q, (aw, bw) in enumerate(self.s.basis):
            r = self.colreduce(aw, bw)
            if r != (aw, bw):
                v = [0]*N
                v[q] = 1
                v[bi[r]] = -1
                Ts.append(v)
        # dressing words: (wa, wb), wa/wb alternating NOT ending in 1
        # (i.e. empty or ending in 0), |wa|+|wb| <= k-2; one word per length.
        def endw(l):        # the unique alternating word of length l ending in 0
            return tuple(((l - 1 - i) % 2) for i in range(l))  # ends with 0
        Ds = []
        self.dress_words = []
        for la in range(0, self.k - 1):
            for lb in range(0, self.k - 1 - la):
                wa, wb = endw(la), endw(lb)
                v = [0]*N
                for da, db, sg in [((), (), 1), ((1,), (), -1),
                                   ((), (1,), -1), ((1,), (1,), 1)]:
                    a, b = self.colreduce(wa + da, wb + db)
                    v[bi[(a, b)]] += sg
                Ds.append(v)
                self.dress_words.append((wa, wb))
        return Ts, Ds


def check_endw():
    # sanity of endw: length-l alternating word ending in 0
    lv = Level(4)
    ws = [(), (0,), (1, 0), (0, 1, 0)]
    got = []
    for l in range(4):
        got.append(tuple(((l - 1 - i) % 2) for i in range(l)))
    return got == ws


# ============================================================ [B] counting
def section_B():
    print("(B) counting laws, k = 2..5:")
    allok = True
    for k in range(2, 6):
        lv = Level(k)
        Ts, Ds = lv.kernel_vectors()
        nT, nD = len(Ts), len(Ds)
        rk = drank(Ts + Ds)
        ok = (lv.N == 2*k*k + 2*k + 1
              and lv.NC == (5*k + 2)*(k + 1)//2
              and nT == (3*k*k + k)//2
              and nD == k*(k - 1)//2
              and rk == 2*k*k)
        allok = allok and ok
        print(f"    k={k}: N={lv.N} NC={lv.NC} triv={nT} dress={nD} "
              f"rank={rk}  {'OK' if ok else 'FAIL'}")
    chk("N=2k^2+2k+1, NC=(5k+2)(k+1)/2, triv=(3k^2+k)/2, dress=k(k-1)/2, "
        "kernel rank 2k^2 (k=2..5)", allok)
    chk("endw convention", check_endw())


# ============================================================ [C] L1 face law
# rational max-rank dual at level 1+AB:  Z0* = W E W^T,  E > 0 (5x5),
# 4 - B0.y = <Z0*, Gamma_{1+AB}(y)>  identically in y.
E_DUAL = [[R(13, 40), 0, R(3, 40), 0, 0],
          [0, R(1, 5), 0, R(1, 20), 0],
          [R(3, 40), 0, R(1, 8), 0, 0],
          [0, R(1, 20), 0, R(1, 8), 0],
          [0, 0, 0, 0, R(1, 4)]]


def onep_relations(lv):
    """the five 1+AB relation vectors W (columns) inside level-k basis."""
    bi = lv.bi
    N = lv.N

    def ev(*pl):
        v = sp.zeros(N, 1)
        for c, w in pl:
            v[bi[w]] += c
        return v
    return sp.Matrix.hstack(
        ev((1, ((0,), ())), (-1, ((), (0,)))),                    # a0 - b0
        ev((1, ((0,), (0,))), (-1, ((), ()))),                    # A0B0 - 1
        ev((1, ((0,), (1,))), (-1, ((), (1,))),
           (-1, ((1,), (0,))), (1, ((1,), ()))),                  # (A0-1)B1-(B0-1)A1
        ev((1, ((0,), ())), (1, ((), (0,))), (-2, ((), ())),
           (1, ((0,), (1,))), (-1, ((), (1,))),
           (1, ((1,), (0,))), (-1, ((1,), ()))),                  # w4
        ev((1, ((), ())), (-1, ((1,), ())),
           (-1, ((), (1,))), (1, ((1,), (1,)))))                  # R1


B0_WORDS = [(((0,), ()), 1), (((), (0,)), 1), (((0,), (0,)), 1),
            (((0,), (1,)), 1), (((1,), (0,)), 1), (((1,), (1,)), -1)]


def is_pd(M):
    return all(M[:m, :m].det() > 0 for m in range(1, M.rows + 1))


def face_value(lv, a, b, delta):
    """face parametrization: m_{i,j}, i = #1's in a, j = #1's in b,
    m_{0,0}=1, m_{0,j}=delta_j, m_{i,j}=delta_i+delta_j-1 (i,j>=1)."""
    i, j = sum(1 for x in a if x == 1), sum(1 for x in b if x == 1)
    if i == 0 and j == 0:
        return sp.Integer(1)
    if i == 0:
        return delta[j]
    if j == 0:
        return delta[i]
    return delta[i] + delta[j] - 1


def section_C():
    print("(C) L1 face law:")
    from npa_general import Scenario
    s1 = Scenario(2, 2)

    class LV1:                                    # minimal Level-like shim
        pass
    lv1 = LV1()
    lv1.bi = {b: i for i, b in enumerate(s1.basis)}
    lv1.N = s1.N
    Wm = onep_relations(lv1)
    E = sp.Matrix(E_DUAL)
    chk("E symmetric PD (exact leading minors), so Z0* = W E W^T is PSD "
        "with range = span W (rank 5)",
        E.T == E and is_pd(E) and Wm.rank() == 5)
    Z0 = Wm*E*Wm.T
    # affine identity: class sums of Z0 = 4 at identity class, -B0 coeff else
    b0w = {s1.words[w]: c for w, c in B0_WORDS}
    w0 = s1.words[((), ())]
    idok = True
    for w in range(s1.NV):
        ssum = sum(Z0[p, q] for p in range(9) for q in range(9)
                   if s1.idx[p, q] == w)
        tgt = 4 if w == w0 else -b0w.get(w, 0)
        idok = idok and (ssum == tgt)
    chk("identity 4 - B0.y == <Z0*, Gamma_{1+AB}(y)> for ALL y "
        "(class sums exact)", idok)
    # linear closure: Gamma.W = 0, y_1 = 1 forces the (betaA,betaB) table
    ys = [sp.Symbol(f'y{c}') for c in range(s1.NV)]
    G = sp.Matrix(9, 9, lambda p, q: ys[s1.idx[p, q]])
    GW = G*Wm
    eqs = [GW[i, j] for i in range(9) for j in range(5)]
    A_, _ = sp.linear_eq_to_matrix(eqs + [ys[w0] - 1], ys)
    chk("linear closure at 1+AB: rank 15 of 17 (free = betaA, betaB only)",
        A_.rank() == 15)
    sol = list(sp.linsolve((A_, sp.Matrix([0]*len(eqs) + [1])), ys))[0]
    bA, bB = sp.symbols('bA bB')
    s2 = sp.solve([sol[s1.words[((1,), ())]] - bA,
                   sol[s1.words[((), (1,))]] - bB],
                  sorted({f for e in sol for f in e.free_symbols}, key=str),
                  dict=True)
    Y = [sp.expand(e.subs(s2[0])) for e in sol]
    tab = {((0,), ()): 1, ((), (0,)): 1, ((0,), (0,)): 1,
           ((0,), (1,)): bB, ((1,), (0,)): bA,
           ((1,), (1,)): bA + bB - 1, ((0, 1), ()): bA,
           ((0, 1), (0, 1)): bA + bB - 1}
    chk("closure forces <A0>=<B0>=<A0B0>=1, <A0B1>=betaB, <A1B0>=betaA, "
        "<A1B1>=betaA+betaB-1, <A0A1>=betaA, <A0A1B0B1>=betaA+betaB-1",
        all(sp.expand(Y[s1.words[w]] - v) == 0 for w, v in tab.items()))

    # level-k affine hull: dim k, parametrization m_ij
    for k in range(2, 5):
        lv = Level(k)
        Ts, Ds = lv.kernel_vectors()
        Nk = sp.Matrix(Ts + Ds).T
        delta = [None] + list(sp.symbols(f'd1:{k+1}'))
        Gam = sp.Matrix(lv.N, lv.N, lambda p, q: 0)
        for p in range(lv.N):
            ap, bp = lv.s.basis[p]
            for q in range(lv.N):
                aq, bq = lv.s.basis[q]
                a, b = lv.inv[lv.s.idx[p, q]]
                Gam[p, q] = face_value(lv, a, b, delta)
        chk(f"k={k}: Gamma(delta).N_k == 0 identically "
            "(m_ij = d_i + d_j - 1 solves the kernel system)",
            sp.expand(Gam*Nk) == sp.zeros(lv.N, Nk.cols))
        # affine hull dim exactly k
        ysym = [sp.Symbol(f'q{c}') for c in range(lv.NC)]
        rows = {}
        for p in range(lv.N):
            for jv, v in enumerate(Ts + Ds):
                key = (p, jv)
                r = [0]*lv.NC
                nz = False
                for q, cq in enumerate(v):
                    if cq:
                        r[lv.cls[lv.s.idx[p, q]]] += cq
                        nz = True
                if nz and any(r):
                    rows[tuple(r)] = 1
        rows = [list(r) for r in rows]
        rk = drank(rows)
        chk(f"k={k}: affine hull of the face has dim EXACTLY k "
            f"(homog rank {rk} = NC-(k+1) = {lv.NC-k-1})",
            rk == lv.NC - k - 1)
        # B0.y == 4 identically on the face
        b0v = sp.Integer(0)
        for w, c in B0_WORDS:
            b0v += c*face_value(lv, w[0], w[1], delta)
        chk(f"k={k}: B0.y == 4 identically on the face",
            sp.expand(b0v - 4) == 0)


# ============================================================ [D] L2 kernel
def residual_words(k):
    """the 2k+1 residual basis words: psi, a_l, b_l (l = 1..k),
    a_l = the unique alternating a-word of length l ending in 1."""
    def endw1(l):
        return tuple(((l - i) % 2) for i in range(l))   # ends with 1
    out = [((), ())]
    for l in range(1, k + 1):
        out.append((endw1(l), ()))
    for l in range(1, k + 1):
        out.append(((), endw1(l)))
    return out


def base_delta(k):
    """the canonical rational base: (1/2, 1/4, 1/2, ..., 1/2)."""
    return [None, R(1, 2), R(1, 4)] + [R(1, 2)]*(k - 2)


def section_D():
    print("(D) L2 kernel law (core Gram + defect realization):")
    # (i) the defect realization reproduces the base family:
    #     delta_j = 1 - p(1 - c_j), p = 1/2, c_j = int cos(j phi) dnu,
    #     nu = sin^2(phi) dphi / pi on [0, 2pi).
    phi = sp.Symbol('phi')
    cs = [sp.integrate(sp.cos(j*phi)*sp.sin(phi)**2, (phi, 0, 2*sp.pi))/sp.pi
          for j in range(1, 7)]
    dl = [1 - R(1, 2)*(1 - c) for c in cs]
    chk("defect realization (p=1/2, nu = sin^2 measure) gives delta = "
        "(1/2, 1/4, 1/2, 1/2, ...) for all j (checked j <= 6)",
        dl[0] == R(1, 2) and dl[1] == R(1, 4)
        and all(d == R(1, 2) for d in dl[2:]))
    for k in range(2, 6):
        lv = Level(k)
        Ts, Ds = lv.kernel_vectors()
        res = residual_words(k)
        chk(f"k={k}: residual words all in basis, and "
            "[residual columns | kernel] spans R^N",
            all(w in lv.bi for w in res)
            and drank([[1 if i == lv.bi[w] else 0 for i in range(lv.N)]
                       for w in res] + Ts + Ds) == lv.N)
        delta = base_delta(k)
        G = sp.Matrix(2*k + 1, 2*k + 1, lambda i, j: 0)
        for i, wi in enumerate(res):
            for j, wj in enumerate(res):
                a, b = lv.inv[lv.s.idx[lv.bi[wi], lv.bi[wj]]]
                G[i, j] = face_value(lv, a, b, delta)
        chk(f"k={k}: core Gram G_k(base) is PD (exact leading minors) "
            "=> ker Gamma_k(base) = relation lattice EXACTLY (dim 2k^2)",
            is_pd(G))
        # sym/antisym ranks (k+1, k): psi & a_l+b_l  /  a_l-b_l
        S = sp.zeros(2*k + 1, k + 1)
        A = sp.zeros(2*k + 1, k)
        S[0, 0] = 1
        for l in range(k):
            S[1 + l, 1 + l] = 1
            S[1 + k + l, 1 + l] = 1
            A[1 + l, l] = 1
            A[1 + k + l, l] = -1
        chk(f"k={k}: sym/antisym base ranks (k+1, k) = "
            f"({k+1}, {k}) (both compressed Grams PD)",
            is_pd(S.T*G*S) and is_pd(A.T*G*A))
        chk(f"k={k}: closed-form entry rules reproduce G_k from the tables",
            core_gram_closed_form(k, delta) == G)
    # larger k: PD of G_k(base) directly from the closed-form entries
    bigok = True
    for k in range(6, 9):
        lvK = None                          # no scenario needed: use rules
        G = core_gram_closed_form(k, base_delta(k))
        bigok = bigok and is_pd(G)
    chk("k=6..8: G_k(base) PD via closed-form entry rules (no scenario "
        "tables)", bigok)


def core_gram_closed_form(k, delta):
    """G_k from the parity entry rules (note, Lemma L2.2):
    [psi,psi]=1, [psi,a_l]=delta_{ceil(l/2)}, a-a same parity: delta_{|l-m|/2},
    a-a opposite: delta_{(l+m+1)/2}, a-b: delta_{ceil(l/2)}+delta_{ceil(m/2)}-1
    (delta_0 := 1)."""
    def d(j):
        return sp.Integer(1) if j == 0 else delta[j]
    n = 2*k + 1
    G = sp.Matrix(n, n, lambda i, j: 0)

    def idx(i):
        # 0 -> psi; 1..k -> a_l; k+1..2k -> b_l
        if i == 0:
            return ('p', 0)
        if i <= k:
            return ('a', i)
        return ('b', i - k)
    for i in range(n):
        ti, li = idx(i)
        for j in range(n):
            tj, lj = idx(j)
            if ti == 'p' and tj == 'p':
                G[i, j] = 1
            elif ti == 'p' or tj == 'p':
                l = li if tj == 'p' else lj
                G[i, j] = d(-(-l//2))
            elif ti == tj:
                l, m = li, lj
                if (l - m) % 2 == 0:
                    G[i, j] = d(abs(l - m)//2)
                else:
                    G[i, j] = d((l + m + 1)//2)
            else:
                G[i, j] = d(-(-li//2)) + d(-(-lj//2)) - 1
    return G


# ============================================================ [E] L3 first order
def u_closed_form(lv):
    """THE CLOSED FORM of the unlock direction u_k (canonical kernel basis
    [trivial reductions in basis order | dressings in (la,lb) order]):
      u[T_(a,b)] = s(a) s(b)   if exactly one of a,b ends in 0, else 0,
                   where s(w) = -1 if w starts with 1, else +1;
      u[D_W]     = 2 (-1)^{|W|}  for W != (), and 0 for W = () (= R1).
    Extends u_{k-1} exactly (the closed form is level-independent)."""
    def s(w):
        return -1 if (w and w[0] == 1) else 1

    def ends0(w):
        return bool(w) and w[-1] == 0
    u = []
    for q, (aw, bw) in enumerate(lv.s.basis):
        if lv.colreduce(aw, bw) != (aw, bw):
            u.append(s(aw)*s(bw) if ends0(aw) != ends0(bw) else 0)
    for (wa, wb) in lv.dress_words:
        L = len(wa) + len(wb)
        u.append(0 if L == 0 else 2*(-1)**L)
    return sp.Matrix(u)


def pairing_rows(lv, K):
    """rows of the map y1 -> P(y1)[a,b] = K_a^T Gamma(y1) K_b, a <= b."""
    nK, N, NC = len(K), lv.N, lv.NC
    rows = []
    for a in range(nK):
        for b in range(a, nK):
            r = [sp.Integer(0)]*NC
            va, vb = K[a], K[b]
            for p in range(N):
                if va[p] == 0:
                    continue
                for q in range(N):
                    if vb[q] == 0:
                        continue
                    r[lv.cls[lv.s.idx[p, q]]] += va[p]*vb[q]
            rows.append(r)
    return rows


def dm(rows, ncols):
    def q(x):
        x = sp.Rational(x)
        return QQ(int(x.p), int(x.q))
    return DomainMatrix([[q(x) for x in r] for r in rows],
                        (len(rows), ncols), QQ)


# certified lam!=0 jet at k=2 (verify_v2_rational_base.py Bv, last row)
Y1_CERT_K2 = [-R(1, 8), -R(7, 8), -R(3, 8), 0, -R(3, 4), -R(7, 4), -R(1, 8),
              -R(1, 2), -1, 0, 0, 0, 0, 0, -1, 0, 0]     # classes 1..17


def section_E():
    print("(E) L3 first-order law (V_k, closed-form u_k, tower, k=5 test):")
    lvs = {}
    for k in range(2, 6):
        lv = Level(k)
        Ts, Ds = lv.kernel_vectors()
        K = Ts + Ds
        lvs[k] = (lv, K)
        u = u_closed_form(lv)
        rows = pairing_rows(lv, K)
        nK, NC = len(K), lv.NC
        # V_k-hat := {(y1, lam): P(y1) = lam u u^T}  (y1_empty free here;
        # the verifiers' V_k = the y1_empty = 0 slice + lam, dim 2k+3)
        aug = []
        i = 0
        for a in range(nK):
            for b in range(a, nK):
                aug.append(rows[i] + [-u[a]*u[b]])
                i += 1
        M = dm(aug, NC + 1)
        ns = M.nullspace().to_Matrix()
        lam_ok = any(ns[i, NC] != 0 for i in range(ns.rows))
        chk(f"k={k}: dim{{(y1,lam): N^T Gamma(y1) N = lam u_k u_k^T}} "
            f"= 2k+4 = {2*k+4} (verifier convention 2k+3) with lam != 0 "
            f"achievable" + (" [NEW DATA POINT]" if k == 5 else ""),
            ns.rows == 2*k + 4 and lam_ok)
        # identities on the y_empty = 0 slice: y1_<A0> = -lam/2 (canonical
        # normalization of u; = -lam/8 at the verifiers' u-scaling),
        # B0.y1 = 0, y1_<A0B0> = 0
        ca = lv.cls[lv.s.words[((0,), ())]]
        cb0 = lv.cls[lv.s.words[((0,), (0,))]]
        cb1 = lv.cls[lv.s.words[((0,), (1,))]]
        cab = lv.cls[lv.s.words[((1,), (1,))]]
        dl = base_delta(k)
        idok = True
        for i in range(ns.rows):
            y = [sp.Rational(ns[i, c]) for c in range(NC)]
            lam = sp.Rational(ns[i, NC])
            y0c = y[0]
            yv = {c: y[c] - y0c*face_value(lv, *lv.crep[c], dl)
                  for c in [ca, cb0, cb1, cab]}
            idok = (idok and yv[ca] == -lam/2 and yv[cb0] == 0
                    and 2*yv[ca] + yv[cb0] + 2*yv[cb1] - yv[cab] == 0)
        chk(f"k={k}: identities on the slice y1_1=0:  y1_<A0> == -lam/2, "
            "y1_<A0B0> == 0, B0.y1 == 0", idok)
        # ker Phi structure: face span (k+1) + length-flats (k+1) + 4 unseen
        prof = [(len(lv.crep[c][0]), len(lv.crep[c][1])) for c in range(NC)]
        seen = set()
        for r in rows:
            for c in range(NC):
                if r[c] != 0:
                    seen.add(c)
        unseen = [c for c in range(NC) if c not in seen]
        maxb = tuple((i % 2) for i in range(1, 2*k))     # (1,0,...,1) len 2k-1
        chk(f"k={k}: exactly 4 unseen classes = the maximal beta_k-core "
            "words", len(unseen) == 4 and
            set(lv.crep[c] for c in unseen) ==
            {((), maxb), ((), (0,) + maxb), (((0,), maxb)), (((1,), maxb))})
        span = []
        for j in range(k + 1):        # face points at delta = e_j vertices
            dlj = [None] + [sp.Integer(1 if i == j else 0)
                            for i in range(1, k + 1)]
            span.append([face_value(lv, *lv.crep[c], dlj)
                         for c in range(NC)])
        # length-flat family: f0 = f1, f_{2i-1} - 2 f_{2i} + f_{2i+1} = 0
        # for i = 1..k-1  ->  explicit basis by free values at odd kinks
        fcond = [[0]*(2*k + 1) for _ in range(k)]
        fcond[0][0], fcond[0][1] = 1, -1
        for i in range(1, k):
            fcond[i][2*i - 1], fcond[i][2*i], fcond[i][2*i + 1] = 1, -2, 1
        fbasis = dm(fcond, 2*k + 1).nullspace().to_Matrix()
        for i in range(fbasis.rows):
            f = [fbasis[i, l] for l in range(2*k + 1)]
            span.append([f[p[0]] + f[p[1]] for p in prof])
        for c in unseen:
            span.append([1 if cc == c else 0 for cc in range(NC)])
        chk(f"k={k}: ker Phi = face-span + length-flats + unseen singles: "
            f"joint rank {2*k+3} = dim ker Phi exactly",
            dm(span, NC).rank() == 2*k + 3
            and dm(rows, NC).rank() == NC - (2*k + 3))
        # length-flats really are flat (each basis y_f solves P(y_f) = 0)
        flok = True
        for i in range(fbasis.rows):
            f = [fbasis[i, l] for l in range(2*k + 1)]
            yf = [f[p[0]] + f[p[1]] for p in prof]
            flok = flok and all(
                sum(r[c]*yf[c] for c in range(NC) if r[c]) == 0
                for r in rows)
        chk(f"k={k}: length-flat law {{f0=f1, f_(2i-1)-2f_(2i)+f_(2i+1)=0, "
            "i<=k-1}} => y_f in ker Phi (all k+1 basis f's)", flok)

    # closed-form u reproduces the CERTIFIED level-2 jet's rank-1 form
    lv, K = lvs[2]
    y1 = [sp.Integer(0)] + Y1_CERT_K2
    nK = len(K)
    P = sp.zeros(nK, nK)
    for a in range(nK):
        for b in range(a, nK):
            va, vb = K[a], K[b]
            val = 0
            for p in range(lv.N):
                if va[p] == 0:
                    continue
                for q in range(lv.N):
                    if vb[q] == 0:
                        continue
                    val += va[p]*vb[q]*y1[lv.cls[lv.s.idx[p, q]]]
            P[a, b] = P[b, a] = val
    u2 = u_closed_form(lv)
    chk("closed-form u_2 matches the certified verify_v2 jet: "
        "P(y1_cert) == (1/4) u_2 u_2^T exactly",
        P == u2*u2.T/4)
    # tower: u_k restricted to the embedded level-(k-1) kernel = u_{k-1}
    towok = True
    for k in [3, 4, 5]:
        lvp, Kp = lvs[k - 1]
        lv, K = lvs[k]
        up, u = u_closed_form(lvp), u_closed_form(lv)
        labp = ([('T', lvp.s.basis[q]) for q, (aw, bw)
                 in enumerate(lvp.s.basis)
                 if lvp.colreduce(aw, bw) != (aw, bw)]
                + [('D', w) for w in lvp.dress_words])
        lab = ([('T', lv.s.basis[q]) for q, (aw, bw)
                in enumerate(lv.s.basis)
                if lv.colreduce(aw, bw) != (aw, bw)]
               + [('D', w) for w in lv.dress_words])
        old = [lab.index(x) for x in labp]
        towok = towok and all(u[old[i]] == up[i] for i in range(len(old)))
    chk("tower u_5|old = u_4|old = u_3|old = u_2 exactly (closed form is "
        "level-independent)", towok)
    # gauge line: r = e_{D_()} + dressing-part(u); (u+r)(u+r)^T still
    # admits a jet (t = 1 member of the projective line), k = 2..4
    lineok = True
    for k in [2, 3, 4]:
        lv, K = lvs[k]
        u = u_closed_form(lv)
        nK, NC = len(K), lv.NC
        nT = sum(1 for q, (aw, bw) in enumerate(lv.s.basis)
                 if lv.colreduce(aw, bw) != (aw, bw))
        r = sp.zeros(nK, 1)
        for i in range(nT, nK):
            r[i] = u[i]/2
        r[nT + lv.dress_words.index(((), ()))] += 1
        v = u + r
        rows = pairing_rows(lv, K)
        aug = []
        i = 0
        for a in range(nK):
            for b in range(a, nK):
                aug.append(rows[i] + [-v[a]*v[b]])
                i += 1
        ns = dm(aug, NC + 1).nullspace().to_Matrix()
        lineok = lineok and any(ns[i, NC] != 0 for i in range(ns.rows))
    chk("gauge line: u + r (r = e_R1 + dressing-part(u)/2, canonical "
        "normalization) also admits a rank-1 unlock jet (t=1), k=2..4",
        lineok)


# ============================================================ [F] L4 Slater
def section_F():
    print("(F) L4 uniform rational Slater: y2* = e_1 - y0(delta=0):")
    for k in range(2, 6):
        lv = Level(k)
        Ts, Ds = lv.kernel_vectors()
        K = Ts + Ds
        nK, N, NC = len(K), lv.N, lv.NC
        zero = [None] + [sp.Integer(0)]*k
        y2 = [(1 if c == 0 else 0) - face_value(lv, *lv.crep[c], zero)
              for c in range(NC)]
        # pins: y2_1 = 0 and y2 at the k face-coordinate classes = 0
        cpins = [0] + [lv.cls[lv.s.words[canon((), tuple((i % 2) for i in
                       range(1, 2*j)))]] for j in range(1, k + 1)]
        chk(f"k={k}: y2* satisfies ALL pins (identity class + k face "
            "coordinates)", all(y2[c] == 0 for c in cpins))
        # compressed second-order form on the kernel = K^T Gamma(y2*) K
        # = K^T K  (since Gamma(e_1) = I and Gamma(y0(0)) K = 0)
        C = sp.zeros(nK, nK)
        for a in range(nK):
            for b in range(a, nK):
                va, vb = K[a], K[b]
                val = 0
                for p in range(N):
                    if va[p] == 0:
                        continue
                    for q in range(N):
                        if vb[q] == 0:
                            continue
                        val += va[p]*vb[q]*y2[lv.cls[lv.s.idx[p, q]]]
                C[a, b] = C[b, a] = val
        Km = sp.Matrix(K).T
        chk(f"k={k}: K^T Gamma(y2*) K == K^T K exactly (Gram of the 2k^2 "
            "independent kernel vectors => PD on the WHOLE kernel, hence "
            "on u-perp x antisym)", C == Km.T*Km)
    print("    (PD of K^T K is automatic: kernel vectors independent, "
          "sec. B)")


if __name__ == "__main__":
    t0 = time.time()
    section_B()
    section_C()
    section_D()
    section_E()
    section_F()
    print(f"\n[{time.time()-t0:.1f}s] " +
          ("ALL CHECKS PASS" if OK else "CHECK FAILURES"))
    raise SystemExit(0 if OK else 1)
