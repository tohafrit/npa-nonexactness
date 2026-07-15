"""Float-free verifier, UPPER-BOUND value for a_{1+AB} = 3/64.

Companion to verify_a1AB.py (which proves the lower bound a_{1+AB} >= 3/64 via a
primal second-order tangent-feasible arc).  Here we pin the dual side.

The optimal dual Z(s)=Z0+sZ1+s^2Z2 of the level-1+AB SDP is complementary to the
optimal primal M(y(s)) and satisfies the class-sum (dual-feasibility) identities
<Z(s),P_v> = -(B_s)_v.  To second order in s this is the LINEAR system

   (0)  Z0 M0 = 0                         [Z0 = N C0 N^T, supported on ker M0]
   (1)  Z0 M1 + Z1 M0 = 0
   (2)  Z0 M2 + Z1 M1 + Z2 M0 = 0
   (cs) <Z0,P_v> = -(B0)_v,  <Z1,P_v> = m_v,  <Z2,P_v> = 0     (all v)

with M0=M(y0), M1=L(y1), M2=L(y2) the exact Q(sqrt3) data of the optimizer's
2nd-order jet.  By SDP strong duality (Slater holds: the maximally-mixed point is
strictly feasible) the dual value <Z(s),1> equals c_{1+AB}(s), so the s^2
coefficient <Z2,1> equals a_{1+AB}.

THEOREM (verified here, symbolically, by field split Q(sqrt3)->Q):
   <Z2, 1_const> = 3/64  on the ENTIRE solution set of (0)-(cs), with NO dependence
   on any free parameter.
Hence the dual value's s^2 coefficient is pinned to 3/64, i.e. a_{1+AB} = 3/64
(together with the lower bound of verify_a1AB.py; a tangent-feasible PSD dual
representative exists, so the value is attained as a genuine dual bound).
Requires sympy.
"""
import sympy as sp

def main():
    r3 = sp.sqrt(3); R = sp.Rational
    sym = {0:'1',1:'a',3:'a',2:'b',4:'b',5:'c',8:'d',6:'e',7:'e',9:'f',12:'f',
           10:'g',13:'g',11:'h',14:'h',15:'p',16:'q'}
    idxtab=[[0,1,2,3,4,5,6,7,8],[1,0,9,5,6,3,4,10,11],[2,9,0,7,8,10,11,3,4],
            [3,5,7,0,12,1,13,2,14],[4,6,8,12,0,13,1,14,2],[5,3,10,1,13,0,12,9,15],
            [6,4,11,13,1,12,0,16,9],[7,10,3,2,14,9,16,0,12],[8,11,4,14,2,15,9,12,0]]
    V=['a','b','c','d','e','f','g','h','p','q']
    y0={'a':1,'b':R(1,2),'c':1,'d':0,'e':R(1,2),'f':R(1,2),'g':R(1,2),'h':0,'p':0,'q':0,'1':1}
    y1={'a':R(-3,32),'b':(-1+r3)/96,'c':0,'d':(-1+r3)/48,'e':(8+r3)/96,'f':(-4+r3)/96,
        'g':(-13+r3)/96,'h':(13+2*r3)/96,'p':(-13+r3)/48,'q':(23+r3)/48,'1':0}
    y2={'a':R(-75,512)+R(3,512)*r3,'b':R(-283,4608)+R(23,1152)*r3,'c':R(-1,32),
        'd':R(-283,2304)+R(23,576)*r3,'e':R(89,1152)+R(65,4608)*r3,'f':R(-91,1152)+R(89,4608)*r3,
        'g':R(-715,4608)+R(29,1152)*r3,'h':R(235,4608)+R(163,4608)*r3,
        'p':R(-643,2304)+R(29,576)*r3,'q':R(905,2304)+R(11,576)*r3,'1':0}
    def Mof(y): return sp.Matrix(9,9,lambda i,j: sp.sympify(y[sym[idxtab[i][j]]]))
    M0,M1,M2 = Mof(y0),Mof(y1),Mof(y2)
    N = sp.Matrix([[-1,1,0,0,0,0,0,0,0],[-1,0,0,1,0,0,0,0,0],[-1,0,0,0,0,1,0,0,0],
                   [0,0,0,0,-1,0,1,0,0],[0,0,-1,0,0,0,0,1,0],[1,0,-1,0,-1,0,0,0,1]]).T
    Pm = {v: sp.Matrix(9,9,lambda i,j: 1 if sym[idxtab[i][j]]==v else 0) for v in V}
    Mc = sp.Matrix(9,9,lambda i,j: 1 if sym[idxtab[i][j]]=='1' else 0)
    B0 = {'a':2,'b':0,'c':1,'d':-1,'e':2}; m = {'a':1}

    # field-split unknowns: each entry = p + q*sqrt3
    def mksym(pre,nn):
        d={}
        for i in range(nn):
            for j in range(i,nn):
                p=sp.Symbol(f'{pre}p{i}_{j}'); q=sp.Symbol(f'{pre}q{i}_{j}')
                d[(i,j)]=p+q*r3; d[(j,i)]=p+q*r3
        return d
    c0d=mksym('c',6); z1d=mksym('u',9); z2d=mksym('w',9)
    C0=sp.Matrix(6,6,lambda i,j:c0d[(i,j)]); Z1=sp.Matrix(9,9,lambda i,j:z1d[(i,j)]); Z2=sp.Matrix(9,9,lambda i,j:z2d[(i,j)])
    Z0=N*C0*N.T
    unk=sorted({s for M in (C0,Z1,Z2) for e in M for s in e.free_symbols}, key=str)

    eqs=[]
    def split_add(expr):
        e=sp.expand(expr); p=e.subs(r3,0); q=sp.expand((e-p)/r3)
        eqs.append(p); eqs.append(q)
    for v in V: split_add(sum(Z0[i,j]*Pm[v][i,j] for i in range(9) for j in range(9))+B0.get(v,0))
    E1=Z0*M1+Z1*M0
    for i in range(9):
        for j in range(9): split_add(E1[i,j])
    E2=Z0*M2+Z1*M1+Z2*M0
    for i in range(9):
        for j in range(9): split_add(E2[i,j])
    for v in V: split_add(sum(Z1[i,j]*Pm[v][i,j] for i in range(9) for j in range(9))-m.get(v,0))
    for v in V: split_add(sum(Z2[i,j]*Pm[v][i,j] for i in range(9) for j in range(9)))

    sol=sp.linsolve(eqs,unk)
    ok_solve=bool(sol)
    V2=None; free=None
    if ok_solve:
        s0=list(sol)[0]; sub=dict(zip(unk,s0))
        V2expr=sp.expand(sum(Z2[i,j]*Mc[i,j] for i in range(9) for j in range(9)).subs(sub))
        free=V2expr.free_symbols
        V2=sp.nsimplify(V2expr.subs({f:0 for f in free}))
    print(f"dual complementarity+class-sum system solvable: {ok_solve}")
    print(f"<Z2, 1_const> = {V2}")
    print(f"free-parameter dependence of <Z2,1>: {free if free else 'NONE (invariant)'}")
    ok = ok_solve and V2==R(3,64) and not free
    print("\na_{1+AB} = 3/64  (dual value's s^2 coefficient pinned) VERIFIED (float-free)" if ok
          else "\nVERIFICATION FAILED")
    return 0 if ok else 1

if __name__ == "__main__":
    raise SystemExit(main())
