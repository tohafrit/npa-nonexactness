"""Float-free verifier: the almost-quantum (level 1+AB) overshoot coefficient of
the symmetric doubly-tilted CHSH functional satisfies  a_{1+AB} >= 3/64,
hence level 1+AB strictly overshoots with curvature at least 3/64 > 0.

Set-up (party-swap-symmetric reduction; 10 moments a=<A0>=<B0>, b=<A1>=<B1>,
c=<A0B0>, d=<A1B1>, e=<A0B1>=<A1B0>, f=<A0A1>=<B0B1>, g=<A0A1B0>=<A0B0B1>,
h=<A0A1B1>=<A1B0B1>, p=<A0A1B0B1>, q=<A0A1B1B0>).  Objective B_s = (2-s)a+c-d+2e,
c_{1+AB}(s)=max{B_s : M(y)>=0}=4-s+a_{1+AB}s^2+O(s^3).

Proof of the lower bound (all exact in Q(sqrt 3)):
  y0 = (1,1/2,1,0,1/2,1/2,1/2,0,0,0)  is a maximizer at s=0 (value 4); M(y0)>=0,
       rank 3, with explicit integer 6-dim kernel N.
  y1, y2 in Q(sqrt3) (the s- and s^2-Taylor coefficients of the optimizer path).
  Checks:
   (0) M0:=M(y0) >= 0 and M0*N = 0                         [y0 optimal, N=ker M0]
   (1) B0.y1 = 0        (first-order objective coeff -1, arc stays optimal)
   (2) N^T M1 N >= 0    (M1=L(y1); first-order tangency to the PSD cone)
   (3) on ker(N^T M1 N): N2^T(M2 - M1 M0^+ M1)N2 >= 0  (M2=L(y2); the exact SDP
       second-order tangent-set condition, Bonnans-Shapiro)
   (4) B0.y2 - m.y1 = 3/64, with m.y1 = (y1)_a
  By (2)+(3) the second-order tangent-set condition holds, so a feasible arc
  y(s)=y0+s y1+s^2 y2+o(s^2) with M(y(s))>=0 exists; by (1)+(4) its objective is
  4 - s + (3/64)s^2 + o(s^2).  Hence c_{1+AB}(s) >= 4-s+(3/64)s^2+o(s^2), i.e.
  a_{1+AB} >= 3/64.  (Numerically a_{1+AB}=3/64 to 15 digits, so the bound is tight.)
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

    def psd(M):  # exact PSD test via eigenvalues (small matrices)
        return all(sp.re(sp.nsimplify(e)) >= 0 for e in M.eigenvals())

    ok = True
    c0 = psd(M0) and (M0*N).is_zero_matrix
    print(f"(0) M0>=0 and M0*N=0                : {c0}")
    B0y1 = 2*y1['a']+y1['c']-y1['d']+2*y1['e']
    c1 = sp.simplify(B0y1)==0
    print(f"(1) B0.y1 = 0                       : {c1}  (={sp.simplify(B0y1)})")
    G1 = sp.simplify(N.T*M1*N)
    c2 = psd(G1)
    print(f"(2) N^T M1 N >= 0                   : {c2}  eig={ {sp.simplify(k):v for k,v in G1.eigenvals().items()} }")
    N2 = N*sp.Matrix.hstack(*G1.nullspace())
    G2 = sp.simplify(N2.T*(M2 - M1*M0.pinv()*M1)*N2)
    c3 = psd(G2)
    print(f"(3) second-order tangent Schur >= 0 : {c3}  eig={ {sp.simplify(k):v for k,v in G2.eigenvals().items()} }")
    coeff = sp.simplify(2*y2['a']+y2['c']-y2['d']+2*y2['e'] - y1['a'])
    c4 = coeff == R(3,64)
    print(f"(4) B0.y2 - m.y1 = 3/64             : {c4}  (={coeff})")
    ok = c0 and c1 and c2 and c3 and c4
    print("\na_{1+AB} >= 3/64 > 0  VERIFIED (float-free, exact in Q(sqrt3))" if ok
          else "\nVERIFICATION FAILED")
    return 0 if ok else 1

if __name__ == "__main__":
    raise SystemExit(main())
