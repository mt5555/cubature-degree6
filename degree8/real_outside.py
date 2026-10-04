# For each real root of f: rational isolating interval for t, interval enclosure of y(t) = -n_i(t)/(c_i q(t)),
# then an exact certificate that some localizing LMI fails: v^T G_g(y) v < 0 for every y in the enclosure.
import sympy as sp, numpy as np
from fractions import Fraction as Fr
from parse_rur import load
from cf8 import K, m8, geigs, G
t=sp.symbols('t'); R=load(); P=R[1][5][1]
poly=lambda d,cs: sp.Poly(list(reversed(cs)),t,domain='ZZ')
f=poly(*P[0]); q=poly(*P[1]); nums=[(poly(*a[0]),int(a[1])) for a in P[2]]
class Iv:
    def __init__(s,lo,hi=None): s.lo=Fr(lo); s.hi=Fr(lo if hi is None else hi)
    def __add__(s,o): o=o if isinstance(o,Iv) else Iv(o); return Iv(s.lo+o.lo,s.hi+o.hi)
    __radd__=__add__
    def __mul__(s,o):
        o=o if isinstance(o,Iv) else Iv(o); p=[s.lo*o.lo,s.lo*o.hi,s.hi*o.lo,s.hi*o.hi]; return Iv(min(p),max(p))
    __rmul__=__mul__
    def inv(s):
        assert s.lo>0 or s.hi<0; return Iv(1/s.hi,1/s.lo)
def peval(p,T):
    acc=Iv(0)
    for c in p.all_coeffs(): acc=acc*T+Iv(int(c))
    return acc
ivs=sp.Poly(f.as_expr(),t).intervals(eps=Fr(1,10**60))
assert len(ivs)==2
fr=lambda x: Fr(int(sp.fraction(x)[0]),int(sp.fraction(x)[1]))
Kq={g:[[fr(K[g][i,j]) for j in range(5)] for i in range(5)] for g in K}
m8q=[fr(v) for v in m8]
def Gq(g,Y):
    M=[[None]*5 for _ in range(5)]
    for i in range(5):
        for j in range(5):
            k=i+j
            if g=='x': d=Y[k]
            elif g=='y': d=Y[k+1]
            else: d=Iv(m8q[k])+Y[k]*(-1)+Y[k+1]*(-1)
            M[i][j]=d+Iv(-Kq[g][i][j])
    return M
for n,((a,b),mult) in enumerate(ivs):
    T=Iv(a,b); qi=peval(q,T).inv()
    Y=[peval(nm,T)*Iv(Fr(-1,c))*qi for nm,c in nums]
    ymid=np.array([float((v.lo+v.hi)/2) for v in Y]); width=max(float(v.hi-v.lo) for v in Y)
    e=geigs(ymid); g=min(e,key=lambda k:e[k].min())
    w,V=np.linalg.eigh(G(ymid)[g]); v=[Fr(float(x)).limit_denominator(10**8) for x in V[:,0]]
    M=Gq(g,Y); val=Iv(0)
    for i in range(5):
        for j in range(5): val=val+M[i][j]*(v[i]*v[j])
    print(f"real root {n}: y enclosure width {width:.1e}; LMI for g={g}: v^T G v in [{float(val.lo):.3e}, {float(val.hi):.3e}]"
          f"  -> certified negative: {val.hi<0}   (whitened min eig {e[g].min():.2f})")
