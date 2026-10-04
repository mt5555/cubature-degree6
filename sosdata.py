# Exact rational data for the certificate, in coordinates t: y = c + hw*t
import numpy as np, sympy as sp, itertools, pickle
from fractions import Fraction as Fr
from cf2 import K, m6, hank_sym, Wm
lo,hi=np.load('box.npy')
n=8
def monos(D): return [m for d in range(D+1) for m in itertools.combinations_with_replacement(range(n),d)]
def key(m): return tuple(sorted(m))
def fr(x): return Fr(x) if not isinstance(x,(sp.Basic,)) else Fr(int(sp.fraction(x)[0]),int(sp.fraction(x)[1]))
c=[Fr(float(v)) for v in (lo+hi)/2]
hw=[Fr(float(v)) for v in (hi-lo)/2]
W={g:[[Fr(float(Wm[g][i,j])) for j in range(4)] for i in range(4)] for g in Wm}
Kq={g:[[fr(K[g][i,j]) for j in range(4)] for i in range(4)] for g in K}
m6q=[fr(v) for v in m6]
def Graw(g, yv):   # D_g(y) - K_g, yv list of 8 'linear forms' given as dict monomial->Fr
    pass
# G_g(y) entries are affine in y: build coefficient form: entry = const + sum_k coef_k * y_k
def Gaff(g):
    const=[[Fr(0)]*4 for _ in range(4)]; lin=[[[Fr(0)]*8 for _ in range(4)] for _ in range(4)]
    for i in range(4):
        for j in range(4):
            k=i+j
            if g=='x': lin[i][j][k]+=1
            elif g=='y': lin[i][j][k+1]+=1
            else: const[i][j]+=m6q[k]; lin[i][j][k]-=1; lin[i][j][k+1]-=1
            const[i][j]-=Kq[g][i][j]
    return const,lin
def matmul(A,B): return [[sum(A[i][l]*B[l][j] for l in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def T(A): return [list(r) for r in zip(*A)]
GT={}   # g -> (B0, [B_k]) exact, Gtilde(t) = B0 + sum t_k B_k
for g in ['x','y','h']:
    const,lin=Gaff(g)
    # substitute y_k = c_k + hw_k t_k
    C0=[[const[i][j]+sum(lin[i][j][k]*c[k] for k in range(8)) for j in range(4)] for i in range(4)]
    Ck=[[[lin[i][j][k]*hw[k] for j in range(4)] for i in range(4)] for k in range(8)]
    Wg=W[g]; B0=matmul(matmul(Wg,C0),T(Wg)); Bk=[matmul(matmul(Wg,Ck[k]),T(Wg)) for k in range(8)]
    GT[g]=(B0,Bk)
ts=sp.symbols('t0:8')
EQ=[]
for h in hank_sym([sp.Rational(c[k].numerator,c[k].denominator)+sp.Rational(hw[k].numerator,hw[k].denominator)*ts[k] for k in range(8)]):
    P=sp.Poly(sp.expand(h),*ts,domain='QQ'); d={}
    for mon,co in P.terms():
        m=tuple(i for i in range(n) for _ in range(mon[i])); d[m]=fr(co)
    mx=max(abs(float(v)) for v in d.values()); s=Fr(1/mx)
    EQ.append({m:v*s for m,v in d.items()})
pickle.dump({'GT':GT,'EQ':EQ,'c':c,'hw':hw,'W':W},open('sosdata.pkl','wb'))
print("B0 min eig per g:",{g:np.linalg.eigvalsh(np.array(GT[g][0],dtype=float)).min() for g in GT})
print("EQ terms:",[len(e) for e in EQ])
