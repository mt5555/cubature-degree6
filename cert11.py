# Certify an 11-point degree-6 rule with all nodes strictly inside T, via the Krawczyk test
# in exact rational interval arithmetic.
import numpy as np, mpmath as mp, pickle, time
from fractions import Fraction as Fr
from math import factorial
from scipy.linalg import qr
t0=time.time()
N=11; mp.mp.dps=70
MON=[(a,d-a) for d in range(7) for a in range(d,-1,-1)]                 # 28 monomials
mom={ab:Fr(factorial(ab[0])*factorial(ab[1]),factorial(ab[0]+ab[1]+2)) for ab in MON}
v=np.load('best11_1.npy')[:3*N]                                         # u = (x_1..x_11, y_1..y_11, w_1..w_11)
def J_float(u):
    x,y,w=u[:N],u[N:2*N],u[2*N:]
    J=np.zeros((28,3*N))
    for r,(a,b) in enumerate(MON):
        J[r,:N]=w*a*x**max(a-1,0)*y**b if a>0 else 0
        J[r,N:2*N]=w*b*x**a*y**max(b-1,0) if b>0 else 0
        J[r,2*N:]=x**a*y**b
    return J
_,_,piv=qr(J_float(v),pivoting=True)
free=sorted(piv[:28]); fixed=sorted(piv[28:])
names=[f"x{i+1}" for i in range(N)]+[f"y{i+1}" for i in range(N)]+[f"w{i+1}" for i in range(N)]
print("fixed (exact rationals):",[names[j] for j in fixed])
fixval={j:Fr(float(v[j])) for j in fixed}                              # exact dyadic values
# --- high-precision Newton on the 28 free unknowns
def full(uf):
    u=[None]*(3*N)
    for j,val in fixval.items(): u[j]=mp.mpf(val.numerator)/val.denominator
    for k,j in enumerate(free): u[j]=uf[k]
    return u
def F_mp(u):
    x,y,w=u[:N],u[N:2*N],u[2*N:]
    return mp.matrix([mp.fsum(w[i]*x[i]**a*y[i]**b for i in range(N))-mp.mpf(mom[(a,b)].numerator)/mom[(a,b)].denominator for a,b in MON])
def J_mp(u):
    x,y,w=u[:N],u[N:2*N],u[2*N:]; J=mp.matrix(28,28)
    for r,(a,b) in enumerate(MON):
        for k,j in enumerate(free):
            i=j%N
            if j<N:     J[r,k]= w[i]*a*x[i]**(a-1)*y[i]**b if a>0 else 0
            elif j<2*N: J[r,k]= w[i]*b*x[i]**a*y[i]**(b-1) if b>0 else 0
            else:       J[r,k]= x[i]**a*y[i]**b
    return J
uf=mp.matrix([mp.mpf(float(v[j])) for j in free])
for it in range(8):
    u=full(uf); Fv=F_mp(u); uf=uf-mp.lu_solve(J_mp(u),Fv)
print("Newton residual: %s  [%.1fs]"%(mp.nstr(mp.norm(F_mp(full(uf))),3),time.time()-t0))
def tofr(x):
    s,man,ex,_=x._mpf_; f=Fr(int(man))*(Fr(2)**ex if ex>=0 else Fr(1,2**(-ex))); return -f if s else f
cen={j:tofr(uf[k]) for k,j in enumerate(free)}; cen.update(fixval)
# --- exact interval arithmetic
class Iv:
    __slots__=('lo','hi')
    def __init__(s,lo,hi=None): s.lo=Fr(lo); s.hi=Fr(lo if hi is None else hi)
    def __add__(s,o): o=o if isinstance(o,Iv) else Iv(o); return Iv(s.lo+o.lo,s.hi+o.hi)
    __radd__=__add__
    def __sub__(s,o): o=o if isinstance(o,Iv) else Iv(o); return Iv(s.lo-o.hi,s.hi-o.lo)
    def __rsub__(s,o): return Iv(o)-s
    def __mul__(s,o):
        o=o if isinstance(o,Iv) else Iv(o); p=[s.lo*o.lo,s.lo*o.hi,s.hi*o.lo,s.hi*o.hi]; return Iv(min(p),max(p))
    __rmul__=__mul__
    def __pow__(s,k):
        r=Iv(1)
        for _ in range(k): r=r*s
        return r
    def mag(s): return max(abs(s.lo),abs(s.hi))
rad=Fr(1,10**30)
X={j:(Iv(cen[j]-rad,cen[j]+rad) if j in free else Iv(cen[j])) for j in range(3*N)}
xc={j:Iv(cen[j]) for j in range(3*N)}
def F_iv(U):
    return [sum((U[2*N+i]*U[i]**a*U[N+i]**b for i in range(N)),Iv(0))-mom[(a,b)] for a,b in MON]
def J_iv(U):
    J=[[None]*28 for _ in range(28)]
    for r,(a,b) in enumerate(MON):
        for k,j in enumerate(free):
            i=j%N
            if j<N:     J[r][k]= U[2*N+i]*a*U[i]**(a-1)*U[N+i]**b if a>0 else Iv(0)
            elif j<2*N: J[r][k]= U[2*N+i]*b*U[i]**a*U[N+i]**(b-1) if b>0 else Iv(0)
            else:       J[r][k]= U[i]**a*U[N+i]**b
    return J
Fc=[f.lo for f in F_iv(xc)]                                             # exact point values
assert all(f.lo==f.hi for f in F_iv(xc))
Jc=np.array([[float(e.lo) for e in row] for row in J_iv(xc)])
Y=[[Fr(float(e)) for e in row] for row in np.linalg.inv(Jc)]           # approximate inverse, exact rationals
JX=J_iv(X)
print("interval Jacobian built  [%.1fs]"%(time.time()-t0))
ok=True; worst=0
for i in range(28):
    YF=sum(Y[i][r]*Fc[r] for r in range(28))
    spread=Fr(0)
    for k in range(28):
        Mik=Iv(1 if i==k else 0)-sum((Y[i][r]*JX[r][k] for r in range(28)),Iv(0))
        spread+=Mik.mag()*rad
    used=abs(YF)+spread                                                 # K_i ⊂ [cen-used, cen+used]
    worst=max(worst,used/rad); ok&= used<rad
print("Krawczyk: K(X) strictly inside X:",ok,"  (worst ratio |K_i - c_i| / r = %.3e)  [%.1fs]"%(float(worst),time.time()-t0))
# --- inside-ness and positivity over the whole box
inside=all(X[i].lo>0 and X[N+i].lo>0 and (Iv(1)-X[i]-X[N+i]).lo>0 for i in range(N))
posw=all(X[2*N+i].lo>0 for i in range(N))
mind=min(min(X[i].lo,X[N+i].lo,(Iv(1)-X[i]-X[N+i]).lo/Fr(1414214,1000000)) for i in range(N))
print("all nodes strictly inside T over the box:",inside," (min edge distance >= %.6f);  all weights > 0:"%float(mind),posw)
print("CERTIFIED 11-POINT INSIDE RULE:",ok and inside and posw)
pickle.dump({'center':cen,'rad':rad,'free':free,'fixed':fixed},open('cert11.pkl','wb'))
print("\n  node    x                     y                     weight (sum = 1/2)")
for i in range(N): print(f"  {i+1:2d}  {float(cen[i]):.17f}  {float(cen[N+i]):.17f}  {float(cen[2*N+i]):.17f}")
