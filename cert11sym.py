# Certify an EXACTLY mirror-symmetric (x <-> y) 11-point degree-6 rule with all nodes inside T.
# Structure: 3 nodes on the diagonal (t_i, t_i) with weights p_i, and 4 mirror pairs (a_j,b_j),(b_j,a_j)
# with weights q_j: 3*2 + 4*3 = 18 unknowns. A mirror-symmetric rule is exact for degree <= 6 iff it is
# exact for the 16 symmetric combinations x^a y^b + x^b y^a (a >= b), so the system is 16 x 18;
# 2 unknowns are fixed at exact rationals and the remaining 16x16 system is certified with Krawczyk.
import numpy as np, sympy as sp, mpmath as mp, time
from fractions import Fraction as Fr
from math import factorial
from scipy.linalg import qr
t0=time.time(); mp.mp.dps=70
MON=[(a,d-a) for d in range(7) for a in range(d,-1,-1) if a>=d-a]          # 16 orbits, a >= b
mom={ab:Fr(factorial(ab[0])*factorial(ab[1]),factorial(sum(ab)+2)) for ab in MON}
U=sp.symbols('u0:18')
t=U[0:3]; p=U[3:6]; a=U[6:10]; b=U[10:14]; q=U[14:18]
F=[]
for (i,j) in MON:
    if i==j:   # x^i y^i: each mirror pair contributes twice
        e=sum(p[k]*t[k]**(2*i) for k in range(3))+sum(2*q[k]*a[k]**i*b[k]**i for k in range(4))
    else:      # x^i y^j (its mirror x^j y^i gives the same equation)
        e=sum(p[k]*t[k]**(i+j) for k in range(3))+sum(q[k]*(a[k]**i*b[k]**j+b[k]**i*a[k]**j) for k in range(4))
    F.append(sp.expand(e-sp.Rational(mom[(i,j)].numerator,mom[(i,j)].denominator)))
J=sp.Matrix(F).jacobian(U)
# starting point from the max-margin rule (best11_1.npy), sorted into diagonal nodes and mirror pairs
v=np.load('best11_1.npy')[:33]; x,y,w=v[:11],v[11:22],v[22:33]
diag=[i for i in range(11) if abs(x[i]-y[i])<1e-6]; rest=[i for i in range(11) if i not in diag]
pairs=[]
for i in rest:
    if any(i in pr for pr in pairs): continue
    j=min((k for k in rest if k!=i), key=lambda k: abs(x[k]-y[i])+abs(y[k]-x[i])); pairs.append((i,j))
assert len(diag)==3 and len(pairs)==4
u0=[x[i] for i in diag]+[w[i] for i in diag]+[x[i] for i,j in pairs]+[y[i] for i,j in pairs]+[(w[i]+w[j])/2 for i,j in pairs]
Ff=sp.lambdify(U,F,'numpy'); Jf=sp.lambdify(U,J,'numpy')
_,_,piv=qr(np.array(Jf(*u0),dtype=float),pivoting=True)
free=sorted(piv[:16]); fixed=sorted(piv[16:])
names=[f"t{k+1}" for k in range(3)]+[f"p{k+1}" for k in range(3)]+[f"a{k+1}" for k in range(4)]+[f"b{k+1}" for k in range(4)]+[f"q{k+1}" for k in range(4)]
fixval={j:Fr(float(u0[j])) for j in fixed}
print("fixed at exact dyadic rationals:",[names[j] for j in fixed])
Fm=sp.lambdify(U,F,'mpmath'); Jm=sp.lambdify(U,J,'mpmath')
def full(uf):
    u=[None]*18
    for j,val in fixval.items(): u[j]=mp.mpf(val.numerator)/val.denominator
    for k,j in enumerate(free): u[j]=uf[k]
    return u
uf=mp.matrix([mp.mpf(float(u0[j])) for j in free])
for it in range(10):
    u=full(uf); Fv=mp.matrix(Fm(*u)); Jv=mp.matrix(Jm(*u)); Jr=mp.matrix(16,16)
    for r in range(16):
        for k,j in enumerate(free): Jr[r,k]=Jv[r,j]
    uf=uf-mp.lu_solve(Jr,Fv)
print("Newton residual: %s"%mp.nstr(mp.norm(mp.matrix(Fm(*full(uf)))),3))
def tofr(z):
    s,man,ex,_=z._mpf_; f=Fr(int(man))*(Fr(2)**ex if ex>=0 else Fr(1,2**(-ex))); return -f if s else f
cen={j:tofr(uf[k]) for k,j in enumerate(free)}; cen.update(fixval)
class Iv:
    __slots__=('lo','hi')
    def __init__(s,lo,hi=None): s.lo=Fr(lo); s.hi=Fr(lo if hi is None else hi)
    def __add__(s,o): o=o if isinstance(o,Iv) else Iv(o); return Iv(s.lo+o.lo,s.hi+o.hi)
    __radd__=__add__
    def __sub__(s,o): o=o if isinstance(o,Iv) else Iv(o); return Iv(s.lo-o.hi,s.hi-o.lo)
    def __mul__(s,o):
        o=o if isinstance(o,Iv) else Iv(o); pr=[s.lo*o.lo,s.lo*o.hi,s.hi*o.lo,s.hi*o.hi]; return Iv(min(pr),max(pr))
    __rmul__=__mul__
    def mag(s): return max(abs(s.lo),abs(s.hi))
def ieval(expr,X):                       # interval evaluation of a polynomial (term by term)
    P=sp.Poly(expr,*U); acc=Iv(0)
    for mon,co in P.terms():
        term=Iv(Fr(int(sp.fraction(co)[0]),int(sp.fraction(co)[1])))
        for k,e in enumerate(mon):
            for _ in range(e): term=term*X[k]
        acc=acc+term
    return acc
rad=Fr(1,10**30)
X=[Iv(cen[j]-rad,cen[j]+rad) if j in free else Iv(cen[j]) for j in range(18)]
Xc=[Iv(cen[j]) for j in range(18)]
Fc=[ieval(f,Xc).lo for f in F]
JX=[[ieval(J[r,j],X) for j in free] for r in range(16)]
Jc=np.array([[float(ieval(J[r,j],Xc).lo) for j in free] for r in range(16)])
Y=[[Fr(float(e)) for e in row] for row in np.linalg.inv(Jc)]
ok=True; worst=0
for i in range(16):
    YF=sum(Y[i][r]*Fc[r] for r in range(16)); spread=Fr(0)
    for k in range(16):
        Mik=Iv(1 if i==k else 0)-sum((Y[i][r]*JX[r][k] for r in range(16)),Iv(0)); spread+=Mik.mag()*rad
    used=abs(YF)+spread; worst=max(worst,used/rad); ok&= used<rad
print("Krawczyk: K(X) strictly inside X:",ok,"(worst ratio %.2e)"%float(worst))
# build the 11 nodes from the box and check inside / positive over the whole box
T=X[0:3]; Pw=X[3:6]; A=X[6:10]; B=X[10:14]; Qw=X[14:18]
nodes=[(T[k],T[k],Pw[k]) for k in range(3)]+[(A[k],B[k],Qw[k]) for k in range(4)]+[(B[k],A[k],Qw[k]) for k in range(4)]
inside=all(xx.lo>0 and yy.lo>0 and (Iv(1)-xx-yy).lo>0 for xx,yy,_ in nodes); pos=all(ww.lo>0 for _,_,ww in nodes)
dmin=min(min(xx.lo,yy.lo,(Iv(1)-xx-yy).lo/Fr(1414214,1000000)) for xx,yy,_ in nodes)
print("all nodes strictly inside over the box:",inside,"(min edge distance >= %.6f); all weights > 0:"%float(dmin),pos)
print("CERTIFIED EXACTLY MIRROR-SYMMETRIC 11-POINT INSIDE RULE:",ok and inside and pos," [%.0fs]"%(time.time()-t0))
mid=lambda I: (I.lo+I.hi)/2
with open('maxmargin11.txt','w') as f:
    f.write("# 11-point degree-6 rule on T = {x>=0,y>=0,x+y<=1}, all nodes inside, exactly mirror-symmetric (x<->y).\n")
    f.write("# Chosen to (approximately) maximize the minimum node-to-edge distance (>= %.4f).\n"%float(dmin))
    f.write("# Certified by cert11sym.py (Krawczyk, exact rational intervals); values below are the certified\n")
    f.write("# box centres, accurate to ~1e-30, printed to 25 digits.  columns: x  y  weight   (weights sum to 1/2)\n")
    for xx,yy,ww in nodes: f.write("%s %s %s\n"%tuple(mp.nstr(mp.mpf(mid(I).numerator)/mid(I).denominator,25) for I in (xx,yy,ww)))
print(open('maxmargin11.txt').read())
