# Degree-9 node search with variable projection: weights = linear least squares given nodes.
# usage: search9.py N mode(inside|free) starts seed
import sys, numpy as np
from math import factorial
from scipy.optimize import least_squares
D=9; MON=[(a,d-a) for d in range(D+1) for a in range(d,-1,-1)]                 # 55 monomials
mom=lambda a,b: factorial(a)*factorial(b)/factorial(a+b+2)
m=np.array([mom(a,b) for a,b in MON])
G=np.array([[mom(a1+a2,b1+b2) for a2,b2 in MON] for a1,b1 in MON]); Linv=np.linalg.inv(np.linalg.cholesky(G))
b=Linv@m
N=int(sys.argv[1]); mode=sys.argv[2]; starts=int(sys.argv[3]); rng=np.random.default_rng(int(sys.argv[4]))
def nodes(z):
    if mode=='inside':
        s=z.reshape(N,3)**2+1e-300; lam=s/s.sum(1,keepdims=True); return lam[:,0],lam[:,1]
    return z[:N],z[N:]
def A_of(x,y): return Linv@np.array([x**p*y**q for p,q in MON])
def res(z):
    x,y=nodes(z); A=A_of(x,y); w=np.linalg.lstsq(A,b,rcond=None)[0]; return A@w-b
out=[]
for t in range(starts):
    lam=rng.dirichlet([1,1,1],N)
    z0=np.sqrt(lam).ravel() if mode=='inside' else np.concatenate([lam[:,0],lam[:,1]])
    r=least_squares(res,z0,method='trf',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=3000)
    x,y=nodes(r.x); A=A_of(x,y); w=np.linalg.lstsq(A,b,rcond=None)[0]
    viol=np.max(np.stack([-x,-y,(x+y-1)/np.sqrt(2)]),axis=0)
    out.append(np.concatenate([[np.linalg.norm(A@w-b),viol.max(),w.min()],x,y,w]))
out=np.array(out); out=out[np.argsort(out[:,0])]
np.save(f'r9_N{N}_{mode}_{sys.argv[4]}.npy',out)
ex=out[:,0]<1e-10
print(f"N={N} {mode}: {starts} starts, exact rules (<1e-10): {ex.sum()}, exact with all w>0 and inside: {(ex&(out[:,1]<=1e-12)&(out[:,2]>0)).sum()}; best residuals:",
      " ".join("%.2e"%v for v in out[:5,0]))
