# sample exact 10-pt rules (free), then locally maximize min whitened eig over the solution set (node space)
import numpy as np, sys
from scipy.optimize import least_squares, minimize
from common import resid
from cf2 import *
N=10; rng=np.random.default_rng(int(sys.argv[1])); vals=[]; asc=[]
def y7of(x,y,w): return np.array([np.sum(w*x**(7-k)*y**k) for k in range(8)])
def eigs_all(u):
    x,y,w=u[:N],u[N:2*N],u[2*N:3*N]; e=geigs(y7of(x,y,w)); return np.concatenate([e['x'],e['y'],e['h']])
f=lambda z: resid(z[:N],z[N:2*N],z[2*N:]**2)
for t in range(int(sys.argv[2])):
    lam=rng.dirichlet([1,1,1],N)
    z0=np.concatenate([lam[:,0],lam[:,1],np.sqrt(np.full(N,0.05))*(0.5+rng.random(N))])
    r=least_squares(f,z0,method='trf',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=4000)
    if np.linalg.norm(r.fun)>1e-10: continue
    u=np.concatenate([r.x[:2*N],r.x[2*N:]**2]); e0=eigs_all(u).min(); vals.append(e0)
    if len(asc) < int(sys.argv[3]):
        v0=np.concatenate([u,[e0]])
        cons=[{'type':'eq','fun':lambda v: resid(v[:N],v[N:2*N],v[2*N:3*N])},
              {'type':'ineq','fun':lambda v: eigs_all(v[:3*N])-v[-1]}]
        rr=minimize(lambda v:-v[-1],v0,method='SLSQP',constraints=cons,options={'ftol':1e-14,'maxiter':1500})
        v=rr.x; ok=np.linalg.norm(resid(v[:N],v[N:2*N],v[2*N:3*N]))<1e-10
        if ok: asc.append((eigs_all(v[:3*N]).min(),e0)); np.save('asc_%s_%d.npy'%(sys.argv[1],len(asc)),v)
vals=np.array(vals)
print("seed",sys.argv[1],"exact rules:",len(vals)," max min-eig over samples %.4f"%vals.max(),
      " ascents:",sorted(["%.4f"%a for a,_ in asc],reverse=True)[:8])
