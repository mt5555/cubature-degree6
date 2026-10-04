# smallest eps such that an exact 10-pt rule exists with nodes in T_eps
# T_eps = {x>=-eps, y>=-eps, x+y<=1+eps*sqrt2}  (offset each edge outward by eps)
import sys, numpy as np
from scipy.optimize import least_squares
from common import *
N=10; rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 and __name__=='__main__' else 0)
def verts(e):
    # vertices of triangle with edges offset by e
    s2=np.sqrt(2); c=1+e*s2+ e  # intersection computations
    return np.array([[-e,-e],[1+e*s2+e,-e],[-e,1+e*s2+e]])
def solve(e, trials):
    Vt=verts(e); best=(1e9,None)
    for t in range(trials):
        z0=np.concatenate([rng.random(3*N), np.sqrt(np.full(N,0.05))*(0.5+rng.random(N))])
        def f(z):
            s=z[:3*N].reshape(N,3)**2+1e-300; lam=s/s.sum(1,keepdims=True)
            P=lam@Vt; return resid(P[:,0],P[:,1],z[3*N:]**2)
        r=least_squares(f,z0,method='trf',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=3000)
        n=np.linalg.norm(r.fun)
        if n<best[0]: best=(n,r.x)
    return best
if __name__=="__main__":
  for e in [float(a) for a in sys.argv[2:]]:
    n,_=solve(e,60); print("eps %.4f  best resid %.3e"%(e,n), flush=True)
