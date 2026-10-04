# node-space search for 15-point degree-8 rules: free (outside allowed) or inside-constrained
import sys, numpy as np
from scipy.optimize import least_squares
from common8 import resid, outside
N=15; mode=sys.argv[1]; trials=int(sys.argv[2]); rng=np.random.default_rng(int(sys.argv[3])); out=[]
def unpack(z):
    if mode=='inside':
        s=z[:3*N].reshape(N,3)**2+1e-300; lam=s/s.sum(1,keepdims=True); x,y=lam[:,0],lam[:,1]
    else: x,y=z[:N],z[N:2*N]
    return x,y,z[-N:]**2
f=lambda z: resid(*unpack(z))
for t in range(trials):
    lam=rng.dirichlet([1,1,1],N)
    if mode=='inside': z0=np.concatenate([np.sqrt(lam).ravel(),np.sqrt(np.full(N,0.5/N))*(0.5+rng.random(N))])
    else: z0=np.concatenate([lam[:,0],lam[:,1],np.sqrt(np.full(N,0.5/N))*(0.5+rng.random(N))])
    r=least_squares(f,z0,method='trf',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=6000)
    x,y,w=unpack(r.x); out.append(np.concatenate([[np.linalg.norm(r.fun),outside(x,y).max()],x,y,w]))
out=np.array(out); out=out[np.argsort(out[:,0])]; np.save(f'ns8_{mode}_{sys.argv[3]}.npy',out)
print(f"{mode}: {trials} starts; exact (<1e-10): {(out[:,0]<1e-10).sum()}; best resid {out[0,0]:.3e}")
print("  best 6:",["%.4f(out %.3f)"%(a,b) for a,b in out[:6,:2]])
