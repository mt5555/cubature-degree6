# find an 11-point degree-6 rule maximizing the minimum node distance to the triangle's edges
import numpy as np, sys
from scipy.optimize import minimize
from common import resid
N=11; r0=np.load('res_11_inside_1.npy'); rng=np.random.default_rng(int(sys.argv[1])); best=None
s2=np.sqrt(2)
def dist(u): x,y=u[:N],u[N:2*N]; return np.concatenate([x,y,(1-x-y)/s2])
for t in range(int(sys.argv[2])):
    z=r0[t % len(r0),3:]; s=z[:3*N].reshape(N,3)**2; lam=s/s.sum(1,keepdims=True)
    if t>=len(r0): lam=0.9*lam+0.1*rng.dirichlet([1,1,1],N)
    u0=np.concatenate([lam[:,0],lam[:,1],z[3*N:]**2]); v0=np.concatenate([u0,[dist(u0).min()]])
    cons=[{'type':'eq','fun':lambda v: resid(v[:N],v[N:2*N],v[2*N:3*N])},
          {'type':'ineq','fun':lambda v: dist(v[:3*N])-v[-1]}]
    r=minimize(lambda v:-v[-1],v0,method='SLSQP',constraints=cons,options={'ftol':1e-14,'maxiter':1000})
    v=r.x; e=np.linalg.norm(resid(v[:N],v[N:2*N],v[2*N:3*N]))
    if e<1e-11 and (best is None or v[-1]>best[-1]): best=v
print("seed",sys.argv[1],"best min edge distance %.6f, min weight %.4e"%(best[-1],best[2*N:3*N].min()))
np.save('best11_%s.npy'%sys.argv[1],best)
