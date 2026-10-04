import numpy as np
from scipy.optimize import least_squares, minimize
from common import *
import minout as M
N=10
M.rng=np.random.default_rng(11)
n,z=M.solve(0.107,40)
Vt=M.verts(0.107); s=z[:30].reshape(N,3)**2; lam=s/s.sum(1,keepdims=True); P=lam@Vt
best=None
for k in range(30):
    rng=np.random.default_rng(k)
    if k>0:
        M.rng=rng; n,z=M.solve(0.106,3); s=z[:30].reshape(N,3)**2; lam=s/s.sum(1,keepdims=True); P=lam@M.verts(0.106)
        w0=z[30:]**2
    else: w0=z[30:]**2
    u0=np.concatenate([P[:,0],P[:,1],w0,[0.107]])
    cons=[{'type':'eq','fun':lambda u: resid(u[:N],u[N:2*N],u[2*N:3*N])},
          {'type':'ineq','fun':lambda u: u[-1]-np.concatenate([-u[:N],-u[N:2*N],(u[:N]+u[N:2*N]-1)/np.sqrt(2)])}]
    r=minimize(lambda u:u[-1],u0,constraints=cons,method='SLSQP',options={'ftol':1e-14,'maxiter':2000})
    u=r.x; e=np.linalg.norm(resid(u[:N],u[N:2*N],u[2*N:3*N]))
    if e<1e-10 and (best is None or u[-1]<best[-1]): best=u
    print(k, "t=%.10f resid=%.1e"%(u[-1],e), flush=True)
u=best; np.save('extremal.npy',u)
np.set_printoptions(precision=12,suppress=True)
print("min t =",u[-1]); print(np.column_stack([u[:N],u[N:2*N],u[2*N:3*N]]))
