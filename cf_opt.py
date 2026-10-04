import numpy as np, sys
from scipy.optimize import minimize
from scipy.linalg import sqrtm
from cf import *
y7L=np.array([mom(a,b) for a,b in D7])
S=[np.real(np.linalg.inv(sqrtm(Mg))) for Mg in locs(y7L)]   # whiten so Lebesgue -> identity
def geigs(y7): return np.concatenate([np.linalg.eigvalsh(s@Mg@s) for s,Mg in zip(S,locs(y7))])
hs=1/mom(4,4)
u=np.load('extremal.npy'); N=10; x,y,w=u[:N],u[N:2*N],u[2*N:3*N]
y7e=np.array([np.sum(w*x**a*y**b) for a,b in D7])
print("extremal rule: min whitened loc eig %.4e"%geigs(y7e).min())
rng=np.random.default_rng(int(sys.argv[1])); res=[]
for t in range(int(sys.argv[2])):
    z0=np.concatenate([y7e if t==0 else y7L*(1+0.2*rng.standard_normal(8)),[-0.5]])
    cons=[{'type':'eq','fun':lambda z: hank(z[:8]*1)*hs},
          {'type':'ineq','fun':lambda z: geigs(z[:8])-z[8]}]
    r=minimize(lambda z:-z[8],z0,constraints=cons,method='SLSQP',options={'maxiter':1000,'ftol':1e-15})
    z=r.x; h=np.abs(hank(z[:8])).max()*hs
    if h<1e-10: res.append((geigs(z[:8]).min(),t))
res.sort(reverse=True)
print("converged %d; top values:"%len(res), ["%.6f"%v for v,_ in res[:10]])
