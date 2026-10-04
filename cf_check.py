import numpy as np
from scipy.optimize import minimize
from cf import *
u=np.load('extremal.npy'); N=10; x,y,w=u[:N],u[N:2*N],u[2*N:3*N]
y7=np.array([np.sum(w*x**a*y**b) for a,b in D7])
print("pairs:",len(pairs)," hankel resid on extremal rule: %.2e"%np.abs(hank(y7)).max())
print("loc min eigs on extremal rule:",[ "%.3e"%np.linalg.eigvalsh(Mg).min() for Mg in locs(y7)])
y7L=np.array([mom(a,b) for a,b in D7])
print("Lebesgue: hankel resid %.2e, loc min eigs"%np.abs(hank(y7L)).max(), ["%.3e"%np.linalg.eigvalsh(Mg).min() for Mg in locs(y7L)])
# maximize s = min localizing eigenvalue subject to hankel(y7)=0
rng=np.random.default_rng(0); best=-1
sc=1/np.linalg.norm(y7L)
for t in range(200):
    z0=np.concatenate([y7L*(1+0.3*rng.standard_normal(8)) if t>0 else y7,[ -0.01]])
    cons=[{'type':'eq','fun':lambda z: hank(z[:8])*1e4},
          {'type':'ineq','fun':lambda z: np.concatenate([np.linalg.eigvalsh(Mg) for Mg in locs(z[:8])])-z[8]}]
    r=minimize(lambda z:-z[8]*1e4,z0,constraints=cons,method='SLSQP',options={'maxiter':500,'ftol':1e-16})
    z=r.x
    if np.abs(hank(z[:8])).max()<1e-12:
        s=min(np.linalg.eigvalsh(Mg).min() for Mg in locs(z[:8]))
        if best==-1 or s>best: best=s; zb=z
print("max over flat extensions of min localizing eig: %.6e"%best)
np.save('cf_best.npy',zb)
