import sys, numpy as np
from scipy.optimize import least_squares
from common import *
N = int(sys.argv[1]); mode = sys.argv[2]; trials = int(sys.argv[3]); seed=int(sys.argv[4])
rng = np.random.default_rng(seed)
def unpack(z):
    if mode == 'inside':
        s = z[:3*N].reshape(N,3)**2 + 1e-300
        lam = s/s.sum(1,keepdims=True)
        x, y = lam[:,0], lam[:,1]
    else:
        x, y = z[:N], z[N:2*N]
    w = z[-N:]**2
    return x, y, w
f = lambda z: resid(*unpack(z))
res=[]
for t in range(trials):
    if mode=='inside':
        z0 = np.concatenate([rng.random(3*N), np.sqrt(np.full(N,0.5/N))*(0.5+rng.random(N))])
    else:
        lam = rng.dirichlet([1,1,1],N)
        z0 = np.concatenate([lam[:,0],lam[:,1], np.sqrt(np.full(N,0.5/N))*(0.5+rng.random(N))])
    r = least_squares(f, z0, method='trf', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
    x,y,w = unpack(r.x)
    res.append((np.linalg.norm(r.fun), outside(x,y).max(), w.min(), r.x))
res.sort(key=lambda t:t[0])
print(f"N={N} mode={mode} trials={trials}")
for k in range(8):
    print("  resid %.3e  max_outside %.3e  wmin %.3e"%res[k][:3])
np.save(f"res_{N}_{mode}_{seed}.npy", np.array([np.concatenate([[a,b,c],z]) for a,b,c,z in res]))
