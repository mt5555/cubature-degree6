import numpy as np, sys
from cf2 import *
from common import resid, outside
print("Lebesgue G eigs (raw):", {g: np.linalg.eigvalsh(GL[g]).min() for g in GL})
print("hank at Lebesgue (nonzero expected):", np.abs(hank(yLf)).max())
u = np.load('extremal.npy'); N=10; x,y,w = u[:N],u[N:2*N],u[2*N:3*N]
y7 = np.array([np.sum(w*x**(7-k)*y**k) for k in range(8)])
print("extremal rule: hank %.2e, whitened min eigs:"%np.abs(hank(y7)).max(), {g:"%.4f"%v.min() for g,v in geigs(y7).items()})
