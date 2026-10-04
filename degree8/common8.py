import numpy as np
from math import factorial
def mono(deg): return [(a, d-a) for d in range(deg+1) for a in range(d, -1, -1)]
def mom(a, b): return factorial(a)*factorial(b)/factorial(a+b+2)
M8 = mono(8)
m = np.array([mom(a, b) for a, b in M8])
G = np.array([[mom(a1+a2, b1+b2) for a2, b2 in M8] for a1, b1 in M8])
L = np.linalg.cholesky(G); Linv = np.linalg.inv(L)
def V(x, y): return np.array([x**a * y**b for a, b in M8])
def resid(x, y, w): return Linv @ (V(x, y) @ w - m)
def weights(x, y):   # least-squares weights for given nodes
    return np.linalg.lstsq(Linv@V(x,y), Linv@m, rcond=None)[0]
def outside(x, y): return np.max(np.stack([-x, -y, (x+y-1)/np.sqrt(2)]), axis=0)
