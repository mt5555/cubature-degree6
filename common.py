import numpy as np
from math import factorial
# reference triangle T = {x>=0, y>=0, x+y<=1}, area 1/2
def mono(deg):
    return [(a, d-a) for d in range(deg+1) for a in range(d, -1, -1)]
def mom(a, b):
    return factorial(a)*factorial(b)/factorial(a+b+2)
M6 = mono(6)                                   # 28 monomials
m = np.array([mom(a, b) for a, b in M6])
G = np.array([[mom(a1+a2, b1+b2) for a2, b2 in M6] for a1, b1 in M6])
L = np.linalg.cholesky(G)
Linv = np.linalg.inv(L)
def V(x, y):
    return np.array([x**a * y**b for a, b in M6])   # 28 x n
def resid(x, y, w):
    # error functional in orthonormal-polynomial coordinates
    return Linv @ (V(x, y) @ w - m)
def outside(x, y):
    # signed distance-ish: max violation of x>=0, y>=0, x+y<=1 (positive = outside)
    return np.max(np.stack([-x, -y, (x+y-1)/np.sqrt(2)]), axis=0)
