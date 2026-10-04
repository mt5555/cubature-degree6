# Compact exact formulation. Unknowns y_k = L(x^(7-k) y^k), k=0..7.
# Localizing:  M_g >= 0  <=>  G_g(y) := D_g(y) - K_g >= 0   (4x4, Schur complement on the Pi_2 block)
#   D_x = Hank4(y_0..y_6), D_y = Hank4(y_1..y_7), D_h = Hank4(m6) - D_x - D_y
# Flatness: C = C0 + R^T S^{-1} R,  R = Hank45(y) - R0, must be Hankel (6 eqs).
import sympy as sp, numpy as np
from math import factorial
Q = sp.Rational
def mom(a, b): return Q(factorial(a)*factorial(b), factorial(a+b+2))
def mono(d): return [(a, k-a) for k in range(d+1) for a in range(k, -1, -1)]
P2 = mono(2); P3c = [(3-i, i) for i in range(4)]; P4c = [(4-j, j) for j in range(5)]
m6 = [mom(6-k, k) for k in range(7)]
yL = [mom(7-k, k) for k in range(8)]          # Lebesgue degree-7 moments
gs = {'x': [((1,0),1)], 'y': [((0,1),1)], 'h': [((0,0),1),((1,0),-1),((0,1),-1)]}
def Lg(g, a, b): return sum(c*mom(a+da, b+db) for (da,db),c in gs[g])
K = {}
for g in gs:
    A = sp.Matrix(6, 6, lambda i,j: Lg(g, P2[i][0]+P2[j][0], P2[i][1]+P2[j][1]))
    F = sp.Matrix(6, 4, lambda i,j: Lg(g, P2[i][0]+P3c[j][0], P2[i][1]+P3c[j][1]))
    K[g] = F.T * A.inv() * F
# flatness pieces
M2 = sp.Matrix(6, 6, lambda i,j: mom(P2[i][0]+P2[j][0], P2[i][1]+P2[j][1]))
F3 = sp.Matrix(6, 4, lambda i,j: mom(P2[i][0]+P3c[j][0], P2[i][1]+P3c[j][1]))
D3 = sp.Matrix(4, 4, lambda i,j: mom(P3c[i][0]+P3c[j][0], P3c[i][1]+P3c[j][1]))
Ba = sp.Matrix(6, 5, lambda i,j: mom(P2[i][0]+P4c[j][0], P2[i][1]+P4c[j][1]))
M2i = M2.inv()
S = D3 - F3.T*M2i*F3
Sinv = S.inv()
R0 = F3.T*M2i*Ba
C0 = Ba.T*M2i*Ba
pairs = []
for s in range(9):
    idx = [(j, s-j) for j in range(5) if 0 <= s-j <= 4 and j <= s-j]
    for t in idx[1:]: pairs.append((idx[0], t))
assert len(pairs) == 6
ys = sp.symbols('y0:8')
def hank4(v, off): return sp.Matrix(4, 4, lambda i,j: v[i+j+off])
def G_sym(y):
    Dx = hank4(y, 0); Dy = hank4(y, 1); Dh = hank4(m6, 0) - Dx - Dy
    return {'x': Dx - K['x'], 'y': Dy - K['y'], 'h': Dh - K['h']}
def hank_sym(y):
    R = sp.Matrix(4, 5, lambda i,j: y[i+j]) - R0
    C = C0 + R.T*Sinv*R
    return [sp.expand(C[p] - C[q]) for p,q in pairs]
# float versions
Kf = {g: np.array(K[g], dtype=float) for g in K}
m6f = np.array(m6, dtype=float); yLf = np.array(yL, dtype=float)
Sinvf = np.array(Sinv, dtype=float); R0f = np.array(R0, dtype=float); C0f = np.array(C0, dtype=float)
def H(v, off, n=4, m=4): return np.array([[v[i+j+off] for j in range(m)] for i in range(n)])
def G(y):
    Dx = H(y,0); Dy = H(y,1); Dh = H(m6f,0) - Dx - Dy
    return {'x': Dx-Kf['x'], 'y': Dy-Kf['y'], 'h': Dh-Kf['h']}
def hank(y):
    R = H(y,0,4,5) - R0f; C = C0f + R.T@Sinvf@R
    return np.array([C[p]-C[q] for p,q in pairs])
GL = G(yLf)
Wm = {g: np.linalg.inv(np.linalg.cholesky(GL[g])) for g in GL}   # whitening: Lebesgue -> I
def geigs(y):
    Gy = G(y); return {g: np.linalg.eigvalsh(Wm[g]@Gy[g]@Wm[g].T) for g in Gy}
def mineig(y): return min(v.min() for v in geigs(y).values())
