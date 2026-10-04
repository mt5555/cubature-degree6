# Degree-8 / 15-point formulation. Unknowns y_k = L(x^(9-k) y^k), k=0..9.
# Localizing:  M4(gL) >= 0  <=>  G_g(y) = D_g(y) - K_g >= 0   (5x5; Schur complement on the fixed P3 block)
#   D_x = Hank5(y_0..y_8), D_y = Hank5(y_1..y_9), D_h = Hank5(m8) - D_x - D_y
# Flatness: C = C0 + R^T S^{-1} R,  R = Hank56(y) - R0  (6x6), must be Hankel: 10 equations.
import sympy as sp, numpy as np
from math import factorial
Q = sp.Rational
def mom(a, b): return Q(factorial(a)*factorial(b), factorial(a+b+2))
def mono(d): return [(a, k-a) for k in range(d+1) for a in range(k, -1, -1)]
P3 = mono(3); D4 = [(4-i, i) for i in range(5)]; D5 = [(5-j, j) for j in range(6)]
m8 = [mom(8-k, k) for k in range(9)]
yL = [mom(9-k, k) for k in range(10)]
gs = {'x': [((1,0),1)], 'y': [((0,1),1)], 'h': [((0,0),1),((1,0),-1),((0,1),-1)]}
def Lg(g, a, b): return sum(c*mom(a+da, b+db) for (da,db),c in gs[g])
K = {}
for g in gs:
    A = sp.Matrix(10, 10, lambda i,j: Lg(g, P3[i][0]+P3[j][0], P3[i][1]+P3[j][1]))
    F = sp.Matrix(10, 5, lambda i,j: Lg(g, P3[i][0]+D4[j][0], P3[i][1]+D4[j][1]))
    K[g] = F.T * A.LUsolve(F)
M3 = sp.Matrix(10, 10, lambda i,j: mom(P3[i][0]+P3[j][0], P3[i][1]+P3[j][1]))
F4 = sp.Matrix(10, 5, lambda i,j: mom(P3[i][0]+D4[j][0], P3[i][1]+D4[j][1]))
DD = sp.Matrix(5, 5, lambda i,j: mom(D4[i][0]+D4[j][0], D4[i][1]+D4[j][1]))
Ba = sp.Matrix(10, 6, lambda i,j: mom(P3[i][0]+D5[j][0], P3[i][1]+D5[j][1]))
M3iF = M3.LUsolve(F4); M3iB = M3.LUsolve(Ba)
S = DD - F4.T*M3iF
Sinv = S.inv()
R0 = F4.T*M3iB
C0 = Ba.T*M3iB
pairs = []
for s in range(11):
    idx = [(j, s-j) for j in range(6) if 0 <= s-j <= 5 and j <= s-j]
    for tt in idx[1:]: pairs.append((idx[0], tt))
assert len(pairs) == 10
ys = sp.symbols('y0:10')
def hank_sym(y):
    R = sp.Matrix(5, 6, lambda i,j: y[i+j]) - R0
    C = C0 + R.T*Sinv*R
    return [sp.expand(C[p] - C[q]) for p,q in pairs]
def G_sym(y):
    Dx = sp.Matrix(5,5,lambda i,j: y[i+j]); Dy = sp.Matrix(5,5,lambda i,j: y[i+j+1])
    Dh = sp.Matrix(5,5,lambda i,j: m8[i+j]) - Dx - Dy
    return {'x': Dx-K['x'], 'y': Dy-K['y'], 'h': Dh-K['h']}
# float versions
Kf = {g: np.array(K[g], dtype=float) for g in K}
m8f = np.array(m8, dtype=float); yLf = np.array(yL, dtype=float)
Sinvf = np.array(Sinv, dtype=float); R0f = np.array(R0, dtype=float); C0f = np.array(C0, dtype=float)
def H(v, off, n, m): return np.array([[v[i+j+off] for j in range(m)] for i in range(n)])
def G(y):
    Dx = H(y,0,5,5); Dy = H(y,1,5,5); Dh = H(m8f,0,5,5) - Dx - Dy
    return {'x': Dx-Kf['x'], 'y': Dy-Kf['y'], 'h': Dh-Kf['h']}
def hank(y):
    R = H(y,0,5,6) - R0f; C = C0f + R.T@Sinvf@R
    return np.array([C[p]-C[q] for p,q in pairs])
GL = G(yLf)
Wm = {g: np.linalg.inv(np.linalg.cholesky(GL[g])) for g in GL}
def geigs(y):
    Gy = G(y); return {g: np.linalg.eigvalsh(Wm[g]@Gy[g]@Wm[g].T) for g in Gy}
def mineig(y): return min(v.min() for v in geigs(y).values())
hscale = np.abs(hank(yLf)).max()
# --- recover nodes/weights of the flat extension from y (P4 basis)
P4 = mono(4)
M4f = np.array([[float(mom(a1+a2,b1+b2)) for a2,b2 in P4] for a1,b1 in P4])
def Lfull(a,b,y): return y[b] if a+b==9 else float(mom(a,b))
def nodes_from_y(y, seed=0):
    Mx = np.array([[Lfull(a1+a2+1,b1+b2,y) for a2,b2 in P4] for a1,b1 in P4])
    My = np.array([[Lfull(a1+a2,b1+b2+1,y) for a2,b2 in P4] for a1,b1 in P4])
    X = np.linalg.solve(M4f, Mx); Y = np.linalg.solve(M4f, My)
    c = np.random.default_rng(seed).standard_normal()
    ev, V = np.linalg.eig(X + c*Y)
    xs = np.array([ (np.linalg.lstsq(V[:,[i]], X@V[:,[i]], rcond=None)[0][0,0]) for i in range(15)])
    yy = np.array([ (np.linalg.lstsq(V[:,[i]], Y@V[:,[i]], rcond=None)[0][0,0]) for i in range(15)])
    return xs, yy
