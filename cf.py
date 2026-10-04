# Curto-Fialkow formulation: unknown degree-7 moments y7 (8 numbers)
import numpy as np
from math import factorial
from fractions import Fraction
from common import mono, mom
P3=mono(3); D4=[(a,4-a) for a in range(4,-1,-1)]
D7=[(a,7-a) for a in range(7,-1,-1)]
idx7={k:i for i,k in enumerate(D7)}
def L(a,b,y7):
    return y7[idx7[(a,b)]] if a+b==7 else mom(a,b)
M3=np.array([[mom(a1+a2,b1+b2) for a2,b2 in P3] for a1,b1 in P3])
M3inv=np.linalg.inv(M3)
def C_of(y7):
    B=np.array([[L(a1+a2,b1+b2,y7) for a2,b2 in D4] for a1,b1 in P3])
    return B.T@M3inv@B
# Hankel consistency: entries of C with equal exponent sum must agree
groups={}
for i,(a1,b1) in enumerate(D4):
    for j,(a2,b2) in enumerate(D4):
        if j>=i: groups.setdefault((a1+a2,b1+b2),[]).append((i,j))
pairs=[(g[0],h) for g in groups.values() for h in g[1:]]
def hank(y7):
    C=C_of(y7); return np.array([C[p]-C[q] for p,q in pairs])
def locs(y7):
    out=[]
    for g in [((1,0,1.),),((0,1,1.),),((0,0,1.),(1,0,-1.),(0,1,-1.))]:
        out.append(np.array([[sum(c*L(a1+a2+da,b1+b2+db,y7) for da,db,c in g) for a2,b2 in P3] for a1,b1 in P3]))
    return out
