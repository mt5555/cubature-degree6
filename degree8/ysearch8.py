# Newton search for real solutions of the 10 flatness equations in the 10 degree-9 moments
import sys, numpy as np, pickle
from scipy.optimize import root
from cf8 import hank, yLf, hscale, mineig, geigs, nodes_from_y
from common8 import resid, weights, outside
rng=np.random.default_rng(int(sys.argv[1])); n=int(sys.argv[2])
f=lambda u: hank(yLf*(1+u))/hscale
sols=[]
for t in range(n):
    s=[0.003,0.01,0.03,0.1,0.3,1.0][t%6]
    u0=s*rng.standard_normal(10)
    r=root(f,u0,method='hybr',options={'xtol':1e-14,'maxfev':4000})
    if np.abs(f(r.x)).max()<1e-9: sols.append(r.x)
pickle.dump(sols,open(f'ys8_{sys.argv[1]}.pkl','wb'))
print("seed",sys.argv[1],": converged",len(sols),"of",n)
