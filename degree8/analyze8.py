import numpy as np, pickle, glob
from cf8 import hank, yLf, hscale, geigs, nodes_from_y, Wm
from common8 import resid, weights, outside
from scipy.optimize import approx_fprime
sols=[s for f in glob.glob('ys8_*.pkl') for s in pickle.load(open(f,'rb'))]
uniq=[]
for u in sols:
    if all(np.abs(u-v).max()>1e-6 for v in uniq): uniq.append(u)
print("distinct real solutions found:",len(uniq))
rows=[]
for u in uniq:
    y=yLf*(1+u); x,yy=nodes_from_y(y)
    imag=max(np.abs(x.imag).max(),np.abs(yy.imag).max()); x,yy=x.real,yy.real
    w=weights(x,yy); r=np.linalg.norm(resid(x,yy,w))
    J=np.array([approx_fprime(u,lambda v,i=i: hank(yLf*(1+v))[i]/hscale,1e-8) for i in range(10)])
    sv=np.linalg.svd(J,compute_uv=False)
    e=geigs(y); me=min(v.min() for v in e.values())
    viol=outside(x,yy); rows.append((viol.max(),(viol>1e-12).sum(),me,w.min(),r,imag,sv[-1]/sv[0],x,yy,w,y))
rows.sort(key=lambda t:t[0])
print(" maxviol  #out  minWhitEig   wmin      noderesid  imag    Jcond^-1")
for t in rows: print(" %7.4f  %3d   %8.4f   %9.2e  %8.1e  %6.0e  %8.1e"%t[:7])
pickle.dump(rows,open('analyzed8.pkl','wb'))
