from common import *
import itertools
orbits = [((0.501426509658179,0.249286745170910,0.249286745170910),0.116786275726379),
          ((0.873821971016996,0.063089014491502,0.063089014491502),0.050844906370207),
          ((0.636502499121399,0.310352451033784,0.053145049844817),0.082851075618374)]
P=[];W=[]
for b,w in orbits:
    for p in set(itertools.permutations(b)):
        P.append(p);W.append(w*0.5)
P=np.array(P);W=np.array(W)
print(len(P), np.linalg.norm(resid(P[:,0],P[:,1],W)))
