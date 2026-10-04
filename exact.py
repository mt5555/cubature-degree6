from fractions import Fraction as Fr
def ldl_pd(M):
    """exact check that symmetric rational matrix M is positive definite; returns (ok, min pivot as float)"""
    n=len(M); A=[[Fr(x) for x in r] for r in M]; mp=None
    for k in range(n):
        p=A[k][k]
        if p<=0: return False, float(p)
        mp=p if mp is None or p<mp else mp
        for i in range(k+1,n):
            if A[i][k]!=0:
                f=A[i][k]/p
                for j in range(k+1,n): A[i][j]-=f*A[k][j]
    return True, float(mp)
def solve_exact(G, r):
    n=len(G); A=[[Fr(x) for x in G[i]]+[Fr(r[i])] for i in range(n)]
    for k in range(n):
        piv=next(i for i in range(k,n) if A[i][k]!=0); A[k],A[piv]=A[piv],A[k]
        for i in range(n):
            if i!=k and A[i][k]!=0:
                f=A[i][k]/A[k][k]; A[i]=[a-f*b for a,b in zip(A[i],A[k])]
    return [A[i][n]/A[i][i] for i in range(n)]
