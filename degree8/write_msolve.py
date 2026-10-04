# write the 10 flatness equations (exact, integer coefficients) in msolve format
import sympy as sp, sys
from cf8 import hank_sym, ys
eqs=hank_sym(list(ys))
lines=[]
for e in eqs:
    num,den=sp.fraction(sp.together(e)); P=sp.Poly(sp.expand(num),*ys)
    c=P.coeffs(); g=sp.gcd_list([sp.Integer(x) for x in c]); P=sp.Poly(P.as_expr()/g,*ys)
    lines.append(str(P.as_expr()).replace('**','^'))
mx=max(len(str(abs(int(c)))) for e in lines for c in sp.Poly(sp.sympify(e.replace('^','**')),*ys).coeffs())
for char,fn in [(0,'flat8_Q.ms'),(1073741827,'flat8_p.ms')]:
    with open(fn,'w') as f:
        f.write(','.join(str(v) for v in ys)+'\n'+str(char)+'\n'+',\n'.join(lines)+'\n')
print("wrote 10 equations; max coefficient digits:",mx, "; terms per eq:",[len(sp.Poly(sp.sympify(l.replace('^','**')),*ys).terms()) for l in lines])
