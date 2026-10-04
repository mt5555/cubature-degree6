import re
from fractions import Fraction
def load(fn='rur_Q.out'):
    s=open(fn).read().strip().rstrip(':').replace('^','**')
    s=re.sub(r'(-?\d+)\s*/\s*(\d+(?:\*\*\d+)?)',lambda m:f'Fraction({m.group(1)},{m.group(2)})',s)
    return eval(s)
if __name__=="__main__":
    R=load()
    def shape(x):
        if isinstance(x,list):
            if len(x)>3 and all(isinstance(e,int) for e in x): return f"ints[{len(x)}]"
            return "["+", ".join(shape(e) for e in x)+"]"
        if isinstance(x,int): return str(x) if abs(x)<10**6 else "big"
        return type(x).__name__ if not isinstance(x,str) else repr(x)
    print(shape(R)[:2500])
