"""Small exact-arithmetic primitives. No floating point in a certificate decision."""
from fractions import Fraction as F
from math import isqrt
import json, hashlib

def frac(x):
    if isinstance(x, bool) or isinstance(x, float):
        raise ValueError('certificate rationals must not be bool or float')
    return F(x)

def need(condition, message):
    if not condition: raise ValueError(message)

def trim(p):
    p=[frac(x) for x in p]
    while len(p)>1 and p[-1]==0: p.pop()
    return p or [F(0)]

def add(a,b):
    p=[F(0)]*max(len(a),len(b))
    for j,x in enumerate(a): p[j]+=x
    for j,x in enumerate(b): p[j]+=x
    return trim(p)

def scale(p,c): return trim([x*c for x in p])

def mul(a,b):
    p=[F(0)]*(len(a)+len(b)-1)
    for j,x in enumerate(a):
        for k,y in enumerate(b):p[j+k]+=x*y
    return trim(p)

def power(p,n):
    need(isinstance(n,int) and n>=0,'power')
    out=[F(1)]
    for _ in range(n):out=mul(out,p)
    return out

def val(p,x):
    out=F(0)
    for a in reversed(p):out=out*x+a
    return out

def product(seq):
    z=F(1)
    for x in seq:z*=x
    return z

def sqrt_bracket(x,bits=200):
    x=frac(x);need(x>=0 and bits>=0,'sqrt arguments')
    d=1<<bits; a=isqrt((x.numerator*d*d)//x.denominator)
    lo=F(a,d);hi=lo if lo*lo==x else F(a+1,d)
    need(lo*lo<=x<=hi*hi,'sqrt bound')
    return lo,hi

def pack(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:pack(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [pack(v) for v in x]
    return x

def digest(x):return hashlib.sha256(json.dumps(pack(x),sort_keys=True,separators=(',',':')).encode()).hexdigest()

def dump(path,x):
    with open(path,'w',encoding='utf8') as f:json.dump(pack(x),f,indent=2,sort_keys=True)

def solve(a,b):
    n=len(b);aa=[[frac(x) for x in row]+[frac(y)] for row,y in zip(a,b)]
    need(len(aa)==n and all(len(row)==n+1 for row in aa),'square solve')
    for j in range(n):
        piv=next((k for k in range(j,n) if aa[k][j]),None)
        need(piv is not None,'singular system')
        aa[j],aa[piv]=aa[piv],aa[j];t=aa[j][j];aa[j]=[x/t for x in aa[j]]
        for k in range(n):
            if k!=j:
                t=aa[k][j]
                if t:aa[k]=[x-t*y for x,y in zip(aa[k],aa[j])]
    return [row[-1] for row in aa]
