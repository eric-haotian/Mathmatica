"""Exact rational primitives and fixed raw RC observation model.
Only Python's standard library is required by the verifier.
Polynomials use ascending coefficients. No floating optimizer is trusted.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial, isqrt


def frac(x):
    if isinstance(x, F): return x
    if isinstance(x, (int, str)): return F(x)
    raise TypeError('rational values must be integer or fraction strings')

def trim(p):
    p=list(map(frac,p)) or [F(0)]
    while len(p)>1 and not p[-1]: p.pop()
    return p

def add(a,b):
    p=[F(0)]*max(len(a),len(b))
    for i,v in enumerate(a): p[i]+=v
    for i,v in enumerate(b): p[i]+=v
    return trim(p)

def scale(p,c): return trim([v*frac(c) for v in p])
def sub(a,b): return add(a,scale(b,-1))
def mul(a,b):
    p=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): p[i+j]+=x*y
    return trim(p)
def power(p,n):
    ans=[F(1)]
    for _ in range(n): ans=mul(ans,p)
    return ans

def deriv(p): return trim([i*p[i] for i in range(1,len(p))])
def val(p,x):
    x=frac(x); y=F(0)
    for v in reversed(p): y=y*x+v
    return y

def compose_affine(p,a,b):
    ans=[F(0)]
    for v in reversed(p): ans=add(mul(ans,[a,b]),[v])
    return ans

def divide(a,b):
    a=trim(a); b=trim(b)
    if b==[0]: raise ZeroDivisionError()
    q=[F(0)]*max(1,len(a)-len(b)+1)
    while len(a)>=len(b) and a!=[0]:
        n=len(a)-len(b); c=a[-1]/b[-1]; q[n]+=c
        a=sub(a,[F(0)]*n+scale(b,c))
    return trim(q),trim(a)

def solve(A,b):
    n=len(A)
    if len(b)!=n or any(len(r)!=n for r in A): raise ValueError('square system required')
    T=[[frac(x) for x in row]+[frac(v)] for row,v in zip(A,b)]
    for j in range(n):
        pivot=next((i for i in range(j,n) if T[i][j]),None)
        if pivot is None: raise ValueError('singular exact system')
        T[j],T[pivot]=T[pivot],T[j]
        z=T[j][j]; T[j]=[x/z for x in T[j]]
        for i in range(n):
            if i!=j and T[i][j]:
                c=T[i][j]; T[i]=[x-c*y for x,y in zip(T[i],T[j])]
    return [r[-1] for r in T]

def bernstein(p,a,b):
    p=compose_affine(p,a,b-a); n=len(p)-1
    return [sum((p[j]*F(comb(i,j),comb(n,j)) for j in range(i+1)),F(0)) for i in range(n+1)]

@lru_cache(None)
def model(ws=(1,2,3,4)):
    ws=tuple(ws); m=len(ws); D=[F(1)]
    for w in ws: D=mul(D,[F(1),F(0),F(w*w)])
    D0=val(D,1); N=[]
    for w in ws:
        q,r=divide(D,[F(1),F(0),F(w*w)]); assert r==[0]
        N.extend([q,mul([F(0),F(w)],q)])
    C=[[p[i] if i<len(p) else F(0) for p in N] for i in range(2*m)]
    lam=[]
    for j in range(2*m-2):
        target=scale(mul([F(0),F(0),F(1)],power([F(-2),F(2)],j)),D0)
        v=solve(C,target+[F(0)]*(2*m-len(target)))
        assert sum(v[::2])==0 and sum((v[2*i+1]*ws[i] for i in range(m)),F(0))==0
        out=[F(0)]
        for c,p in zip(v,N): out=add(out,scale(p,c))
        assert out==trim(target)
        lam.append(v)
    return {'ws':ws,'D':D,'D0':D0,'numerators':N,'lambda':lam}

def lam(p,ws=(1,2,3,4)):
    mo=model(tuple(ws)); p=trim(p)
    if len(p)>len(mo['lambda']): raise ValueError('polynomial exceeds calibrated observation span')
    return [sum((p[j]*mo['lambda'][j][i] for j in range(len(p))),F(0)) for i in range(2*len(ws))]

def l1(v): return sum(map(abs,v),F(0))
def dot(a,b): return sum((x*y for x,y in zip(a,b)),F(0))
def rho(x,ws=(1,2,3,4)):
    mo=model(tuple(ws)); t=1+frac(x)/2
    return val(mo['D'],t)/(mo['D0']*t*t)

def weight(x,ws=(1,2,3,4)): return 1/rho(x,ws)
def phi(x,ws=(1,2,3,4)):
    t=1+frac(x)/2; out=[]
    for w in ws: out.extend([1/(1+w*w*t*t),w*t/(1+w*w*t*t)])
    return out

def raw(atoms,theta,ws=(1,2,3,4)):
    out=[F(0)]*(2*len(ws))
    for x,mass in atoms:
        if mass<0: raise ValueError('negative spectrum mass')
        out=[a+mass*b for a,b in zip(out,phi(x,ws))]
    if min(theta)<0: raise ValueError('negative calibration')
    for i,w in enumerate(ws): out[2*i]+=theta[0]; out[2*i+1]-=w*theta[1]
    return out

def physical(nu,ws=(1,2,3,4)):
    out={}
    for x,p in nu:
        out[x]=out.get(x,F(0))+p*rho(x,ws)
    return sorted((x,p) for x,p in out.items() if p)

def ref_nu(a):
    a=frac(a)
    return [(F(0),F(1))] if not a else [(-a,F(1,2)),(a,F(1,2))]

def Wpoly(u): return [u*u,F(0),-2*u,F(0),F(1)]

def hermite_resolvent(u,b):
    W=Wpoly(u); den=(b*b-u)**2
    H=[(b**3-2*u*b)/den,(b*b-2*u)/den,b/den,1/den]
    assert mul([b,-F(1)],H)==sub([F(1)],scale(W,1/den))
    return H,den

def laurent_deriv_bound(poly,shift=0,order=0):
    # sum_j poly[j] t^(j-shift), d/dx=(1/2)d/dt; t in [1/2,3/2]
    ans=F(0)
    for j,c in enumerate(poly):
        exponent=j-shift; mult=1
        for r in range(order): mult*=exponent-r
        e=exponent-order
        ans+=abs(c*mult)*((F(3,2)**e) if e>=0 else (F(1,2)**e))/2**order
    return ans

@lru_cache(None)
def norms(ws=(1,2,3,4)):
    mo=model(tuple(ws)); D=mo['D']; Dp=deriv(D); D0=mo['D0']
    rho_bounds=[laurent_deriv_bound(scale(D,1/D0),2,k) for k in range(9)]
    numerator=[F(0),F(0),D0]; wb=[]
    # Exact Bernstein rational derivative bounds, cells cover t in [1/2,3/2].
    for j in range(5):
        bound=F(0)
        for cell in range(32):
            a=F(1,2)+F(cell,32); b=a+F(1,32)
            nb=max(map(abs,bernstein(numerator,a,b)))
            # D has nonnegative coefficients and t>0, hence monotone increasing.
            bound=max(bound,nb/(val(D,a)**(j+1)))
        wb.append(bound)
        numerator=scale(sub(mul(deriv(numerator),D),scale(mul(numerator,Dp),j+1)),F(1,2))
    S2=max(wb[0]/16,wb[1]/16+3*wb[0],wb[2]/16+6*wb[1]+18*wb[0])
    S2=F((S2.numerator+S2.denominator-1)//S2.denominator)
    # f=w/(2-x); unscaled C4 norm bound.
    B4=max(sum((F(comb(j,k))*wb[k]*factorial(j-k) for k in range(j+1)),F(0)) for j in range(5))
    S4=F((B4.numerator+B4.denominator-1)//B4.denominator)
    C3=max(4*wb[0],4*wb[1]+12*wb[0],4*wb[2]+24*wb[1]+24*wb[0],
           4*wb[3]+36*wb[2]+72*wb[1]+24*wb[0])
    S3=F((C3.numerator+C3.denominator-1)//C3.denominator)
    M4=max(laurent_deriv_bound(scale(p,1/D0),2,4) for p in mo['numerators'])
    return {'rho':rho_bounds,'weight':wb,'S2':S2,'S3':S3,'S4':S4,'raw_u_derivative4':M4}

def jsonable(x):
    if isinstance(x,F): return str(x)
    if isinstance(x,dict): return {k:jsonable(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [jsonable(v) for v in x]
    return x


def bumpG(x,a):
    if not a: return F(0)
    if abs(x)>=2*a: return F(0)
    return a*a*(1-x*x/(4*a*a))**3


def original_target_bump(x,a): return weight(x)*bumpG(x,a)/norms()['S2']


def floor_sqrt_dyadic(q,bits=128):
    q=frac(q)
    if q<0: raise ValueError('negative square root')
    n=isqrt((q.numerator << (2*bits))//q.denominator)
    r=F(n,2**bits)
    assert r*r<=q<(r+F(1,2**bits))**2
    return r

def C3_shape(x,b):
    z=abs(frac(x)); b=frac(b)
    if z<=b: return z**4/b
    h=z-b
    return b**3+4*b*b*h+6*b*h*h+4*h**3

def original_target_C3(x,b):return weight(x)*C3_shape(x,b)/norms()['S3']
