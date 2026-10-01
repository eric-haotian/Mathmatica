"""Generate rational calibration-aware identifiability counterexamples.
No noise realization, optimization or empirical coverage is claimed.
"""
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
from exact import add,mul,scale,power,val,product,solve,dump,digest
ROOT=Path(__file__).resolve().parent

def response(atoms,theta,ws):
    R,L=theta
    out=[]
    for w in ws:
        rr=R; ii=-w*L
        for x,a in atoms:
            t=1+x/2; den=1+(w*t)**2
            rr+=a/den; ii+=a*w*t/den
        out.extend([rr,ii])
    return out

def denominator(ws):
    return product_poly([[F(1),F(0),w*w] for w in ws])

def product_poly(seq):
    p=[F(1)]
    for q in seq:p=mul(p,q)
    return p

def coefficients_in_t(px):
    q=[F(0)]
    for j,c in enumerate(px):q=add(q,scale(power([F(-2),F(2)],j),c))
    return q

def reproduce(W,ws,known):
    # Numerators of Re k and -Im k form a basis of P_{2m-1} in t.
    nums=[]
    for i,w in enumerate(ws):
        d=product_poly([[F(1),F(0),v*v] for j,v in enumerate(ws) if j!=i])
        nums += [d,scale(mul([F(0),F(1)],d),w)]
    d1=product(1+w*w for w in ws)
    rhs=scale(coefficients_in_t(W),d1)
    if not known:rhs=mul([F(0),F(0),F(1)],rhs)
    n=2*len(ws)
    a=[[nums[j][i] if i<len(nums[j]) else F(0) for j in range(n)] for i in range(n)]
    b=[rhs[i] if i<len(rhs) else F(0) for i in range(n)]
    assert len(rhs)<=n
    return solve(a,b)

def build(K,design,known):
    m=K if known else K+1
    if design=='uniform':ws=[F(j) for j in range(1,m+2)]
    elif design=='halfshift':ws=[F(2*j-1,2) for j in range(1,m+2)]
    elif design=='dyadic':ws=[F(2)**j for j in range(m+1)]
    else:raise ValueError(design)
    omega=ws[:-1]; next_w=ws[-1]
    h=F(1,8); xs=[F(j-K)*h for j in range(2*K+1)]
    d1=product(1+w*w for w in omega)
    def rho(x):
        t=1+x/2
        d=product(1+(w*t)**2 for w in omega)
        return d/d1 if known else d/(d1*t*t)
    un=[F(comb(2*K,j),2**(2*K-1))*rho(x) for j,x in enumerate(xs)]
    normalizer=sum(un[1::2])
    mu0=[(x,un[j]/normalizer) for j,x in enumerate(xs) if j%2]
    mu1=[(x,un[j]/normalizer) for j,x in enumerate(xs) if not j%2]
    theta0=[F(1,10),F(1,100)]
    diff=[v-u for u,v in zip(response(mu0,[0,0],omega),response(mu1,[0,0],omega))]
    theta1=list(theta0) if known else [theta0[0]-diff[0],theta0[1]+diff[1]/omega[0]]
    assert min(theta1)>0
    y=response(mu0,theta0,omega)
    assert y==response(mu1,theta1,omega)
    W=product_poly([power([-x,F(1)],2) for x,a in mu0])
    bounds=[]
    for q in range(3):
        bounds.append(sum(abs(W[j])*F(factorial(j),factorial(j-q)) for j in range(q,len(W))))
    S=max(F(1),*bounds)
    f=scale(W,1/S)
    gap=sum(a*val(f,x) for x,a in mu1)
    assert gap>0 and sum(a*val(f,x) for x,a in mu0)==0
    ext=omega+[next_w]
    yy0=response(mu0,theta0,ext); yy1=response(mu1,theta1,ext)
    extra=[yy1[-2]-yy0[-2],yy1[-1]-yy0[-1]]
    assert any(extra)
    lam=reproduce(W,ext,known)
    row={
        'schema':'calibration_threshold_v1','K':K,'mode':'known' if known else 'unknown',
        'frequencies':omega,'additional_frequency':next_w,'reference':mu0,'alternative':mu1,
        'reference_calibration':theta0,'alternative_calibration':theta1,
        'raw_data':y,'error_radius':F(0),'target_coefficients':f,
    }
    cert={'input_sha256':digest(row),'annihilator':W,'normalizer':S,
          'target_gap':gap,'additional_row_difference':extra,
          'extended_design_annihilator_dual':lam}
    ident=f"{row['mode']}_K{K}_{design}"
    dump(ROOT/'inputs'/f'{ident}.json',row);dump(ROOT/'certificates'/f'{ident}.json',cert)
    return ident

if __name__=='__main__':
    ids=[]
    for K in range(1,5):
        for design in ['uniform','halfshift','dyadic']:ids.append(build(K,design,False))
        ids.append(build(K,'uniform',True))
    dump(ROOT/'results/generation.json',{'cases':ids,'count':len(ids),'description':'16 exact noiseless pairs; not 16 optimizations'})
    print(len(ids),'exact pairs generated')
