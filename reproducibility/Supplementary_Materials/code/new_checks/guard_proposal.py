"""One finite Phase-I grid under the revised quadratic guard.
The LP engine is inherited unchanged; its output is an exact proposal, not truth.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import sys, json, hashlib, time
ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent/'inherited/fs2/EIS_PEM_Math_FS2_Sequential_Cost'
sys.path.insert(0,str(BASE/'vendor/software'))
from exact_lp import solve_lp

def q(x): return F(str(x))
def pack(x):
    if isinstance(x,F): return str(x)
    if isinstance(x,dict): return {k:pack(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [pack(v) for v in x]
    return x

def ceilf(x): return -(-x.numerator//x.denominator)
def csqrt(x):
    n=isqrt(x.numerator//x.denominator)
    return n if n*n*x.denominator==x.numerator else n+1

def main(inp_path):
    raw=Path(inp_path).read_bytes(); inp=json.loads(raw)
    ws=list(map(q,inp['model']['parameters']));y=list(map(q,inp['y']));ds=list(map(q,inp['delta']))
    assert inp['model']['kind']=='rc' and ws[0]==1 and len(set(ds))==1
    b=ds[0]/3; M=F(13,4)*(y[0]+ds[0]);assert M>0 and b>0
    n=max(1,csqrt(4*(M+1)/b)); oldn=ceilf(8*(M+1)/b)
    xs=[F(-1)+F(2*i,n) for i in range(n+1)]
    A=[]
    for w in ws:
        A.append([1/(1+w*w*(1+x/2)**2) for x in xs]+[F(1),F(0)])
        A.append([w*(1+x/2)/(1+w*w*(1+x/2)**2) for x in xs]+[F(0),-w])
    nv=len(xs)+2
    AA=[r+[F(1)] for r in A]+[[-z for z in r]+[F(1)] for r in A]+[[F(0)]*nv+[F(1)]]
    bb=[z+d for z,d in zip(y,ds)]+[d-z for z,d in zip(y,ds)]+[min(ds)]
    cc=[F(0)]*nv+[F(1)]
    # First seek a feasible point in a strictly tightened finite box.  The
    # binary64 vector is rounded once and accepted only after exact checking
    # of the ORIGINAL box.  No LP optimality is needed for a strict anchor.
    import numpy as np
    from scipy.optimize import linprog
    t0=time.perf_counter()
    tight=ds[0]-b/2
    ar=A+[[-z for z in row] for row in A]
    br=[z+tight for z in y]+[tight-z for z in y]
    rr=linprog(np.zeros(nv),A_ub=np.array(ar,dtype=float),
               b_ub=np.array(br,dtype=float),bounds=(0,None),method='highs-ds',
               options={'primal_feasibility_tolerance':1e-9,
                        'dual_feasibility_tolerance':1e-9})
    if rr.success:
        scale=1<<70
        x=[F(max(0,round(float(z)*scale)),scale) for z in rr.x]
        margin=min(d-abs(sum(aa*zz for aa,zz in zip(row,x))-yy)
                   for row,yy,d in zip(A,y,ds))
    else:
        margin=F(-1)
    if margin>0:
        sol=x+[margin]
        stats={'method':'tightened_feasibility_proposal_exact_primal_check',
               'float_status':int(rr.status),'optimality_claimed':False}
    else:
        res,stats=solve_lp(AA,bb,cc)
        assert res is not None,'no finite-grid optimum'
        sol,dual=res
    elapsed=time.perf_counter()-t0
    assert sol[-1]>0,'nonpositive margin'
    assert min(sol)>=0
    assert all(sum(a*z for a,z in zip(r,sol))<=v for r,v in zip(AA,bb))
    atoms=[(x,w) for x,w in zip(xs,sol[:len(xs)]) if w]
    theta=sol[len(xs):len(xs)+2]
    out={'format':'quadratic_guard_v1','input_sha256':hashlib.sha256(raw).hexdigest(),
         'input_name':Path(inp_path).name,'b':b,'physical_mass_bound':M,
         'cells':n,'nodes':n+1,'old_lipschitz_nodes':oldn+1,'max_cell':F(2,n),
         'atoms':atoms,'theta':theta,'margin':sol[-1],
         'solver_stats':stats,'proposal_seconds':elapsed}
    path=ROOT/'certificates'/('guard_'+Path(inp_path).name)
    path.write_text(json.dumps(pack(out),indent=2))
    print(json.dumps({'input':Path(inp_path).name,'nodes':n+1,'old_nodes':oldn+1,
                     'margin_over_b':float(sol[-1]/b),'seconds':elapsed,'stats':stats}),flush=True)

if __name__=='__main__': main(sys.argv[1])
