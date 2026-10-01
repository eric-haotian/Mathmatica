"""Exact algebra and a-posteriori width checks for manuscript v1.2.

Uses only the standard library. These are finite validation checks, not a proof
assistant or new optimization experiments. Endpoint outputs must first pass the
inherited continuous verifiers; the wrapper replay_all.py performs that step.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
from math import comb, factorial
import json, hashlib, copy

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'closure'/'results'

def solve(a: list[list[F]], b: list[F]) -> list[F]:
    n=len(b); a=[r[:]+[b[i]] for i,r in enumerate(a)]
    for j in range(n):
        k=next((i for i in range(j,n) if a[i][j]),None)
        if k is None: raise ValueError('singular exact system')
        a[j],a[k]=a[k],a[j]; v=a[j][j]; a[j]=[x/v for x in a[j]]
        for i in range(n):
            if i!=j:
                v=a[i][j]; a[i]=[x-v*y for x,y in zip(a[i],a[j])]
    return [r[-1] for r in a]

def cadd(a,b): return (a[0]+b[0],a[1]+b[1])
def cmul(a,b): return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def cdiv(a,b):
    d=b[0]*b[0]+b[1]*b[1]
    if not d: raise ZeroDivisionError('complex denominator')
    return ((a[0]*b[0]+a[1]*b[1])/d,(a[1]*b[0]-a[0]*b[1])/d)
def cpow(a,n):
    z=(F(1),F(0))
    for _ in range(n):z=cmul(z,a)
    return z

def width_bracket(plus: dict, minus: dict, *, atol: F, rtol: F) -> dict:
    if atol<0 or rtol<0:raise ValueError('nonnegative tolerances required')
    # These are endpoint enclosures [lower, upper], already including E.
    lp,up=F(plus['lower']),F(plus['upper'])
    lm,um=F(minus['lower']),F(minus['upper'])
    if up<lp or um<lm:raise ValueError('inverted endpoint enclosure')
    lo=max(F(0),lp+lm);hi=up+um
    if hi<lo:raise ValueError('inconsistent two-endpoint certificates')
    gap=hi-lo
    return dict(lower=lo,upper=hi,excess_upper=gap,
        relative_excess_upper=None if lo==0 else gap/lo,
        absolute_tolerance=atol,relative_tolerance=rtol,
        accepted=gap<=atol+rtol*lo)

def canon(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,list):return [canon(v) for v in x]
    if isinstance(x,dict):return {k:canon(v) for k,v in x.items()}
    return x

def main():
    checks=[]
    def check(name,ok,detail=None):
        if not ok:raise AssertionError(name)
        checks.append(dict(name=name,status='PASS',detail=detail))
    # One case per K; each verifies every required derivative identity.
    reflection=[]
    for k in range(1,7):
        n=2*k+1
        a=solve([[F((-ell)**q) for ell in range(1,n+1)] for q in range(n)],[F(1)]*n)
        defects=[sum(a[ell-1]*(-ell)**q for ell in range(1,n+1))-1 for q in range(n)]
        check(f'endpoint_extension_K{k}',all(v==0 for v in defects),{'identities':n})
        reflection.append(dict(K=k,coefficients=a))
    # Response identity, not sampled continuum nonnegativity. Exact values
    # are regression tests of the symbolic identity proved in the manuscript.
    for r in range(1,6):
        h=F(1,32);c=F(1,8)
        ok=True
        for w in (F(1,7),F(1),F(5,2),F(40)):
            direct=(F(0),F(0));prod=(F(1),F(0))
            for ell in range(2*r+1):
                x=c+(ell-r)*h; t=1+x/2
                den=(F(1),w*t);prod=cmul(prod,den)
                atom=F((-1)**ell*comb(2*r,ell),2**(2*r-1))
                direct=cadd(direct,cdiv((atom,F(0)),den))
            num=cpow((F(0),w*h/2),2*r)
            factor=F(factorial(2*r),2**(2*r-1))
            formula=cdiv((num[0]*factor,num[1]*factor),prod)
            ok=ok and direct==formula
        check(f'full_response_identity_r{r}',ok,{'rational_frequencies':4})
        moments=[sum(F((-1)**ell*comb(2*r,ell),2**(2*r-1))*(c+(ell-r)*h)**j for ell in range(2*r+1)) for j in range(2*r)]
        check(f'binomial_moment_cancellation_r{r}',all(x==0 for x in moments),{'moments':2*r})
    for q in (1,2,4,8,12):
        ok=all(F(q)/z-F(q+1)*z/(1+z*z)==(q-z*z)/(z*(1+z*z)) for z in (F(1,3),F(2),F(5)))
        check(f'derivative_maximum_stationarity_order_{q}',ok)
    # Fixed algebra in the summable acquisition-cost bound.
    for q in (F(1,100),F(1,20),F(1,8)):
        check(f'cost_tail_constant_q_{q}',(1+4*q)/(1-4*q)**3<=12 and 1/(1-q)**2<=F(64,49))
    # Numerical budget identity including a nonzero finite-node gap.
    check('three_term_endpoint_budget',F(1,32)+F(1,32)+2*F(1,16)==F(3,16))
    # Flat target: this valid shifted certificate sequence cannot pass pure
    # relative tests with zero true width, but passes a declared absolute floor.
    zero_plus={'lower':'0','upper':'1/64'}
    z=width_bracket(zero_plus,zero_plus,atol=F(0),rtol=F(1,50))
    check('zero_target_rejects_pure_relative_shifted_pair',not z['accepted'] and z['lower']==0)
    z2=width_bracket(zero_plus,zero_plus,atol=F(1,32),rtol=F(1,50))
    check('zero_target_accepts_declared_absolute_floor',z2['accepted'])
    # Generic approximation arithmetic E > 0, not just rational targets E=0.
    E=F(1,100);Pp=F(3);Pm=F(-1);Dp=F(31,10);Dm=F(-9,10)
    pp={'lower':str(Pp-E),'upper':str(Dp+E)}
    mm={'lower':str(Pm-E),'upper':str(Dm+E)}
    br=width_bracket(pp,mm,atol=F(0),rtol=F(1))
    check('proxy_errors_in_both_endpoint_directions',br['lower']==F(99,50) and br['upper']==F(111,50) and br['excess_upper']==F(6,25))
    # Reuse checked original traces. They are not new datasets or solves.
    base=ROOT/'inherited/fs2/EIS_PEM_Math_FS2_Sequential_Cost'
    widths=[]
    for path in sorted((base/'results').glob('sequential_*.json')):
        trace=json.loads(path.read_text())
        for stage in trace['stages']:
            a,b=[e['result'] for e in stage['endpoints']]
            row=width_bracket(a,b,atol=F(0),rtol=F(1,50))
            expected_upper=F(stage['certified_length'])
            check(f'relative_width_{trace["case"]}_stage{stage["stage"]}',row['upper']==expected_upper and row['excess_upper']<=F(stage['numerical_excess_budget']))
            row.update(case=trace['case'],stage=stage['stage'],
                source_input=stage['input'],source_sha256=stage['input_hash'])
            widths.append(row)
    # Explicit rejection cases in this small postprocessing layer.
    negative=[]
    for name,args in [
       ('negative_absolute_tolerance',({'lower':'0','upper':'1'},{'lower':'0','upper':'1'},F(-1),F(1))),
       ('negative_relative_tolerance',({'lower':'0','upper':'1'},{'lower':'0','upper':'1'},F(0),F(-1))),
       ('inverted_positive_endpoint',({'lower':'2','upper':'1'},{'lower':'0','upper':'1'},F(0),F(1))),
       ('inconsistent_width_pair',({'lower':'-3','upper':'-2'},{'lower':'-3','upper':'-2'},F(0),F(1))),
    ]:
        try:width_bracket(args[0],args[1],atol=args[2],rtol=args[3])
        except ValueError: negative.append(dict(name=name,status='REJECTED'))
        else:raise AssertionError(name)
    report=dict(status='PASS',positive_checks=len(checks),rejection_checks=len(negative),
                checks=checks,rejections=negative,reflection_coefficients=reflection,
                relative_widths=widths,
                scope='Exact finite identities and six inherited two-endpoint arithmetic checks; no new endpoint optimization or statistical coverage experiment.')
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'closeout_checks.json').write_text(json.dumps(canon(report),indent=2))
    print(json.dumps({'status':'PASS','positive_checks':len(checks),'rejection_checks':len(negative),'relative_widths':[dict(case=w['case'],stage=w['stage'],ratio=float(w['relative_excess_upper']),accepted=w['accepted']) for w in widths]},indent=2))

if __name__=='__main__':main()
