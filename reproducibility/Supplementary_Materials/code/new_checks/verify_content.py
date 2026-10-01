"""Standard-library exact checks for v1.1 content revisions.
Independent of the LP/design proposers and the historical optimization modules.
"""
from fractions import Fraction as F
from pathlib import Path
from math import comb,factorial,isqrt
import json,hashlib,sys

def need(c,msg):
    if not c: raise ValueError(msg)
def q(x): return F(str(x))
def ceilf(x): return -(-x.numerator//x.denominator)
def csqrt(x):
    a=isqrt(x.numerator//x.denominator)
    return a if a*a*x.denominator==x.numerator else a+1

def physical(atoms,theta,ws):
    y=[]
    for w in ws:
        a,b=theta[0],-w*theta[1]
        for x,m in atoms:
            t=1+x/2; z=1+w*w*t*t
            a+=m/z; b+=m*w*t/z
        y.extend([a,b])
    return y

def guard(inp_bytes,cert):
    inp=json.loads(inp_bytes);ws=list(map(q,inp['model']['parameters']))
    need(cert['input_sha256']==hashlib.sha256(inp_bytes).hexdigest(),'input binding')
    need(inp['model']['kind']=='rc' and ws[0]==1 and len(set(ws))==len(ws),'RC baseline')
    y=list(map(q,inp['y']));ds=list(map(q,inp['delta']))
    need(len(y)==len(ds)==2*len(ws) and len(set(ds))==1 and ds[0]>0,'raw dimensions')
    b=ds[0]/3;M=F(13,4)*(y[0]+ds[0]);need(M>0,'positive bound')
    n=max(1,csqrt(4*(M+1)/b));old=ceilf(8*(M+1)/b)+1
    need(q(cert['b'])==b and q(cert['physical_mass_bound'])==M,'derived scales')
    need(cert['cells']==n and cert['nodes']==n+1 and cert['old_lipschitz_nodes']==old,'grid count')
    h=F(2,n);need(q(cert['max_cell'])==h and h*h<=b/(M+1),'quadratic threshold')
    need(M*h*h/4<=b/4,'whole-grid analytic interpolation budget')
    atoms=[(q(x),q(m)) for x,m in cert['atoms']];theta=list(map(q,cert['theta']))
    need(len(theta)==2 and min(theta)>=0,'nonnegative nuisance')
    need(all(-1<=x<=1 and m>0 for x,m in atoms),'positive spectrum')
    need(len(set(x for x,m in atoms))==len(atoms),'unique nodes')
    need(all(((x+1)*n/2).denominator==1 for x,m in atoms),'anchor on declared grid')
    pred=physical(atoms,theta,ws);margin=min(d-abs(v-z) for d,v,z in zip(ds,pred,y))
    need(margin>0 and margin==q(cert['margin']),'strict raw anchor')
    need(sum(m for _,m in atoms)<=M,'mass bound')
    return {'input':cert['input_name'],'nodes':n+1,'old_nodes':old,
            'margin_over_b':str(margin/b),'support_count':len(atoms),
            'old_to_new_count_ratio':str(F(old,n+1))}

def mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out

def design(cert):
    r=cert['r'];need(type(r) is int and 1<=r<=4,'r range')
    n=2*r;h=q(cert['h']);c=q(cert['center']);need(h>0,'positive h')
    xs=[c+(i-r)*h for i in range(n+1)];ts=[1+x/2 for x in xs]
    need(all(-1<=x<=1 for x in xs),'compact support')
    ref=[(q(x),q(m)) for x,m in cert['reference']];alt=[(q(x),q(m)) for x,m in cert['alternative']]
    expect0=[(xs[i],F(comb(n,i),2**(n-1))) for i in range(n+1) if i%2]
    expect1=[(xs[i],F(comb(n,i),2**(n-1))) for i in range(n+1) if i%2==0]
    need(ref==expect0 and alt==expect1,'positive pair structure')
    need(sum(m for x,m in ref)==sum(m for x,m in alt)==1,'unit mass')
    need(all(sum(m*x**j for x,m in ref)==sum(m*x**j for x,m in alt) for j in range(n)),'moment cancellation')
    # Rebuild the complex-frequency numerator as a polynomial in z=i*omega.
    coeff=[F(0)]*(n+1)
    for i in range(n+1):
        p=[F(1)]
        for j,t in enumerate(ts):
            if j!=i:p=mul(p,[F(1),t])
        signed=F((-1)**i*comb(n,i),2**(n-1))
        for k,a in enumerate(p):coeff[k]+=signed*a
    lead=F(factorial(n),2**(n-1))*(h/2)**n
    need(coeff==[F(0)]*n+[lead],'continuum numerator identity')
    def Q(v):
        z=F(1)
        for t in ts:z*=1+t*t*v
        return z
    L=q(cert['v_lo']);U=q(cert['v_hi'])
    A=lambda v:sum(t*t*v/(1+t*t*v) for t in ts)
    need(0<L<U and A(L)<=n<=A(U),'maximizer bracket')
    mid=(L+U)/2;low=lead*lead*mid**n/Q(mid);high=lead*lead*U**n/Q(L)
    need(low==q(cert['max_gap_square_lower']) and high==q(cert['max_gap_square_upper']),'full norm bounds')
    need(0<low<=high and high/low-1<F(1,2**90),'norm accuracy')
    result={'r':r,'h':str(h),'efficiencies':{}}
    for name,expected_ws in [('baseline',list(map(F,range(1,7)))),('high_band',[F(16*i) for i in range(1,7)])]:
        d=cert['designs'][name];ws=list(map(q,d['frequencies']))
        need(ws==expected_ws,'fixed, unchanged design')
        p=physical(ref,[F(0),F(0)],ws);a=physical(alt,[F(0),F(0)],ws)
        vals=[(a[2*j]-p[2*j])**2+(a[2*j+1]-p[2*j+1])**2 for j in range(6)]
        need(vals==list(map(q,d['point_gap_squares'])),'raw pair differences')
        mean=sum(vals)/6;need(mean==q(d['mean_gap_square']),'balanced effort average')
        lo=mean/high;hi=mean/low
        need(lo==q(d['oracle_efficiency_lower']) and hi==q(d['oracle_efficiency_upper']),'efficiency enclosure')
        need(0<lo<=hi<=1,'efficiency range')
        t0=1+c/2
        maxlimit=F(n**n,(n+1)**(n+1))/t0**(2*n)
        avg_limit=sum(w**(2*n)/(1+w*w*t0*t0)**(n+1) for w in ws)/6
        need(avg_limit/maxlimit==q(d['limit_efficiency']),'nonzero limiting constant')
        result['efficiencies'][name]=[str(lo),str(hi)]
    return result

# A finite-budget Bernstein checker used to exhibit and avoid an inner-loop stall.
def bernstein_on_unit(p):
    n=len(p)-1
    return [sum(p[k]*F(comb(j,k),comb(n,k)) for k in range(j+1)) for j in range(n+1)]
def casteljau_half(b):
    left=[b[0]];right=[b[-1]];row=b[:]
    while len(row)>1:
        row=[(a+c)/2 for a,c in zip(row,row[1:])]
        left.append(row[0]);right.append(row[-1])
    return left,right[::-1]
def bounded_positive(p,depth):
    # p is in the power basis on [0,1]. Each call has a finite tree budget.
    need(type(depth)is int and depth>=0,'depth')
    stack=[(bernstein_on_unit(list(map(q,p))),0)];pending=False;visits=0
    while stack:
        b,d=stack.pop();visits+=1
        if b[0]<0 or b[-1]<0:return {'status':'NEGATIVE_WITNESS','visits':visits}
        if min(b)>=0:continue
        if d==depth:pending=True;continue
        a,c=casteljau_half(b);stack.extend([(a,d+1),(c,d+1)])
    return {'status':'UNRESOLVED' if pending else 'CERTIFIED','visits':visits}

def mechanism_tests():
    out=[]
    p=[F(1,9),F(-2,3),F(1)]
    for d in [0,2,4,8,16]:
        r=bounded_positive(p,d);need(r['status']=='UNRESOLVED','interior non-dyadic double zero not accepted')
        out.append({'name':f'zero_contact_depth_{d}',**r})
    pshift=p[:];pshift[0]+=F(1,2**16)
    first=next((d for d in range(20) if bounded_positive(pshift,d)['status']=='CERTIFIED'),None)
    need(first is not None,'positive shift eventual certification')
    out.append({'name':'positive_shift_finite_completion','first_success_depth':first})
    need(bounded_positive([F(-1,4),F(0),F(1)],6)['status']=='NEGATIVE_WITNESS','negative polynomial rejection')
    out.append({'name':'negative_polynomial_rejected','passed':True})
    # A verified, NONOPTIMAL finite LP pair in a mass-box inverse problem.
    tau=F(1,1024);M=F(9,8);P=M-tau/64;gamma=F(1);delta_lp=M-P
    a=tau/(32*M);D=(gamma+a)*M
    need(delta_lp>0 and delta_lp<=tau/32 and D-P==3*tau/64,'nonzero finite gap budget')
    need(D-P<tau/2 and P>F(7,8),'valid finite witness')
    out.append({'name':'nonzero_node_gap_is_allowed','node_gap_over_tau':str(delta_lp/tau),
                'full_gap_over_tau':str((D-P)/tau)})
    need(F(1,32)+F(1,32)+2*F(1,16)==F(3,16),'end-to-end exact budget')
    out.append({'name':'proxy_node_continuum_budget','bound_over_tau':'3/16'})
    # Error-spending / expected-cost summability accounting.
    for alpha in [F(1,4),F(1,20),F(1,100)]:
        risk=[alpha/F(2**j) for j in range(1,21)]
        need(sum(risk)<alpha and 4*(alpha/2)<1,'all-stage risk and tail')
    out.append({'name':'risk_spending_and_cost_tail','passed':True})
    return out

def main():
    if sys.argv[1]=='guard':r=guard(Path(sys.argv[2]).read_bytes(),json.loads(Path(sys.argv[3]).read_text()))
    elif sys.argv[1]=='design':r=design(json.loads(Path(sys.argv[2]).read_text()))
    elif sys.argv[1]=='mechanisms':r=mechanism_tests()
    else:raise ValueError('unknown task')
    print(json.dumps(r))
if __name__=='__main__':main()
