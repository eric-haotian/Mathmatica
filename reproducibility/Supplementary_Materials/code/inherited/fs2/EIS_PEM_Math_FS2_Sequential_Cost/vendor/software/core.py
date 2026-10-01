"""R15 exact model/proxy/certificate primitives. Verification uses stdlib only.
All numbers are rational strings or integers. Polynomial coefficients ascend.
Historical algebra primitives are imported unchanged, with explicit provenance.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
from functools import lru_cache
import sys,json,hashlib
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'inherited'))
from r13_core import (frac,trim,add,sub,scale,mul,power,deriv,val,compose_affine,
                      divide,solve,bernstein,dot,l1,jsonable,model as oldmodel)

def need(ok,msg):
    if not ok: raise ValueError(msg)
def save(path,obj):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(jsonable(obj),ensure_ascii=False,indent=2)+'\n')
def digest(obj):
    return hashlib.sha256(json.dumps(jsonable(obj),sort_keys=True,separators=(',',':')).encode()).hexdigest()
def load(path): return json.loads(Path(path).read_text())
def mat_det(A):
    A=[list(map(frac,row)) for row in A]; n=len(A); out=F(1)
    for j in range(n):
        k=next((i for i in range(j,n) if A[i][j]),None)
        if k is None:return F(0)
        if k!=j:A[k],A[j]=A[j],A[k];out=-out
        z=A[j][j];out*=z
        for i in range(j+1,n):
            q=A[i][j]/z
            for l in range(j+1,n):A[i][l]-=q*A[j][l]
    return out

@lru_cache(None)
def get_model(kind,parameters=()):
    parameters=tuple(frac(x) for x in parameters)
    if kind=='rc':
        ws=parameters or tuple(map(F,(1,2,3,4)));mo=oldmodel(ws)
        den=scale(power([F(1),F(1,2)],2),mo['D0'])
        nums=[compose_affine(p,F(1),F(1,2)) for p in mo['numerators']]
        rhon=compose_affine(mo['D'],F(1),F(1,2));rhod=den
        B=[]
        for w in ws:B.extend([[F(1),F(0)],[F(0),-w]])
        lam=mo['lambda'][:5]
    elif kind=='resolvent':
        al=parameters or tuple(map(F,(2,3,4,5,6,7,8)))
        need(len(al)>=7 and len(set(al))==len(al) and min(al)>1,'resolvent poles')
        D=[F(1)]
        for a in al:D=mul(D,[a,F(1)])
        d0=val(D,0);nums=[]
        for a in al:
            p,r=divide(D,[a,F(1)]);need(r==[0],'resolvent factor');nums.append(p)
        den=[d0];rhon=D;rhod=den;B=[[F(1),a] for a in al]
        C=[[p[j] if j<len(p) else F(0) for p in nums] for j in range(len(al))]
        lam=[]
        for j in range(5):
            rhs=[F(0)]*len(al);rhs[j]=d0;lam.append(solve(C,rhs))
    elif kind=='moments':
        den=[F(1)];nums=[[F(0)]*i+[F(1)] for i in range(5)]
        rhon=rhod=[F(1)];B=[[] for _ in range(5)]
        lam=[[F(i==j) for i in range(5)] for j in range(5)]
    elif kind=='mass_only':
        den=[F(1)];nums=[[F(1)]];rhon=rhod=[F(1)];B=[[]];lam=[[F(1)]]
    else:raise ValueError('unknown model')
    m=len(nums);k=len(B[0])
    for j,L in enumerate(lam):
        out=[F(0)]
        for c,p in zip(L,nums):out=add(out,scale(p,c))
        need(out==mul(den,[F(0)]*j+[F(1)]),'polynomial reproduction')
        need(all(dot(L,[B[i][l] for i in range(m)])==0 for l in range(k)),'calibration annihilation')
    need(min(bernstein(den,-1,1))>0,'positive common denominator')
    return dict(kind=kind,parameters=parameters,den=den,nums=nums,B=B,lam=lam,rhon=rhon,rhod=rhod,m=m,k=k)

def from_spec(s):
    need(set(s)<= {'kind','parameters'},'unsupported model metadata')
    return get_model(s['kind'],tuple(s.get('parameters',[])))
def obs(mo,x):return [val(p,x)/val(mo['den'],x) for p in mo['nums']]
def raw(mo,atoms,theta):
    need(len(theta)==mo['k'],'calibration dimensions')
    out=[F(0)]*mo['m']
    for x,w in atoms:
        need(-1<=x<=1 and w>=0,'nonnegative measure on declared interval')
        out=[a+w*b for a,b in zip(out,obs(mo,x))]
    need(all(t>=0 for t in theta),'nonnegative calibration')
    return [x+dot(row,theta) for x,row in zip(out,mo['B'])]
def lambda_poly(mo,p):
    need(len(p)<=len(mo['lam']),'degree outside observed moment coordinates')
    return [sum((c*mo['lam'][j][i] for j,c in enumerate(p)),F(0)) for i in range(mo['m'])]
def support(g,y,de):return dot(g,y)+dot(list(map(abs,g)),de)
def data_context(inp):
    need(set(inp)=={'format','model','y','delta','target','tau'},'raw input must not contain truth/reference/anchor fields')
    need(inp['format']=='R15_raw_data_v1','format')
    mo=from_spec(inp['model']);y=list(map(frac,inp['y']));de=list(map(frac,inp['delta']));tau=frac(inp['tau'])
    need(len(y)==len(de)==mo['m'] and min(de)>0 and tau>0,'positive raw radii and tolerance')
    M=support(mo['lam'][0],y,de);need(M>0,'positive mass budget')
    return mo,y,de,tau,M

def localization(mo,y,de):
    need(len(mo['lam'])>=5,'five observed moments needed')
    mm=[dot(z,y) for z in mo['lam'][:5]];rr=[dot(list(map(abs,z)),de) for z in mo['lam'][:5]]
    eta=max(sum((rr[i+j] for j in range(3)),F(0)) for i in range(3))
    need(eta>0,'positive moment uncertainty')
    Hb=[[mm[i+j]+(2*eta if i==j else 0) for j in range(3)] for i in range(3)]
    need(all(mat_det([row[:j] for row in Hb[:j]])>0 for j in range(1,4)),'positive inflated moment matrix')
    # Exact convex quadratic minimization on b in [-1/4,1/4], a in [-1/2,1/2].
    clip=lambda z,lo,hi:max(lo,min(hi,z))
    cand=[]
    ba=solve([row[:2] for row in Hb[:2]],[-Hb[0][2],-Hb[1][2]])
    if -F(1,4)<=ba[0]<=F(1,4) and -F(1,2)<=ba[1]<=F(1,2):cand.append(ba)
    for b in [-F(1,4),F(1,4)]:
        a=clip(-(Hb[1][0]*b+Hb[1][2])/Hb[1][1],-F(1,2),F(1,2));cand.append([b,a])
    for a in [-F(1,2),F(1,2)]:
        b=clip(-(Hb[0][1]*a+Hb[0][2])/Hb[0][0],-F(1,4),F(1,4));cand.append([b,a])
    qform=lambda z:dot(z+[F(1)],[dot(row,z+[F(1)]) for row in Hb])
    ba=min(cand,key=qform);Qbar=qform(ba);T=ba+[F(1)];W=mul(T,T)
    Qraw=support(lambda_poly(mo,W),y,de);need(Qraw>=0 and Qbar>=Qraw,'valid polynomial localization')
    return dict(moments=mm,radii=rr,eta=eta,Hbar=Hb,T=T,W=W,Qbar=Qbar,Q=Qraw,candidates=len(cand))

def proxy(mo,y,de,tau,M,target):
    kind=target['kind'];E=F(0);loc=None
    if kind=='poly':
        need(set(target)=={'kind','coefficients'},'poly target schema');p=trim(target['coefficients']);deg=len(p)-1
    elif kind=='resolvent':
        need(set(target)=={'kind','beta','scale'},'resolvent target schema')
        beta=frac(target['beta']);S=frac(target['scale']);need(beta>=2 and S>0,'resolvent target bounds')
        loc=localization(mo,y,de);W=loc['W'];Wb=val(W,beta);need(Wb>0,'no denominator root at target pole')
        H,rem=divide(sub([F(1)],scale(W,1/Wb)),[beta,F(-1)]);need(rem==[0],'root-free quotient identity')
        n=0
        while loc['Q']/(Wb*S*(beta-1)*beta**(n+1))>tau/16:n+=1
        tn=[F(1)/beta**(k+1) for k in range(n+1)]
        p=scale(add(H,scale(mul(W,tn),1/Wb)),1/S)
        E=loc['Q']/(Wb*S*(beta-1)*beta**(n+1));deg=len(p)-1
    elif kind=='exp':
        need(set(target)=={'kind','scale'},'exp target schema');S=frac(target['scale']);need(S>0,'scale')
        n=0
        while 3*M/(S*factorial(n+1))>tau/16:n+=1
        p=[F(1,S*factorial(k)) for k in range(n+1)];E=3*M/(S*factorial(n+1));deg=len(p)-1
    else:raise ValueError('unsupported target')
    need(E<=tau/16,'target approximation budget')
    degree=max(max(len(x)-1 for x in mo['nums']),len(mo['den'])-1+deg)
    return dict(poly=p,E=E,E_ar=F(0),degree=deg,moment_scalars=degree+1,localization=loc)

def positive(p,maxdepth=48):
    stack=[(F(-1),F(1),0)];leaves=deep=0
    while stack:
        a,b,j=stack.pop();bb=bernstein(p,a,b)
        if min(bb)>=0:leaves+=1;deep=max(deep,j);continue
        c=(a+b)/2
        for x in (a,c,b):
            if val(p,x)<0:raise ValueError('negative polynomial at '+str(x))
        if j>=maxdepth:raise ValueError('positivity unresolved at depth cap')
        stack.extend([(a,c,j+1),(c,b,j+1)])
    return dict(leaves=leaves,max_depth=deep)

def dual_numerator(mo,g,p,sign):
    out=[F(0)]
    for c,z in zip(g,mo['nums']):out=add(out,scale(z,c))
    return sub(out,scale(mul(mo['den'],p),sign))

def derivative_bound(mo,order):
    D=mo['den'];Dp=deriv(D);bounds=[]
    for N0 in mo['nums']:
        N=N0
        for j in range(order):N=sub(mul(deriv(N),D),scale(mul(N,Dp),j+1))
        bd=F(0)
        for i in range(32):
            a=-1+F(i,16);b=a+F(1,16);dl=min(bernstein(D,a,b));need(dl>0,'denominator bound')
            bd=max(bd,max(map(abs,bernstein(N,a,b)))/dl**(order+1))
        bounds.append(bd)
    return bounds

def verify(inp,cert):
    mo,y,de,tau,M=data_context(inp)
    need(cert['format']=='R15_endpoint_v1' and cert['input_hash']==digest(inp),'input binding')
    sg=cert['sign'];need(sg in[-1,1],'endpoint sign')
    pr=proxy(mo,y,de,tau,M,inp['target'])
    need(cert['proxy_hash']==digest(pr),'proxy binding')
    atoms=[[frac(x),frac(w)] for x,w in cert['atoms']];th=list(map(frac,cert['theta']))
    z=raw(mo,atoms,th);need(all(abs(v-b)<=d for v,b,d in zip(z,y,de)),'raw primal box')
    aa=[[frac(x),frac(w)] for x,w in cert['anchor_atoms']];at=list(map(frac,cert['anchor_theta']))
    az=raw(mo,aa,at);margin=min(d-abs(v-b) for v,b,d in zip(az,y,de));need(margin>0,'strict data anchor')
    nodes=sorted(set(map(frac,cert['nodes'])));need(nodes[0]==-1 and nodes[-1]==1 and all(x in nodes for x,w in atoms),'finite node domain')
    gamma=list(map(frac,cert['gamma']));shift=frac(cert['shift']);need(len(gamma)==mo['m'] and shift>=0,'dual dimensions')
    need(all(dot(gamma,[row[k] for row in mo['B']])>=0 for k in range(mo['k'])),'dual calibration cone')
    need(all(dot(gamma,obs(mo,x))>=sg*val(pr['poly'],x) for x in nodes),'finite dual inequalities')
    corrected=[a+shift*b for a,b in zip(gamma,mo['lam'][0])]
    pp=dual_numerator(mo,corrected,pr['poly'],sg);need(isinstance(cert['maxdepth'],int) and cert['maxdepth']>=1,'positivity proof depth');proof=positive(pp,maxdepth=cert['maxdepth'])
    P=sum((w*sg*val(pr['poly'],x) for x,w in atoms),F(0));D0=support(gamma,y,de);D=support(corrected,y,de)
    need(D0>=P,'finite weak duality')
    # Triangle support bound suffices even if J(corrected)-J(gamma) is negative.
    need(D<=D0+shift*M,'mass repair support bound')
    E=pr['E'];lo=P-E;hi=D+E;gap=hi-lo
    need(gap>=0 and gap<=tau/2,'complete endpoint stopping budget')
    return dict(sign=sg,lower=lo,upper=hi,gap=gap,gap_over_tau=gap/tau,
                finite_gap=D0-P,continuous_dual_repair_bound=shift*M,
                interpolation_bound=E,coefficient_bound=pr['E_ar'],
                strict_anchor_margin=margin,mass_bound=M,proxy_degree=pr['degree'],
                moment_scalars=pr['moment_scalars'],primal_atoms=len(atoms),
                finite_nodes=len(nodes),positivity=proof,
                posterior_Q=None if pr['localization'] is None else pr['localization']['Q'])
