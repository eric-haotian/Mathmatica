"""Rational LP certification with finite exact Bland-simplex fallback.
Maximize c*x subject to A*x<=b, x>=0. Caller guarantees bounded objective.
The fast numerical stage is ONLY a proposer; no floating result is accepted.
A second exhaustive mode enumerates all standard-form bases for small tests. It is deliberately
not represented as polynomial-time or suitable for large instances.
"""
from fractions import Fraction as F
from itertools import combinations,islice
from core import frac,solve,dot,need

def check(A,b,c,x,u):
    if min(x,default=0)<0 or min(u,default=0)<0:return False
    if any(dot(row,x)>bb for row,bb in zip(A,b)):return False
    if any(dot(u,[row[j] for row in A])<cc for j,cc in enumerate(c)):return False
    return dot(c,x)==dot(b,u)

def fast_candidate(A,b,c):
    import numpy as np
    from scipy.optimize import linprog
    Af=np.array(A,dtype=float);bf=np.array(b,dtype=float);cf=np.array(c,dtype=float)
    scale=max(float(np.max(np.abs(cf))),1e-10)
    ans=linprog(-cf/scale,A_ub=Af,b_ub=bf,bounds=(0,None),method='highs-ds',
                options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    if not ans.success:return None,dict(float_status=int(ans.status),message=ans.message)
    r=len(b);n=len(c);res=bf-Af@ans.x;du=-ans.ineqlin.marginals*scale;tested=0
    for threshold in [1e-8,1e-11,1e-14,0.0]:
        J=[j for j,z in enumerate(ans.x) if z>threshold];k=len(J)
        if k==0:
            x=[F(0)]*n;u=[F(0)]*r
            if check(A,b,c,x,u):return (x,u),dict(method='exact_zero',tested=tested)
            continue
        if k>r:continue
        E1=[i for i,v in enumerate(du) if v>1e-12]
        E2=sorted(range(r),key=lambda i:abs(res[i]))
        pool=list(dict.fromkeys(E1+E2[:min(r,k+3)]))
        trials=[tuple(E1)] if len(E1)==k else []
        trials+=list(islice(combinations(pool,k),2000))
        for E in trials:
            tested+=1
            try:
                M=[[A[i][j] for j in J] for i in E];sol=solve(M,[b[i] for i in E])
                uv=solve(list(map(list,zip(*M))),[c[j] for j in J])
            except ValueError:continue
            if min(sol)<0 or min(uv)<0:continue
            x=[F(0)]*n;u=[F(0)]*r
            for j,v in zip(J,sol):x[j]=v
            for i,v in zip(E,uv):u[i]=v
            if check(A,b,c,x,u):return (x,u),dict(method='floating_proposal_exact_basis',tested=tested,float_status=0)
    return None,dict(method='unresolved_basis_proposal',tested=tested,float_status=0)

def exhaustive(A,b,c):
    r=len(b);n=len(c);cols=[[row[j] for row in A] for j in range(n)]
    cols += [[F(i==j) for i in range(r)] for j in range(r)]
    costs=c+[F(0)]*r;count=0;feasible=False
    for J in combinations(range(n+r),r):
        count+=1;M=[[cols[j][i] for j in J] for i in range(r)]
        try:sol=solve(M,b)
        except ValueError:continue
        if min(sol)<0:continue
        feasible=True
        u=solve(list(map(list,zip(*M))),[costs[j] for j in J])
        x=[F(0)]*n
        for j,v in zip(J,sol):
            if j<n:x[j]=v
        if check(A,b,c,x,u):return (x,u),dict(method='exhaustive_rational_bases',tested=count)
    if not feasible:return None,dict(method='exhaustive_rational_bases',tested=count,infeasible=True)
    raise ValueError('no bounded optimum basis: caller assumptions violated')

def solve_lp(A,b,c,complete=False,force_exact=False,warm=None):
    A=[list(map(frac,row)) for row in A];b=list(map(frac,b));c=list(map(frac,c))
    need(len(A)==len(b) and all(len(row)==len(c) for row in A),'LP dimensions')
    if not force_exact:
        try:
            result,stats=fast_candidate(A,b,c)
        except Exception as ex:
            result,stats=None,dict(method='proposal_exception',message=str(ex))
        if result is not None:return result,stats
        # Exact simplex is the production fallback; exhaustive bases remain testable.
        if warm is not None:
            try: result,st=warm_simplex(A,b,c,warm)
            except ValueError: result,st=bland_simplex(A,b,c)
        else:result,st=bland_simplex(A,b,c)
        st['failed_float_proposal']=stats
        return result,st
    return bland_simplex(A,b,c)


def bland_simplex(A,b,c):
    """Exact two-phase simplex. Bland pivoting; no floating tolerance.
    Slack-column reduced costs reconstruct the dual in ORIGINAL inequality signs.
    This is a classical finite LP algorithm, not an R15 novelty claim.
    """
    R=len(b);N=len(c);neg=[i for i,bb in enumerate(b) if bb<0]
    nc=N+R+len(neg);rows=[];rhs=[];base=[];art={i:N+R+k for k,i in enumerate(neg)}
    for i,(row,bb) in enumerate(zip(A,b)):
        sg=F(1) if bb>=0 else F(-1)
        rr=[sg*x for x in row]+[F(0)]*(nc-N);rr[N+i]=sg
        if i in art:rr[art[i]]=F(1);base.append(art[i])
        else:base.append(N+i)
        rows.append(rr);rhs.append(abs(bb))
    pivots=0
    def pivot(i,j):
        nonlocal pivots
        pp=rows[i][j];need(pp!=0,'simplex nonzero pivot')
        rows[i]=[v/pp for v in rows[i]];rhs[i]/=pp
        for k in range(len(rows)):
            if k==i:continue
            f=rows[k][j]
            if f:
                rows[k]=[a-f*v for a,v in zip(rows[k],rows[i])];rhs[k]-=f*rhs[i]
        base[i]=j;pivots+=1
    def optimize(cost):
        while True:
            reduced=cost[:]
            for i,j in enumerate(base):
                if cost[j]:reduced=[x-cost[j]*v for x,v in zip(reduced,rows[i])]
            entering=next((j for j,z in enumerate(reduced) if z>0),None)
            if entering is None:return reduced
            candidates=[(rhs[i]/row[entering],base[i],i) for i,row in enumerate(rows) if row[entering]>0]
            if not candidates:raise ValueError('unbounded exact LP')
            leaving=min(candidates)[2];pivot(leaving,entering)
    phase1=[F(0)]*(N+R)+[F(-1)]*len(neg)
    optimize(phase1)
    objective=dot([phase1[j] for j in base],rhs)
    if objective<0:return None,dict(method='exact_bland_simplex',pivots=pivots,infeasible=True)
    need(objective==0,'phase I objective')
    for i in range(R):
        if base[i]>=N+R:
            need(rhs[i]==0,'zero artificial basic')
            used=set(base)
            j=next((j for j in range(N+R) if j not in used and rows[i][j]),None)
            need(j is not None,'full-row-rank slack formulation');pivot(i,j)
    rows[:]=[row[:N+R] for row in rows]
    cost=c+[F(0)]*R;reduced=optimize(cost)
    x=[F(0)]*N
    for i,j in enumerate(base):
        if j<N:x[j]=rhs[i]
    u=[-reduced[N+i] for i in range(R)]
    need(check(A,b,c,x,u),'exact simplex primal dual self-check')
    return (x,u),dict(method='exact_bland_simplex',pivots=pivots)


def warm_simplex(A,b,c,warm):
    """Exact primal simplex reoptimization from a previous feasible vertex.
    New proposal nodes enter with zero mass. This accelerates, not changes, the LP.
    """
    R=len(b);N=len(c);slack=[bb-dot(row,warm) for row,bb in zip(A,b)]
    need(min(warm,default=0)>=0 and min(slack)>=0,'warm start primal feasibility')
    cols=[[row[j] for row in A] for j in range(N)]+[[F(i==j) for i in range(R)] for j in range(R)]
    values=warm+slack;positive_cols=[j for j,v in enumerate(values) if v>0]
    order=positive_cols+[j for j in range(N,N+R) if j not in positive_cols]+[j for j in range(N) if j not in positive_cols]
    chosen=[];ec=[]
    for j in order:
        v=cols[j][:]
        for k,ee in sorted(ec):
            q=v[k]
            if q:v=[x-q*y for x,y in zip(v,ee)]
        k=next((i for i,x in enumerate(v) if x),None)
        if k is None:
            need(j not in positive_cols,'warm point is not a basic feasible point')
            continue
        z=v[k];v=[x/z for x in v];ec.append((k,v));chosen.append(j)
        if len(chosen)==R:break
    need(len(chosen)==R and all(j in chosen for j in positive_cols),'full warm basis')
    aug=[[cols[j][i] for j in chosen]+[cols[j][i] for j in range(N+R)]+[b[i]] for i in range(R)]
    for j in range(R):
        k=next(i for i in range(j,R) if aug[i][j]);aug[k],aug[j]=aug[j],aug[k]
        z=aug[j][j];aug[j]=[x/z for x in aug[j]]
        for i in range(R):
            if i!=j and aug[i][j]:
                z=aug[i][j];aug[i]=[x-z*y for x,y in zip(aug[i],aug[j])]
    rows=[row[R:-1] for row in aug];rhs=[row[-1] for row in aug];base=chosen;cost=c+[F(0)]*R;pivots=0
    need(min(rhs)>=0,'warm basis remains feasible')
    while True:
        red=cost[:]
        for i,j in enumerate(base):
            if cost[j]:red=[x-cost[j]*v for x,v in zip(red,rows[i])]
        entering=next((j for j,z in enumerate(red) if z>0),None)
        if entering is None:break
        candidates=[(rhs[i]/row[entering],base[i],i) for i,row in enumerate(rows) if row[entering]>0]
        need(bool(candidates),'bounded warm LP');i=min(candidates)[2];pp=rows[i][entering]
        rows[i]=[x/pp for x in rows[i]];rhs[i]/=pp
        for k in range(R):
            if k!=i and rows[k][entering]:
                q=rows[k][entering];rows[k]=[a-q*v for a,v in zip(rows[k],rows[i])];rhs[k]-=q*rhs[i]
        base[i]=entering;pivots+=1
    x=[F(0)]*N
    for i,j in enumerate(base):
        if j<N:x[j]=rhs[i]
    u=[-red[N+i] for i in range(R)];need(check(A,b,c,x,u),'warm exact LP certificate')
    return (x,u),dict(method='exact_warm_bland_simplex',pivots=pivots)
