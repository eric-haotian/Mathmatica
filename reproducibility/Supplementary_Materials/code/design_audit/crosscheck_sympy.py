"""Independent algebra rebuild; does not import the production verifier or exact.py.
SymPy is needed only for this secondary check, not for frozen certificate replay.
"""
from pathlib import Path
import json
import sympy as s
R=Path(__file__).resolve().parent
x,t=s.symbols('x t', real=True)
Q=s.Rational
out=[]
for p in sorted((R/'inputs').glob('*.json')):
 d=json.loads(p.read_text());c=json.loads((R/'certificates'/p.name).read_text())
 atoms=[[ (Q(a),Q(b)) for a,b in d[k] ] for k in ['reference','alternative']]
 th=[[Q(z) for z in d[k]] for k in ['reference_calibration','alternative_calibration']]
 ws=[Q(z) for z in d['frequencies']]
 def obs(mu,tt,w):
  return [tt[0]+sum(a/(1+w*w*(1+xx/2)**2) for xx,a in mu),
          -w*tt[1]+sum(a*w*(1+xx/2)/(1+w*w*(1+xx/2)**2) for xx,a in mu)]
 raw=[]
 for w in ws:
  v0=obs(atoms[0],th[0],w);v1=obs(atoms[1],th[1],w)
  assert v0==v1
  raw+=v0
 assert raw==list(map(Q,d['raw_data']))
 f=sum(Q(a)*x**j for j,a in enumerate(d['target_coefficients']))
 vals=[s.factor(sum(a*f.subs(x,xx) for xx,a in mu)) for mu in atoms]
 assert vals[0]==0 and vals[1]==Q(c['target_gap'])>0
 W=s.prod((x-xx)**2 for xx,a in atoms[0]);assert s.expand(f*Q(c['normalizer'])-W)==0
 ext=ws+[Q(d['additional_frequency'])];D=s.prod(1+w*w*t*t for w in ext)
 lam=list(map(Q,c['extended_design_annihilator_dual']))
 N=sum((lam[2*i]+lam[2*i+1]*w*t)*s.div(D,1+w*w*t*t,t)[0] for i,w in enumerate(ext))
 tar=s.prod(1+w*w for w in ext)*W.subs(x,2*t-2)
 if d['mode']=='unknown':
  tar*=t*t
  assert sum(lam[::2])==0 and sum(w*lam[2*i+1] for i,w in enumerate(ext))==0
 assert s.Poly(s.expand(N-tar),t).is_zero
 assert obs(atoms[0],th[0],ext[-1])!=obs(atoms[1],th[1],ext[-1])
 out.append({'case':p.stem,'status':'PASS','raw_rows':len(raw),'gap':str(vals[1])})
result={'status':'PASS','cases':out,'count':len(out),'scope':'Independent symbolic rebuild of raw equality, target, and extended-design polynomial identity. Not independent human proof review.'}
(R/'results/sympy_crosscheck.json').write_text(json.dumps(result,indent=2))
print(json.dumps({'status':'PASS','count':len(out)}))
