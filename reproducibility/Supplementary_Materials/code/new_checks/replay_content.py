"""Cold replay of new exact certificates, plus intentional invalid-certificate tests."""
from pathlib import Path
import subprocess,sys,tempfile,shutil,json,time,copy
import verify_content as V
ROOT=Path(__file__).resolve().parent

def run_one(kind,cert):
    with tempfile.TemporaryDirectory(prefix='sinum_v11_') as tmp:
        d=Path(tmp);shutil.copy2(ROOT/'verify_content.py',d/'verify_content.py')
        shutil.copy2(cert,d/'certificate.json')
        args=[sys.executable,'-I',str(d/'verify_content.py'),kind]
        if kind=='guard':
            obj=json.loads(cert.read_text());src=ROOT/'inputs'/obj['input_name']
            shutil.copy2(src,d/'input.json');args.extend([str(d/'input.json'),str(d/'certificate.json')])
        else:args.append(str(d/'certificate.json'))
        proc=subprocess.run(args,check=True,capture_output=True,text=True,timeout=40)
        return json.loads(proc.stdout)

def main():
    t0=time.perf_counter();results=[]
    for cert in sorted((ROOT/'certificates').glob('*.json')):
        kind='guard' if cert.name.startswith('guard_') else 'design'
        result=run_one(kind,cert);results.append({'certificate':cert.name,'kind':kind,'result':result})
    mechanisms=V.mechanism_tests();negative=[]
    gc=json.loads(next((ROOT/'certificates').glob('guard_*')).read_text())
    ib=(ROOT/'inputs'/gc['input_name']).read_bytes()
    mutants=[]
    a=copy.deepcopy(gc);a['atoms'][0][1]='-1';mutants.append(('negative_guard_mass',a))
    a=copy.deepcopy(gc);a['theta'][0]='-1';mutants.append(('negative_guard_nuisance',a))
    a=copy.deepcopy(gc);a['margin']='1';mutants.append(('inflated_guard_margin',a))
    a=copy.deepcopy(gc);a['cells']=1;mutants.append(('undersized_guard_grid',a))
    for name,a in mutants:
        try:V.guard(ib,a)
        except (ValueError,AssertionError):negative.append({'name':name,'rejected':True})
        else:raise AssertionError(name)
    dc=json.loads(next((ROOT/'certificates').glob('information_*')).read_text())
    for name in ['fake_oracle_efficiency','changed_same_budget_design']:
        a=copy.deepcopy(dc)
        if name.startswith('fake'):a['designs']['baseline']['oracle_efficiency_lower']='1'
        else:a['designs']['baseline']['frequencies'][0]='2'
        try:V.design(a)
        except (ValueError,AssertionError):negative.append({'name':name,'rejected':True})
        else:raise AssertionError(name)
    out={'fresh_process_tasks':len(results),'guard_tasks':sum(x['kind']=='guard' for x in results),
         'design_tasks':sum(x['kind']=='design' for x in results),'results':results,
         'mechanisms':mechanisms,'invalid_certificate_tests':negative,
         'seconds':time.perf_counter()-t0,
         'scope':'Exact finite certificate checks, not theorem proof or Monte Carlo coverage.'}
    (ROOT/'results/replay.json').write_text(json.dumps(out,indent=2))
    print(json.dumps({k:v for k,v in out.items() if k not in ['results','mechanisms','invalid_certificate_tests']}))
if __name__=='__main__':main()
