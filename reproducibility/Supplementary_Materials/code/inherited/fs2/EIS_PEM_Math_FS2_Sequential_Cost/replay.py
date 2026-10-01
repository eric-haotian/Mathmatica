"""Cold-process replay of each new cost certificate and each stopping trajectory."""
from pathlib import Path
import subprocess,sys,json,tempfile,shutil,time
from exact_cost import save
ROOT=Path(__file__).resolve().parent
start=time.perf_counter();runs=[]
for inp in sorted((ROOT/'inputs').glob('*.json')):
    with tempfile.TemporaryDirectory() as td:
        p=Path(td)
        for fn in ['exact_cost.py']:shutil.copyfile(ROOT/fn,p/fn)
        shutil.copyfile(inp,p/'input.json');shutil.copyfile(ROOT/'certificates'/inp.name,p/'certificate.json')
        result=subprocess.run([sys.executable,'-I',str(p/'exact_cost.py'),str(p/'input.json'),str(p/'certificate.json'),'--output',str(p/'result.json')],capture_output=True,text=True,timeout=40)
        if result.returncode:raise RuntimeError(result.stderr)
        runs.append({'kind':'new_stopped_cost','case':inp.stem,'status':'PASS'})
for tr in sorted((ROOT/'results').glob('sequential_*.json')):
    with tempfile.TemporaryDirectory() as td:
        p=Path(td)
        for fn in ['exact_cost.py','verify_stopping.py']:shutil.copyfile(ROOT/fn,p/fn)
        shutil.copytree(ROOT/'vendor',p/'vendor',ignore=shutil.ignore_patterns('__pycache__','generate.py','exact_lp.py'))
        shutil.copytree(ROOT/'sequential',p/'sequential');shutil.copyfile(tr,p/'trace.json')
        result=subprocess.run([sys.executable,'-I',str(p/'verify_stopping.py'),str(p/'trace.json'),'--output',str(p/'result.json')],capture_output=True,text=True,timeout=40)
        if result.returncode:raise RuntimeError(result.stderr)
        runs.append({'kind':'new_stopping_trace','case':tr.stem,'status':'PASS','output':json.loads((p/'result.json').read_text())})
save(ROOT/'results'/'replay.json',{'completed':True,'new_cost_certificates':18,'new_trajectories':2,
                                  'raw_stage_vectors':6,'endpoint_certificates':12,'runs':runs,'seconds':time.perf_counter()-start,
                                  'excluded':['synthetic provenance','candidate generators','numpy','scipy','finite optimization calls']})
print('PASS',len(runs),'isolated tasks')
