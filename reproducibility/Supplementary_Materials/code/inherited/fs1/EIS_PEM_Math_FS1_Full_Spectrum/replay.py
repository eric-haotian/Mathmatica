"""Fresh-process replay of frozen certificates, with no generator/optimizer present."""
from pathlib import Path
import tempfile,subprocess,shutil,json,time,sys
from exact import dump,digest
ROOT=Path(__file__).resolve().parent

def main():
    out=[];start=time.perf_counter()
    for inp in sorted((ROOT/'inputs').glob('*.json')):
        model=inp.name.startswith('model_')
        with tempfile.TemporaryDirectory(prefix='fs1_') as td:
            d=Path(td)
            for py in ['exact.py','verify.py']:shutil.copy2(ROOT/py,d/py)
            shutil.copy2(inp,d/'input.json')
            cmd=[sys.executable,'verify.py','input.json']
            if model:cmd+=['--model']
            else:
                shutil.copy2(ROOT/'certificates'/inp.name,d/'certificate.json');cmd+=['certificate.json']
            cmd+=['--output','checked.json']
            t=time.perf_counter();r=subprocess.run(cmd,cwd=d,capture_output=True,text=True,timeout=25)
            if r.returncode:raise RuntimeError(inp.name+': '+r.stdout+r.stderr)
            result=json.loads((d/'checked.json').read_text());seconds=time.perf_counter()-t
            out.append(dict(case=inp.stem,type='operator' if model else 'positive_pair',result=result,seconds=seconds))
    report=dict(status='PASS',pair_certificates=sum(a['type']=='positive_pair' for a in out),
                operator_certificates=sum(a['type']=='operator' for a in out),
                target_pair_evaluations=3*sum(a['type']=='positive_pair' for a in out),
                total_seconds=time.perf_counter()-start,processes=len(out),results=out,
                scope='Exact finite instances; no simulated confidence coverage, estimator, endpoint optimization, or formal proof assistant.')
    dump(ROOT/'results/replay.json',report);print({k:v for k,v in report.items() if k!='results'})
if __name__=='__main__':main()
