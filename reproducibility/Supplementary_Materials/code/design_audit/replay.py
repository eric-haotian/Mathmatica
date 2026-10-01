"""Cold-process replay of the new certificates, then semantic mutation tests."""
from pathlib import Path
import sys,json,subprocess,tempfile,shutil,copy,time
from exact import dump,digest
from verify import verify
ROOT=Path(__file__).resolve().parent
start=time.perf_counter();results=[]
for ip in sorted((ROOT/'inputs').glob('*.json')):
    cp=ROOT/'certificates'/ip.name
    with tempfile.TemporaryDirectory() as td:
        p=Path(td)
        for f in ['exact.py','verify.py']:shutil.copy2(ROOT/f,p/f)
        shutil.copy2(ip,p/'input.json');shutil.copy2(cp,p/'certificate.json')
        run=subprocess.run([sys.executable,'-I','-c',"import sys,runpy;sys.path.insert(0,'.');sys.argv=['verify.py','input.json','certificate.json'];runpy.run_path('verify.py',run_name='__main__')"],cwd=p,capture_output=True,text=True,timeout=30)
        if run.returncode:raise RuntimeError(ip.name+'\n'+run.stderr)
        row=json.loads(run.stdout);row['case']=ip.stem;results.append(row)
base=json.loads((ROOT/'inputs/unknown_K2_uniform.json').read_text())
basec=json.loads((ROOT/'certificates/unknown_K2_uniform.json').read_text())
mutations=[]
def run_bad(name,fn):
    d,c=copy.deepcopy(base),copy.deepcopy(basec);fn(d,c);c['input_sha256']=digest(d)
    try:verify(d,c)
    except (ValueError,KeyError,TypeError,ZeroDivisionError) as e:
        mutations.append({'name':name,'status':'REJECTED','reason':str(e)})
    else:raise RuntimeError('Invalid certificate accepted: '+name)
run_bad('negative_mass',lambda d,c:d['alternative'][0].__setitem__(1,'-1'))
run_bad('negative_calibration',lambda d,c:d['alternative_calibration'].__setitem__(0,'-1'))
run_bad('changed_raw_data',lambda d,c:d['raw_data'].__setitem__(0,'0'))
run_bad('known_calibration_relabel',lambda d,c:d.__setitem__('mode','known'))
run_bad('inflated_gap',lambda d,c:c.__setitem__('target_gap','100'))
run_bad('scaled_wrong_annihilator',lambda d,c:c['annihilator'].__setitem__(0,'100'))
run_bad('duplicate_frequency',lambda d,c:d['frequencies'].__setitem__(1,d['frequencies'][0]))
run_bad('false_additional_row',lambda d,c:c.__setitem__('additional_row_difference',['0','0']))
run_bad('zero_dual',lambda d,c:c.__setitem__('extended_design_annihilator_dual',['0']*len(c['extended_design_annihilator_dual'])))
run_bad('introduced_noise',lambda d,c:d.__setitem__('error_radius','1/1000'))
report={'case_count':len(results),'certificate_results':results,'mutation_count':len(mutations),'mutation_results':mutations,
        'elapsed_seconds':time.perf_counter()-start,'scope':'new exact finite examples only; no optimization, simulation, all-parameter proof verification or external peer review'}
dump(ROOT/'results/replay.json',report)
print('PASS',len(results),'isolated certificates;',len(mutations),'semantic rejections')
