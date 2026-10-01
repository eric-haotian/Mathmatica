#!/usr/bin/env python3
"""Verify the frozen manuscript evidence in a disposable, independent work tree.

Python standard library only. No optimizer, network request, experiment generation,
repository account, or source-truth file is required by the frozen verifier runs.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent
FS1 = Path('inherited/fs1/EIS_PEM_Math_FS1_Full_Spectrum')
FS2 = Path('inherited/fs2/EIS_PEM_Math_FS2_Sequential_Cost')

def read_json(path: Path):
    return json.loads(path.read_text(encoding='utf-8'))

def validate_inventory(code: Path) -> None:
    expected = [
        (FS1/'inputs', '*.json', 26), (FS1/'certificates', '*.json', 22),
        (FS2/'inputs', '*.json', 18), (FS2/'certificates', '*.json', 18),
        (FS2/'sequential/inputs', '*.json', 6),
        (FS2/'sequential/certificates', '*.json', 12),
        (FS2/'results', 'sequential_*.json', 2),
        (Path('new_checks/certificates'), 'guard_*.json', 6),
        (Path('new_checks/certificates'), 'information_*.json', 12),
        (Path('design_audit/inputs'), '*.json', 16),
        (Path('design_audit/certificates'), '*.json', 16),
    ]
    for directory, pattern, count in expected:
        found = len(list((code/directory).glob(pattern)))
        if found != count:
            raise RuntimeError(f'Incomplete inventory: {directory}/{pattern}: {found}, expected {count}')

def validate_hashes() -> int:
    manifest = ROOT/'MANIFEST_SHA256.json'
    if not manifest.exists():
        raise RuntimeError('MANIFEST_SHA256.json is missing')
    entries = read_json(manifest)['files']
    for name, expected in entries.items():
        p = ROOT/name
        if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest() != expected:
            raise RuntimeError(f'Checksum mismatch or missing file: {name}')
    return len(entries)

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, default=ROOT/'replay_output')
    ap.add_argument('--timeout', type=int, default=900, help='Seconds per suite; inherited per-case limits also apply.')
    ap.add_argument('--skip-hashes', action='store_true', help='For package assembly only; normal use verifies the manifest.')
    args = ap.parse_args()
    if args.timeout <= 0: ap.error('--timeout must be positive')
    output = args.output.resolve()
    if output == ROOT or ROOT/'code' in output.parents:
        ap.error('Choose an output directory outside the frozen code tree')
    output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    count = None if args.skip_hashes else validate_hashes()
    validate_inventory(ROOT/'code')
    commands = [
        ('full_spectrum', FS1/'replay.py'),
        ('full_spectrum_mechanisms', FS1/'check_mechanisms.py'),
        ('sequential', FS2/'replay.py'),
        ('sequential_mechanisms', FS2/'check_mechanisms.py'),
        ('feasibility_and_design', Path('new_checks/replay_content.py')),
        ('identification_threshold', Path('design_audit/replay.py')),
        ('relative_width_and_algebra', Path('closure/check_closeout.py')),
        ('tables', Path('regenerate_tables.py')),
        ('design_tables', Path('new_checks/make_tables.py')),
        ('full_spectrum_table', FS1/'report_results.py'),
    ]
    runs = []
    with tempfile.TemporaryDirectory(prefix='sinum_submission_replay_') as tmp:
        code = Path(tmp)/'code'
        shutil.copytree(ROOT/'code', code, ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        (code/'manuscript').mkdir(exist_ok=True)
        (code/FS1/'documents').mkdir(exist_ok=True)
        for name, rel in commands:
            script = code/rel
            boot = 'import runpy,sys; from pathlib import Path; p=Path(sys.argv[1]); sys.path.insert(0,str(p.parent)); sys.argv=[str(p)]; runpy.run_path(str(p),run_name="__main__")'
            t = time.perf_counter()
            result = subprocess.run([sys.executable, '-I', '-c', boot, str(script)],
                cwd=script.parent, text=True, encoding='utf-8', capture_output=True,
                timeout=args.timeout, env={**os.environ, 'PYTHONHASHSEED':'0'})
            (output/(name+'.stdout.txt')).write_text(result.stdout, encoding='utf-8')
            (output/(name+'.stderr.txt')).write_text(result.stderr, encoding='utf-8')
            row = {'command':str(rel), 'returncode':result.returncode,
                   'seconds':time.perf_counter()-t}
            runs.append(row)
            if result.returncode:
                (output/'FAILED.json').write_text(json.dumps({'failed':name,'runs':runs},indent=2))
                raise RuntimeError(f'{name} failed; see {output}')
            print(f'PASS {name}', flush=True)
        reports = {
            'full_spectrum.json': code/FS1/'results/replay.json',
            'full_spectrum_mechanisms.json': code/FS1/'results/mechanism_checks.json',
            'sequential.json': code/FS2/'results/replay.json',
            'sequential_mechanisms.json': code/FS2/'results/mechanism_checks.json',
            'feasibility_and_design.json': code/'new_checks/results/replay.json',
            'identification_threshold.json': code/'design_audit/results/replay.json',
            'relative_width_and_algebra.json': code/'closure/results/closeout_checks.json',
        }
        for name, path in reports.items(): shutil.copy2(path,output/name)
        tab = output/'tables'; tab.mkdir(exist_ok=True)
        for path in (code/'manuscript').glob('*.tex'): shutil.copy2(path,tab/path.name)
        for path in (code/FS1/'results').glob('summary.*'): shutil.copy2(path,tab/('full_spectrum_'+path.name))
        full = read_json(reports['full_spectrum.json']); seq=read_json(reports['sequential.json'])
        feas=read_json(reports['feasibility_and_design.json']); ident=read_json(reports['identification_threshold.json'])
        close=read_json(reports['relative_width_and_algebra.json'])
        tasks=full['processes']+len(seq['runs'])+feas['fresh_process_tasks']+ident['case_count']
        fs1m=read_json(reports['full_spectrum_mechanisms.json']);fs2m=read_json(reports['sequential_mechanisms.json'])
        checks=fs1m['count']+fs2m['count']+len(feas['mechanisms'])+len(feas['invalid_certificate_tests'])+ident['mutation_count']+close['positive_checks']+close['rejection_checks']
        rejections=fs1m['rejections']+fs2m['rejection_checks']+len(feas['invalid_certificate_tests'])+ident['mutation_count']+close['rejection_checks']
        if tasks != 80 or checks != 154 or rejections != 53:
            raise RuntimeError(f'Unexpected accounting: {tasks}, {checks}, {rejections}')
    report={'status':'PASS','python':sys.version,'manifest_files_checked':count,
        'independent_certificate_tasks':tasks, 'checks_including_rejections':checks,
        'rejection_checks':rejections,'continuous_endpoint_certificates_inside_two_trace_tasks':12,
        'new_experiments':0,'new_optimization_runs':0,'commands':runs,
        'seconds':time.perf_counter()-started,
        'scope':'Replay of supplied exact inputs and certificates, including complete table regeneration; internal computational verification, not formal proof verification or a coverage simulation.'}
    (output/'SUMMARY.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(f'PASS: {tasks} certificate tasks; {checks} checks including {rejections} required rejections. {output}',flush=True)

if __name__=='__main__':
    try: main()
    except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as exc:
        print(f'FAILED: {exc}', file=sys.stderr); sys.exit(1)
