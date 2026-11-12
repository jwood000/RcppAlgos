#!/usr/bin/env python3
"""Build and run the frozen matched-cutoff diagnostics outside the package."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import platform
import random
import shutil
import statistics
import subprocess

HERE=Path(__file__).resolve().parent

def save(path,value): path.write_text(json.dumps(value,indent=2)+'\n')
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def invoke(cmd):
    p=subprocess.run(list(map(str,cmd)),capture_output=True,text=True,timeout=300)
    if p.returncode or p.stderr: raise RuntimeError(f'{cmd}\n{p.returncode}\n{p.stdout}\n{p.stderr}')
    return p.stdout

def build(root):
    target=root/'build'
    invoke(['python3',HERE/'instrument_p2.py',root/'source',HERE.parent/'P2Hybrid.cpp',target])
    invoke(['python3',HERE/'instrument_phi.py',target/'phi-counts','--source-root',root/'source'])
    for name in ('trace.h','driver.cpp','phi_metrics.h'):
        shutil.copy2(HERE/name,target/name)
    flags=['clang++','-std=c++11','-O3','-g','-fno-omit-frame-pointer','-pthread','-Wall','-Wextra','-Wpedantic',
           '-I'+str(root/'source/inst/include'),'-I'+str(target),'-DP2_HYBRID_MULTIPLIER=16']
    commands=[]
    for kind in ('plain','timeline','counts','asan','tsan'):
        defs=[] if kind=='plain' else ['-DMATCHED_TRACE']
        if kind in ('counts','asan','tsan'): defs+=['-DMATCHED_COUNTERS']
        extra=['-O1','-fsanitize=address,undefined'] if kind=='asan' else ['-O1','-fsanitize=thread'] if kind=='tsan' else []
        phi=target/'phi-counts/Phi.cpp' if kind in ('counts','asan','tsan') else root/'source/src/Phi.cpp'
        cmd=flags+defs+extra+[target/'driver.cpp',target/'P2Production.cpp',phi,root/'source/src/Eratosthenes.cpp','-o',root/('profile-'+kind)]
        invoke(cmd); commands.append(list(map(str,cmd))); print('built '+kind,flush=True)
    save(root/'build-commands.json',commands)
    save(root/'build-metadata.json',dict(compiler=invoke(['clang++','--version']),platform=platform.platform(),
        inputs={str(p.relative_to(root)):sha(p) for p in list((root/'source').rglob('*'))+list(target.rglob('*')) if p.is_file()},
        binaries={kind:sha(root/('profile-'+kind)) for kind in ('plain','timeline','counts','asan','tsan')}))

def expected(mode,cutoff):
    if mode=='prefix': return {'zero':9887491223,'meissel':718298196}[cutoff]
    if mode.startswith('tail-'): return 1336311383
    return {'zero':11238880589,'meissel':2069687562}[cutoff]

def cases():
    result=[]
    for threads in (1,2):
        for cutoff in ('zero','meissel'):
            for mode in ('baseline','hybrid','split-phi','prefix'):
                result.append(dict(mode=mode,cutoff=cutoff,threads=threads,history='none'))
        for mode in ('tail-phi','tail-stream'):
            for history in ('none','zero-prefix','meissel-prefix','tail'):
                result.append(dict(mode=mode,cutoff='meissel',threads=threads,history=history))
    return result

def verify(root):
    m=json.loads((root/'build-metadata.json').read_text())
    for name,h in m['inputs'].items():
        if sha(root/name)!=h: raise RuntimeError('build input changed '+name)
    for kind,h in m['binaries'].items():
        if sha(root/('profile-'+kind))!=h: raise RuntimeError('executable changed '+kind)

def command(root,kind,case,repeats):
    return [root/('profile-'+kind),case['mode'],case['cutoff'],case['threads'],case['history'],repeats,expected(case['mode'],case['cutoff'])]

def collect(root,kind,blocks,repeats,seed,configs):
    verify(root)
    dest=root/(kind+'.jsonl.gz')
    if dest.exists(): raise RuntimeError('refusing overwrite '+str(dest))
    rng=random.Random(seed); runs=[]; rows=0
    with gzip.open(dest,'wt') as f:
        for block in range(blocks):
            configs=list(configs); rng.shuffle(configs)
            for case in configs:
                cmd=command(root,kind,case,repeats)
                parsed=[json.loads(line) for line in invoke(cmd).splitlines()]
                if len(parsed)!=repeats: raise RuntimeError('row count')
                for i,row in enumerate(parsed):
                    if row['repeat']!=i or row['result']!=expected(case['mode'],case['cutoff']): raise RuntimeError('result')
                    row.update(block=block,build=kind); f.write(json.dumps(row,separators=(',',':'))+'\n'); rows+=1
                runs.append(dict(block=block,command=list(map(str,cmd))))
    save(root/(kind+'-metadata.json'),dict(blocks=blocks,repeats=repeats,seed=seed,runs=runs,rows=rows,sha256=sha(dest)))
    print(kind,rows,'matched rows',flush=True)

def controls(root):
    old=Path((root/'prior-artifact.txt').read_text().strip())
    variants=('baseline','hybrid16','selector-v2','selector-v3')
    configs=[(v,c) for v in variants for c in ('zero','meissel')]
    dest=root/'controls.jsonl.gz'
    if dest.exists(): raise RuntimeError('refusing overwrite')
    rng=random.Random(2026110801); runs=[]; rows=0
    metadata=json.loads((old/'build-metadata.json').read_text())
    with gzip.open(dest,'wt') as f:
        for block in range(9):
            rng.shuffle(configs)
            for variant,cutoff in configs:
                exe=old/('profile-'+variant)
                if sha(exe)!=metadata['binary_sha256'][variant]: raise RuntimeError('old executable changed')
                cmd=[exe,'p2',cutoff,80000000000,2,8,expected('baseline',cutoff),6,5]
                parsed=list(csv.DictReader(io.StringIO(invoke(cmd))))
                if len(parsed)!=8: raise RuntimeError('old row count')
                for row in parsed:
                    if int(row['result'])!=expected('baseline',cutoff): raise RuntimeError('old result mismatch')
                    row.update(block=block,variant=variant); f.write(json.dumps(row,separators=(',',':'))+'\n');rows+=1
                runs.append(dict(block=block,command=list(map(str,cmd))))
    save(root/'controls-metadata.json',dict(rows=rows,runs=runs,sha256=sha(dest)))
    print('controls',rows,'matched rows',flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=('build','validate','counts','timeline','timings','controls'));p.add_argument('root',type=Path);a=p.parse_args()
    if a.action=='build': build(a.root)
    elif a.action=='validate':
        for kind in ('plain','timeline','counts','asan','tsan'):
            collect(a.root,kind,1,3,2026110701,cases())
            (a.root/(kind+'.jsonl.gz')).rename(a.root/('validation-'+kind+'.jsonl.gz'))
            (a.root/(kind+'-metadata.json')).rename(a.root/('validation-'+kind+'-metadata.json'))
    elif a.action=='counts': collect(a.root,'counts',2,3,2026110702,cases())
    elif a.action=='timeline': collect(a.root,'timeline',5,6,2026110703,[c for c in cases() if c['threads']==2 and c['mode'] in ('baseline','hybrid','split-phi')])
    elif a.action=='timings': collect(a.root,'plain',9,8,2026110704,cases())
    else: controls(a.root)

if __name__=='__main__': main()
