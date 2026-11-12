#!/usr/bin/env python3
"""Temporary initial-chunk experiment; production and prior artifacts stay intact."""
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
import subprocess
import time

HERE=Path(__file__).resolve().parent
VARIANTS=['baseline','head1','head2','head4','head8','head32','hybrid16']

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,v): p.write_text(json.dumps(v,indent=2)+'\n')
def invoke(cmd,timeout=900):
    p=subprocess.run(list(map(str,cmd)),capture_output=True,text=True,timeout=timeout)
    if p.returncode or p.stderr: raise RuntimeError(f'{cmd}\nexit={p.returncode}\n{p.stdout}\n{p.stderr}')
    return p.stdout

def suite(root):
    cases=[]
    def add(mode,x,cutoffs,threads,numerator=1,denominator=1):
        for cutoff in cutoffs:
            for t in threads:
                cases.append(dict(case_id=f'{mode}-{cutoff}-{x}-t{t}-L{numerator}over{denominator}',
                    mode=mode,cutoff=cutoff,x=x,threads=t,numerator=numerator,denominator=denominator,
                    blocks=3 if x==10**14 else 5,repeats=4 if x==10**14 else 6))
    add('p2',8*10**10,('zero','meissel','lehmer'),(1,2,4,8),6,5)
    for x in (10**10,10**12):
        add('p2',x,('meissel','lehmer'),(1,2,4,8))
        add('p2',x,('zero',),(2,8))
        add('phases',x,('meissel','lehmer'),(1,2,4,8))
    for mode in ('p2','phases'):
        add(mode,10**14,('meissel','lehmer'),(1,8))
    assert len(cases)==56 and len({c['case_id'] for c in cases})==56
    if (root/'suite.json').exists(): raise RuntimeError('suite already frozen')
    save(root/'suite.json',dict(variants=VARIANTS,cases=cases,seed=2026110705,
        design='first32 in K=1/2/4/8 then32; head32 mapping control; baseline and fixed16 unchanged',
        large_cohort='x=1e14 uses3blocks/4calls, predefined before results; others5blocks/6calls',
        suite_generator_sha256=sha(Path(__file__))))
    print('froze56 configurations',flush=True)

def definitions(variant):
    if variant=='hybrid16': return ['-DP2_HYBRID_MULTIPLIER=16']
    return ['-DP2_SELECTOR_BASELINE']+(['-DP2_INITIAL_CHUNK='+variant[4:]] if variant.startswith('head') else [])

def build(root):
    d=root/'build';source=root/'source'
    if not (d/'P2Initial.cpp').exists(): invoke(['python3',HERE/'generate.py',source,d])
    local=['selector_v3_profile.cpp','selector_v3_cases.h','P2Selector.h','P2Hybrid.cpp']
    for n in local: shutil.copy2(HERE.parent/n,d/n)
    harness=d/'selector_v3_profile.cpp';s=harness.read_text()
    anchor='    const char* BuildKind() {\n#if defined(P2_SELECTOR_BASELINE)'
    assert s.count(anchor)==1
    harness.write_text(s.replace(anchor,'    const char* BuildKind() {\n#if defined(P2_INITIAL_CHUNK)\n        return "initial-chunks";\n#elif defined(P2_SELECTOR_BASELINE)'))
    flags=['clang++','-std=c++11','-O3','-g','-fno-omit-frame-pointer','-pthread','-Wall','-Wextra','-Wpedantic',
           '-I'+str(source/'inst/include'),'-I'+str(d)]
    commands=[];common=[]
    for n in ('Phi','PhiParallel','P3','Eratosthenes'):
        p=d/(n+'.o');cmd=flags+['-c',source/('src/'+n+'.cpp'),'-o',p];invoke(cmd);commands.append(list(map(str,cmd)));common.append(p)
    for variant in VARIANTS:
        native=source/'src/P2.cpp' if variant=='baseline' else d/'P2Hybrid.cpp' if variant=='hybrid16' else d/'P2Initial.cpp'
        cmd=flags+definitions(variant)+[harness,native]+common+['-o',root/('profile-'+variant)]
        invoke(cmd);commands.append(list(map(str,cmd)));print('built '+variant,flush=True)
    inputs=[p for p in source.rglob('*') if p.is_file()]+[d/n for n in local+['P2Initial.cpp','InitialChunks.h','generation.json']]+[root/'suite.json']
    save(root/'build-commands.json',commands)
    save(root/'build-metadata.json',dict(compiler=invoke(['clang++','--version']),platform=platform.platform(),
        commit=(root/'commit.txt').read_text().strip(),inputs={str(p.relative_to(root)):sha(p) for p in inputs},
        binaries={v:sha(root/('profile-'+v)) for v in VARIANTS},bench_sha256=sha(Path(__file__))))

def verify(root):
    m=json.loads((root/'build-metadata.json').read_text())
    assert sha(Path(__file__))==m['bench_sha256'],'runner changed afterbuild'
    for n,h in m['inputs'].items(): assert sha(root/n)==h,n
    for v,h in m['binaries'].items(): assert sha(root/('profile-'+v))==h,v

def command(root,variant,c,repeats,expected,mode=None):
    return [root/('profile-'+variant),mode or c['mode'],c['cutoff'],c['x'],c['threads'],repeats,expected,c['numerator'],c['denominator']]

def references(root):
    verify(root);s=json.loads((root/'suite.json').read_text());expect=[];known={};p2={};refs=[]
    oldterms=HERE.parent/'results-matched-profile/reference-terms.csv.gz'
    with gzip.open(oldterms,'rt') as f: terms=list(csv.DictReader(f))
    oldmatched={}
    for cutoff in ('zero','meissel','lehmer'):
        if cutoff=='zero': a=0
        else:
            power=3 if cutoff=='meissel' else 4
            lo,hi=0,8*10**10
            while lo<hi:
                mid=(lo+hi+1)//2
                if mid**power<=8*10**10:lo=mid
                else:hi=mid-1
            a=sum(int(t['prime'])<=lo for t in terms)
        oldmatched[cutoff]=(a,sum(int(t['term']) for t in terms if int(t['i'])>a))
    cli=shutil.which('primecount');assert cli
    for x in sorted({c['x'] for c in s['cases'] if c['mode']=='phases'}):
        values=[]
        for method in ('--gourdon','--deleglise-rivat'):
            cmd=[cli,x,method,'--threads=1'];v=int(invoke(cmd).strip());values.append(v);refs.append(dict(command=list(map(str,cmd)),result=v))
        assert values[0]==values[1];known[x]=values[0]
    for c in s['cases']:
        if c['mode']=='phases': value=known[c['x']];provenance='primecount-gourdon-and-deleglise-rivat'
        else:
            key=(c['x'],c['cutoff'],c['numerator'],c['denominator'])
            if key not in p2:
                cc=dict(c,threads=1);row=list(csv.DictReader(io.StringIO(invoke(command(root,'baseline',cc,1,-1,'reference')))))[0]
                value=int(row['result']);provenance='production-P2-differential'
                if c['x']==8*10**10:
                    a,independent=oldmatched[c['cutoff']];assert int(row['a'])==a and value==independent
                    provenance='independent-saved-sieve-and-primecount-terms'
                p2[key]=(value,provenance)
            value,provenance=p2[key]
        expect.append(dict(case_id=c['case_id'],expected=value,provenance=provenance))
    save(root/'expectations.json',expect);save(root/'references-metadata.json',dict(independent_counts=refs,unique_p2=len(p2),
        matched_term_sha256=sha(oldterms),cli=cli,cli_sha256=sha(Path(cli).resolve())))
    print('prepared56references',flush=True)

def check_row(row,c,expected,repeat,variant):
    assert int(row['result'])==expected and row['verification']=='expected-match'
    assert (row['method'],int(row['x']),int(row['threads']),int(row['repeat']),int(row['table_numerator']),int(row['table_denominator'])) == (
        c['cutoff'],c['x'],c['threads'],repeat,c['numerator'],c['denominator'])
    assert row['build_kind']==('initial-chunks' if variant.startswith('head') else 'fixed' if variant=='hybrid16' else 'baseline')

def preflight(root):
    verify(root);records=[]
    for cutoff in ('zero','meissel','lehmer'):
        c=dict(mode='p2',cutoff=cutoff,x=10**8,threads=1,numerator=1,denominator=1)
        oracle=list(csv.DictReader(io.StringIO(invoke(command(root,'baseline',c,1,-1,'oracle')))))[0]
        expected=int(oracle['result'])
        for t in (1,2,4,8):
            cc=dict(c,threads=t)
            for v in VARIANTS:
                rows=list(csv.DictReader(io.StringIO(invoke(command(root,v,cc,2,expected)))))
                for i,row in enumerate(rows):check_row(row,cc,expected,i,v)
                records.append(dict(cutoff=cutoff,threads=t,variant=v,expected=expected,calls=len(rows)))
    save(root/'harness-preflight.json',dict(records=records,calls=sum(r['calls'] for r in records)))
    print('harness168independent matches',flush=True)

def timings(root):
    verify(root);s=json.loads((root/'suite.json').read_text());refs={r['case_id']:r for r in json.loads((root/'expectations.json').read_text())}
    output=root/'timings.csv.gz'
    if output.exists():raise RuntimeError('refusing overwrite timings')
    rng=random.Random(s['seed']);runs=[];writer=None;row_count=0;started=time.time();process=0
    with gzip.open(output,'wt',newline='') as f:
        for block in range(max(c['blocks'] for c in s['cases'])):
            cases=[c for c in s['cases'] if block<c['blocks']];rng.shuffle(cases)
            for position,c in enumerate(cases):
                variants=list(s['variants']);rng.shuffle(variants);ref=refs[c['case_id']]
                for order,v in enumerate(variants):
                    cmd=command(root,v,c,c['repeats'],ref['expected']);rows=list(csv.DictReader(io.StringIO(invoke(cmd))))
                    assert len(rows)==c['repeats']
                    for i,row in enumerate(rows):
                        check_row(row,c,ref['expected'],i,v)
                        row.update(case_id=c['case_id'],variant=v,block=block,order=order,process=process,
                                   initial_k=int(v[4:]) if v.startswith('head') else 0,reference_provenance=ref['provenance'])
                        if writer is None:writer=csv.DictWriter(f,fieldnames=list(row));writer.writeheader()
                        writer.writerow(row);row_count+=1
                    runs.append(dict(block=block,case_id=c['case_id'],variant=v,process=process,command=list(map(str,cmd))));process+=1
                f.flush()
                print(f'block{block+1} case{position+1}/{len(cases)} {c["case_id"]} elapsed={time.time()-started:.1f}s',flush=True)
    verify(root)
    save(root/'timings-metadata.json',dict(rows=row_count,processes=process,seed=s['seed'],runs=runs,sha256=sha(output),elapsed_seconds=time.time()-started))
    print(f'PASS {row_count} timed results',flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=('suite','build','references','preflight','timings'));p.add_argument('root',type=Path);a=p.parse_args()
    globals()[a.action](a.root)

if __name__=='__main__':main()
