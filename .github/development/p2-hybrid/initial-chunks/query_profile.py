#!/usr/bin/env python3
"""Generate and run first-query diagnostics on original P2 scheduling only.

No production file is modified. Timing fields are diagnostic elapsed intervals,
not substitutes for the separate uninstrumented variant comparisons.
"""
import argparse
from array import array
import collections
import ctypes
import gzip
import hashlib
import json
import math
from pathlib import Path
import random
import shutil
import statistics
import subprocess

HERE=Path(__file__).resolve().parent

def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def save(path,value): path.write_text(json.dumps(value,indent=2)+'\n')
def invoke(command):
    p=subprocess.run(list(map(str,command)),text=True,capture_output=True,timeout=900)
    if p.returncode or p.stderr:
        raise RuntimeError(f'{command}\nexit={p.returncode}\n{p.stdout}\n{p.stderr}')
    return p.stdout

def once(text,old,new):
    if text.count(old)!=1: raise RuntimeError('source drift at '+old[:100])
    return text.replace(old,new)

def instrument(text):
    text='#include "query_trace.h"\n'+text
    text=once(text,'            const std::int64_t piQuotient = Phi(quotient, c, tables) + c - 1;',
        '            const double queryBegin = QueryProfile::queryBegin(i);\n'
        '            const std::int64_t piQuotient = Phi(quotient, c, tables) + c - 1;\n'
        '            QueryProfile::queryEnd(i, primes[i], quotient, c, piQuotient, queryBegin);')
    text=once(text,'        if (last == first) return sum;',
        '        QueryProfile::range(first, last);\n        if (last == first) return sum;')
    text=once(text,'        if (workers == 1) return sum + countP2Range(x, first, last, tables);',
        '        QueryProfile::current->workerCount = workers;\n'
        '        if (workers == 1) {\n'
        '            QueryProfile::WorkerScope worker(0);\n'
        '            QueryProfile::ChunkScope chunk(first, last);\n'
        '            return sum + countP2Range(x, first, last, tables);\n        }')
    text=once(text,'            batch.launch([&, slot] {',
        '            batch.launch([&, slot] {\n                QueryProfile::WorkerScope worker(slot);')
    text=once(text,'                    subtotal += countP2Range(x, begin,',
        '                    QueryProfile::ChunkScope chunk(begin, std::min(begin + CHUNK_SIZE, last));\n'
        '                    subtotal += countP2Range(x, begin,')
    return text

def build(root):
    directory=root/'query-build'
    if directory.exists(): raise RuntimeError('refusing existing query-build')
    source=root/'source'
    generated=instrument((source/'src/P2.cpp').read_text())
    directory.mkdir()
    (directory/'P2Query.cpp').write_text(generated)
    for name in ('query_trace.h','query_profile.cpp'):
        shutil.copy2(HERE/name,directory/name)
    flags=['clang++','-std=c++11','-O3','-g','-fno-omit-frame-pointer','-pthread',
           '-Wall','-Wextra','-Wpedantic','-I'+str(source/'inst/include'),'-I'+str(directory)]
    commands=[]
    for kind in ('normal','asan','tsan'):
        extra=[] if kind=='normal' else ['-O1','-fsanitize='+('address,undefined' if kind=='asan' else 'thread')]
        command=flags+extra+[directory/'query_profile.cpp',directory/'P2Query.cpp',source/'src/Phi.cpp',
            source/'src/Eratosthenes.cpp','-o',root/('query-profile'+('' if kind=='normal' else '-'+kind))]
        output=invoke(command)
        (root/('query-'+kind+'-build.log')).write_text(output)
        commands.append(list(map(str,command)))
    save(root/'query-build-commands.json',commands)
    paths=[p for p in source.rglob('*') if p.is_file()]+[p for p in directory.iterdir() if p.is_file()]
    save(root/'query-build-metadata.json',dict(
        source_sha256={str(p.relative_to(root)):digest(p) for p in paths},
        binary_sha256={name:digest(root/name) for name in ('query-profile','query-profile-asan','query-profile-tsan')},
        generator_sha256=digest(Path(__file__)),compiler=invoke(['clang++','--version']),
        definition='Original P2 scheduling; clocks only first32 Phi queries, first64 chunks, and worker lifetimes.'))
    print('Built original-scheduler query diagnostic: normal, ASan/UBSan, TSan.',flush=True)

def verify(root):
    m=json.loads((root/'query-build-metadata.json').read_text())
    for name,h in m['source_sha256'].items():
        if digest(root/name)!=h: raise RuntimeError('query build source changed: '+name)
    for name,h in m['binary_sha256'].items():
        if digest(root/name)!=h: raise RuntimeError('query executable changed: '+name)

def check(row):
    eps=0.00001
    assert row['result']==row['expected']
    first,last,workers=row['first'],row['last'],row['workers']
    assert len(row['worker_records'])==workers
    assert sum(w['queries'] for w in row['worker_records'])==last-first
    assert len(row['queries'])==min(32,last-first)
    assert [q['index'] for q in row['queries']]==list(range(first,min(first+32,last)))
    work={w['slot']:w for w in row['worker_records']}
    assert len(work)==workers
    chunks=row['chunks']
    assert len(chunks)==(min(64,(last-first+31)//32) if workers>1 else int(last>first))
    for index,c in enumerate(chunks):
        assert c['first']==first+(32*index if workers>1 else 0)
        assert c['last']==(min(c['first']+32,last) if workers>1 else last)
        w=work[c['worker']]
        assert -eps<=w['begin_us']<=c['begin_us']+eps<=c['end_us']+2*eps<=w['end_us']+3*eps
        assert w['end_us']<=row['elapsed_us']+eps
    for q in row['queries']:
        assert q['quotient']==row['x']//q['prime']
        assert q['term']==q['pi']-(q['index']-1)
        c=next(c for c in chunks if c['first']<=q['index']<c['last'])
        assert c['worker']==q['worker']
        assert c['begin_us']-eps<=q['begin_us']<=q['end_us']+eps<=c['end_us']+2*eps
    assert sum(w['chunks'] for w in row['worker_records'])==((last-first+31)//32 if workers>1 else int(last>first))

def command(root,case,repeats,kind='normal'):
    return [root/('query-profile'+('' if kind=='normal' else '-'+kind)),case['cutoff'],case['x'],
            case['threads'],repeats,case['expected'],case.get('table_numerator',case.get('numerator')),
            case.get('table_denominator',case.get('denominator'))]

def prepare(root):
    output=root/'query-suite.json'
    if output.exists(): raise RuntimeError('refusing to overwrite query suite')
    suite=json.loads((root/'suite.json').read_text())
    expected={r['case_id']:r for r in json.loads((root/'expectations.json').read_text())}
    cases=[]
    for case in suite['cases']:
        if case['mode']!='p2':continue
        selected=(case['x']==80000000000 or
            (case['x']==1000000000000 and case['cutoff']!='zero' and case['threads'] in (1,8)) or
            (case['x']==100000000000000 and case['threads']==8))
        if not selected:continue
        value=dict(case);value['expected']=expected[case['case_id']]['expected'];cases.append(value)
    if len(cases)!=18:raise RuntimeError('unexpected diagnostic subset size: '+str(len(cases)))
    save(output,dict(cases=cases,source_suite_sha256=digest(root/'suite.json'),
        expectations_sha256=digest(root/'expectations.json'),blocks=3,repeats=4))
    print('Prepared eighteen diagnostic configurations.',flush=True)

def small_references(x):
    flags=bytearray(b'\x01')*(x+1);flags[0:2]=b'\x00\x00'
    for p in range(2,math.isqrt(x)+1):
        if flags[p]:flags[p*p::p]=b'\x00'*((x-p*p)//p+1)
    pi=array('I',[0])*(x+1);count=0
    for i,f in enumerate(flags):count+=f;pi[i]=count
    root=math.isqrt(x);primes=[0]+[p for p in range(2,root+1) if flags[p]]
    cbrt=0
    while (cbrt+1)**3<=x:cbrt+=1
    cutoffs=dict(zero=0,meissel=pi[cbrt],lehmer=pi[math.isqrt(root)])
    return {c:sum(pi[x//primes[i]]-(i-1) for i in range(a+1,len(primes))) for c,a in cutoffs.items()}

def library():
    path=Path('/opt/homebrew/lib/libprimecount.dylib').resolve(strict=True)
    api=ctypes.CDLL(str(path));api.primecount_pi.argtypes=[ctypes.c_int64];api.primecount_pi.restype=ctypes.c_int64
    api.primecount_set_num_threads.argtypes=[ctypes.c_int];api.primecount_set_num_threads.restype=None
    api.primecount_set_num_threads(1)
    return api,path

def validate(root):
    verify(root);cases=[]
    for cutoff,expected in small_references(100000).items():
        for threads in (1,2,8):cases.append(dict(cutoff=cutoff,x=100000,threads=threads,expected=expected,table_numerator=1,table_denominator=1))
    for cutoff,expected in (('zero',11238880589),('meissel',2069687562)):
        for threads in (1,2,8):cases.append(dict(cutoff=cutoff,x=80000000000,threads=threads,expected=expected,table_numerator=6,table_denominator=5))
    api,path=library();known={};rows=[];commands=[]
    for kind in ('normal','asan','tsan'):
        for case in cases:
            cmd=command(root,case,2,kind);commands.append(list(map(str,cmd)))
            records=[json.loads(s) for s in invoke(cmd).splitlines()];assert len(records)==2
            for row in records:
                check(row)
                for q in row['queries']:
                    n=q['quotient']
                    if n not in known:known[n]=api.primecount_pi(n)
                    assert q['pi']==known[n]
                row['build']=kind;rows.append(row)
    with gzip.open(root/'query-validation.jsonl.gz','wt') as f:
        for row in rows:f.write(json.dumps(row,separators=(',',':'))+'\n')
    save(root/'query-validation.json',dict(rows=len(rows),queries=sum(len(r['queries']) for r in rows),
        independent_pi_queries=len(known),library=str(path),library_sha256=digest(path),
        commands=commands,per_build=len(cases)*2,all_trace_and_reference_checks_pass=True,
        raw_sha256=digest(root/'query-validation.jsonl.gz')))
    print('Query diagnostics passed',len(rows),'rows across normal/ASan/TSan.',flush=True)

def collect(root,suite,blocks,repeats,seed):
    verify(root);configs=json.loads(suite.read_text());configs=configs['cases'] if isinstance(configs,dict) else configs
    output=root/'query-profile.jsonl.gz'
    if output.exists(): raise RuntimeError('refusing to overwrite query matrix')
    rng=random.Random(seed);runs=[];row_count=0
    with gzip.open(output,'wt') as f:
        for block in range(blocks):
            cases=list(configs);rng.shuffle(cases)
            for case in cases:
                cmd=command(root,case,repeats);parsed=[json.loads(s) for s in invoke(cmd).splitlines()]
                assert len(parsed)==repeats
                for repeat,row in enumerate(parsed):
                    assert row['repeat']==repeat;check(row)
                    row.update(block=block,case=case.get('case_id',case.get('id','')))
                    f.write(json.dumps(row,separators=(',',':'))+'\n');row_count+=1
                runs.append(dict(block=block,command=list(map(str,cmd))))
    save(root/'query-profile-metadata.json',dict(rows=row_count,blocks=blocks,repeats=repeats,seed=seed,
        suite=str(suite),suite_sha256=digest(suite),runs=runs,raw_sha256=digest(output)))
    print('Collected',row_count,'original-scheduler query traces.',flush=True)

def audit(root):
    verify(root);path=root/'query-profile.jsonl.gz';meta=json.loads((root/'query-profile-metadata.json').read_text())
    assert digest(path)==meta['raw_sha256']
    with gzip.open(path,'rt') as f:rows=[json.loads(line) for line in f]
    assert len(rows)==meta['rows'];api,lib=library();known={}
    grouped=collections.defaultdict(list)
    for row in rows:
        check(row)
        for q in row['queries']:
            n=q['quotient']
            if n not in known:known[n]=api.primecount_pi(n)
            assert known[n]==q['pi']
        durations=[q['end_us']-q['begin_us'] for q in row['queries']]
        total=sum(durations);first=durations[0] if durations else 0
        metrics=dict(first_query_us=first,first32_sum_us=total,
            first_share=first/total if total else 0,first4_share=sum(durations[:4])/total if total else 0,
            first8_share=sum(durations[:8])/total if total else 0,
            first32_max_us=max(durations) if durations else 0,elapsed_us=row['elapsed_us'],
            first32_query_sum_share=total/row['elapsed_us'] if row['elapsed_us'] else 0)
        for i,value in enumerate(durations):metrics['query'+str(i+1)+'_us']=value
        if row['chunks']:
            c=row['chunks'][0];metrics['initial_chunk_us']=c['end_us']-c['begin_us']
            metrics['first_query_chunk_share']=first/metrics['initial_chunk_us'] if metrics['initial_chunk_us'] else 0
            metrics['initial_chunk_worker_finishes_last']=int(row['worker_records'][c['worker']]['end_us']==max(w['end_us'] for w in row['worker_records']))
        key=(row['x'],row['cutoff'],row['threads'],row['table_numerator'],row['table_denominator'],
             'first' if row['repeat']==0 else 'warm',row['block'])
        grouped[key].append(metrics)
    process={k:{name:statistics.median(r[name] for r in values) for name in values[0]} for k,values in grouped.items()}
    aggregates=collections.defaultdict(lambda:collections.defaultdict(list))
    for key,metrics in sorted(process.items()):
        for name,value in metrics.items():aggregates[key[:-1]][name].append(value)
    result=[]
    for key,metrics in sorted(aggregates.items()):
        rec=dict(zip(('x','cutoff','threads','table_numerator','table_denominator','state'),key))
        rec['metrics']={n:dict(median=statistics.median(v),min=min(v),max=max(v),values=v) for n,v in metrics.items()}
        result.append(rec)
    save(root/'query-analysis.json',dict(rows=len(rows),independent_pi_queries=len(known),
        all_pi_and_trace_checks_pass=True,primecount_library=str(lib),primecount_library_sha256=digest(lib),
        raw_sha256=digest(path),definitions=dict(times='diagnostic microseconds; query windows cover Phi+c-1 only',
            aggregation='median warmed queries within process, then median/min/max across blocks',
            caution='Counterfactual scheduling speedups cannot be inferred by summing or dividing these costs; parallel cache histories differ.'),aggregates=result))
    for r in result:
        m=r['metrics'];print(r['x'],r['cutoff'],r['threads'],r['state'],
            'first',round(m['first_query_us']['median'],3),'first32',round(m['first32_sum_us']['median'],3),
            'share',round(m['first_share']['median']*100,2))

def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=('build','validate','prepare','collect','audit'));p.add_argument('root',type=Path)
    p.add_argument('--suite',type=Path);p.add_argument('--blocks',type=int,default=3);p.add_argument('--repeats',type=int,default=4)
    p.add_argument('--seed',type=int,default=2026110711);a=p.parse_args()
    if a.action=='build':build(a.root)
    elif a.action=='validate':validate(a.root)
    elif a.action=='prepare':prepare(a.root)
    elif a.action=='collect':collect(a.root,a.suite or a.root/'query-suite.json',a.blocks,a.repeats,a.seed)
    else:audit(a.root)

if __name__=='__main__':main()
