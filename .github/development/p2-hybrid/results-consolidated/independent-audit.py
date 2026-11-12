#!/usr/bin/env python3
"""Independent raw-data recomputation; does not import the experiment analyzer."""
import collections, csv, gzip, hashlib, json, math, statistics, sys
from pathlib import Path
R=Path(__file__).resolve().parent
P=Path('/Users/owacasa/RcppAlgos')
E=P/'.github/development/p2-hybrid/consolidated'
FIELDS=('setup_ms','phi_ms','p2_ms','p3_ms','teardown_ms','total_ms')
MED=statistics.median

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    with p.open(newline='') as f:return list(csv.DictReader(f))
def near(actual,expected):
    if expected is None:assert actual in ('',None), (actual,expected)
    elif isinstance(expected,bool):assert str(actual)==str(expected),(actual,expected)
    else:assert math.isclose(float(actual),expected,rel_tol=1e-12,abs_tol=1e-10),(actual,expected)
def check_fields(record,expected):
    for k,v in expected.items():near(record[k],v)

suite=json.loads((R/'suite.json').read_text())
build=json.loads((R/'build-metadata.json').read_text())
assert sha(E/'P2SelectorCandidate.h')==build['candidate_sha256']==suite['candidate_sha256']
assert sha(E/'bench.py')==build['runner_sha256']==suite['runner_sha256']
assert sha(E/'plan.json')==suite['plan_sha256']
for p,h in build['inputs'].items():assert sha(R/p)==h,p
for v,h in build['binaries'].items():assert sha(R/('profile-'+v))==h,v
production=0
for p,h in build['inputs'].items():
    if p.startswith(('source/src/','source/inst/','source/R/','source/tests/','source/man/')) or p in ('source/DESCRIPTION','source/NAMESPACE'):
        current=P/p[len('source/'):]
        if current.exists():assert sha(current)==h,str(current);production+=1
metadata_refs=json.loads((R/'references-metadata.json').read_text())
assert sha(R/'expectations.json')==metadata_refs['expectations_sha256']
references={a['case_id']:a for a in json.loads((R/'expectations.json').read_text())}
assert set(references)=={a['case_id'] for a in suite['cases']}
all_p2={};all_prime={};audits=[]
for cohort in sys.argv[1:] or ['focus']:
    cases={a['case_id']:a for a in suite['cases'] if a['cohort']==cohort}
    m=json.loads((R/(cohort+'-metadata.json')).read_text())
    raw=R/(cohort+'.csv.gz')
    assert sha(raw)==m['sha256'];assert sha(R/'suite.json')==m['suite_sha256']
    assert m['candidate_sha256']==suite['candidate_sha256']
    with gzip.open(raw,'rt',newline='') as f:rows=list(csv.DictReader(f))
    wanted={(cid,v,b,i) for cid,c in cases.items() for v in c['variants'] for b in range(c['blocks']) for i in range(c['repeats'])}
    seen=set();groups=collections.defaultdict(list);pids={};orders={};geometry={}
    for a in rows:
        cid,v,b,i=a['case_id'],a['variant'],int(a['block']),int(a['repeat']);c=cases[cid]
        key=(cid,v,b,i);assert key in wanted and key not in seen;seen.add(key)
        assert (a['mode'],a['method'],int(a['x']),int(a['threads']),int(a['table_numerator']),int(a['table_denominator']))==(c['mode'],c['cutoff'],c['x'],c['threads'],c['numerator'],c['denominator'])
        assert a['cohort']==cohort and a['cache_state']==('first' if i==0 else 'warm')
        assert int(a['initial_k'])==(32 if v=='baseline' else 1)
        bk='baseline' if v=='baseline' else 'phi-initial-1' if v=='phi1' else 'fixed' if v.startswith('fixed') else 'candidate' if v=='candidate' else 'selector-'+v
        assert a['build_kind']==bk
        assert a['verification']=='expected-match' and a['measurement']=='timing'
        assert int(a['expected'])==int(a['result'])==references[cid]['expected']
        assert a['reference_provenance']==references[cid]['provenance']
        pid=int(a['process']);identity=(cid,v,b)
        assert pids.setdefault(pid,identity)==identity
        assert orders.setdefault(identity,(pid,int(a['order'])))==(pid,int(a['order']))
        numeric=(int(a['x']),int(a['a']))
        assert all_p2.setdefault(numeric,int(a['p2']))==int(a['p2'])
        if a['prime_count']:
            assert int(a['phi'])+int(a['a'])-1-int(a['p2'])-int(a['p3'])==int(a['prime_count'])
            assert all_prime.setdefault(int(a['x']),int(a['prime_count']))==int(a['prime_count'])
        assert int(a['table_limit'])==math.isqrt(int(a['x']))*c['numerator']//c['denominator']
        assert int(a['phi_queries'])+int(a['stream_queries'])==int(a['non_table_queries'])
        g=tuple(a[k] for k in ('selection_reason','selected_q','selected_k','phi_queries','stream_queries','direct_queries','phi_worker_budget','phi_workers','stream_worker_budget','stream_workers','candidate_q','candidate_stream_first','candidate_stream_queries','candidate_stream_span','candidate_stream_worker_budget'))
        assert geometry.setdefault((cid,v),g)==g
        vals={f:float(a[f]) for f in FIELDS};assert all(math.isfinite(n) and n>=0 for n in vals.values())
        near(vals['total_ms'],sum(vals[f] for f in FIELDS[:-1]))
        groups[(cid,v,b,a['cache_state'])].append(vals)
    assert wanted==seen and len(rows)==m['rows']
    assert len(pids)==m['processes']==len(m['runs'])
    assert len({a['process'] for a in m['runs']})==len(m['runs'])
    assert len({(a['case_id'],a['variant'],a['block']) for a in m['runs']})==len(m['runs'])
    for run in m['runs']:
        cid,v,b=run['case_id'],run['variant'],run['block'];c=cases[cid]
        assert pids[run['process']]==(cid,v,b)
        cmd=run['command'];assert Path(cmd[0]).resolve()==(R/('profile-'+v)).resolve()
        assert cmd[1:]==[c['mode'],c['cutoff'],str(c['x']),str(c['threads']),str(c['repeats']),str(references[cid]['expected']),str(c['numerator']),str(c['denominator'])]
    for cid,c in cases.items():
        for b in range(c['blocks']):assert sorted(orders[cid,v,b][1] for v in c['variants'])==list(range(len(c['variants'])))
    process={k:{f:MED(r[f] for r in rs) for f in FIELDS} for k,rs in groups.items()}
    procrows=read(R/(cohort+'-process-medians.csv'));assert len(procrows)==len(process)
    procseen=set()
    for a in procrows:
        k=(a['case_id'],a['variant'],int(a['block']),a['cache_state']);assert k not in procseen;procseen.add(k)
        assert int(a['samples'])==len(groups[k])==(1 if k[-1]=='first' else cases[k[0]]['repeats']-1)
        check_fields(a,process[k])
    summaries=read(R/(cohort+'-summary.csv'));sumseen=set()
    for a in summaries:
        cid,v,s,f=a['case_id'],a['variant'],a['cache_state'],a['field'];k=(cid,v,s,f);assert k not in sumseen;sumseen.add(k)
        vals=[process[cid,v,b,s][f] for b in range(cases[cid]['blocks'])]
        check_fields(a,{'median_ms':MED(vals),'minimum_ms':min(vals),'maximum_ms':max(vals),'blocks':len(vals)})
    assert len(summaries)==sum(len(c['variants'])*2*len(FIELDS) for c in cases.values())
    comparisons=read(R/(cohort+'-comparisons.csv'));compseen=set();expectedflags=[]
    expectedpairs={}
    for cid,c in cases.items():expectedpairs[cid]=[('candidate',v) for v in c['variants'] if v!='candidate']+[(v,'phi1') for v in ('v1','v2','v3') if v in c['variants']]+[('phi1','baseline')]
    for a in comparisons:
        cid,l,r,s,f=a['case_id'],a['left'],a['right'],a['cache_state'],a['field'];k=(cid,l,r,s,f);assert k not in compseen;compseen.add(k);c=cases[cid]
        assert (l,r) in expectedpairs[cid]
        lv=[process[cid,l,b,s][f] for b in range(c['blocks'])];rv=[process[cid,r,b,s][f] for b in range(c['blocks'])]
        ds=[x-y for x,y in zip(lv,rv)];ps=[100*(x/y-1) for x,y in zip(lv,rv) if y>0]
        expect={'blocks':len(ds),'left_median_ms':MED(lv),'right_median_ms':MED(rv),'delta_median_ms':MED(ds),'delta_min_ms':min(ds),'delta_max_ms':max(ds),'percent_median':MED(ps) if ps else None,'percent_min':min(ps) if ps else None,'percent_max':max(ps) if ps else None,'percent_pairs':len(ps),'improved_blocks':sum(d<0 for d in ds),'regressed_blocks':sum(d>0 for d in ds),'all_improve':all(d<0 for d in ds),'all_regress':all(d>0 for d in ds)}
        check_fields(a,expect)
        if l=='candidate' and r in ('baseline','phi1') and f==('p2_ms' if c['mode']=='p2' else 'total_ms'):
            short=s=='warm' and geometry[cid,'candidate'][0]=='short-full-range';allow=0.02 if short else 0.05
            if all(d>0 for d in ds) and MED(ds)>allow and (short or (ps and MED(ps)>5)):
                expectedflags.append((cid,r,s,f,expect,allow))
    assert len(comparisons)==sum(len(v)*2*len(FIELDS) for v in expectedpairs.values())
    flags=json.loads((R/(cohort+'-flags.json')).read_text());assert len(flags)==len(expectedflags)
    for flag,(cid,r,s,f,vals,allow) in zip(flags,expectedflags):
        assert (flag['case_id'],flag['right'],flag['cache_state'],flag['field'])==(cid,r,s,f);check_fields(flag,vals);near(flag['absolute_allowance_ms'],allow)
    bests=read(R/(cohort+'-best-controls.csv'));assert len(bests)==2*len(cases);bestseen=set()
    for a in bests:
        cid,s=a['case_id'],a['cache_state'];assert (cid,s) not in bestseen;bestseen.add((cid,s));c=cases[cid];f='p2_ms' if c['mode']=='p2' else 'total_ms'
        opts=[v for v in c['variants'] if v in ('baseline','phi1') or v.startswith('fixed')]
        vm={v:MED(process[cid,v,b,s][f] for b in range(c['blocks'])) for v in opts};best=min(opts,key=vm.get);cv=MED(process[cid,'candidate',b,s][f] for b in range(c['blocks']))
        assert a['best_control']==best
        check_fields(a,{'candidate_ms':cv,'best_control_ms':vm[best],'candidate_delta_ms':cv-vm[best],'candidate_percent':100*(cv/vm[best]-1) if vm[best]>0 else None})
    audits.append(dict(cohort=cohort,rows=len(rows),cases=len(cases),processes=len(pids),process_distributions=len(process)*len(FIELDS),summaries=len(summaries),comparisons=len(comparisons),best_controls=len(bests),flags=len(flags),raw_sha256=sha(raw),all_verified=True))
out=dict(auditor='independent read/recompute; analyzer not imported',auditor_sha256=sha(Path(__file__)),cohorts=audits,frozen_inputs_verified=len(build['inputs']),binaries_verified=len(build['binaries']),production_files_matched_to_snapshot=production,candidate_sha256=build['candidate_sha256'],numeric_p2_values_consistent_across_modes_tables_threads_variants=len(all_p2),prime_counts_consistent=len(all_prime),all_verified=True)
name='independent-'+('-'.join(sys.argv[1:] or ['focus']))+'-audit.json'
(R/name).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
