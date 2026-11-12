from pathlib import Path
import csv,hashlib,json,random,statistics,subprocess,sys,time
root=Path(sys.argv[1]).resolve()
out=root/'worker-cap-timing'
out.mkdir()
source=root/'native-validation/source'
units=('Phi','PhiParallel','P3','PiMeissel','PiLehmer','Eratosthenes')
flags=['clang++','-std=c++11','-O3','-pthread','-Wall','-Wextra','-Wpedantic','-I'+str(source/'inst/include')]
commands=[]
def invoke(cmd,log):
    result=subprocess.run(list(map(str,cmd)),capture_output=True,text=True,timeout=300)
    log.write_text(result.stdout+result.stderr)
    commands.append(list(map(str,cmd)))
    if result.returncode or result.stderr: raise RuntimeError(str(cmd)+'\n'+result.stdout+result.stderr)
    return result.stdout
objects=[]
for name in units:
    obj=out/(name+'.o')
    invoke(flags+['-c',source/'src'/(name+'.cpp'),'-o',obj],out/(name+'-build.log'))
    objects.append(obj)
harness=r"""#include "VariantP2.cpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
int main(int argc, char** argv) {
    if (argc != 2) return 2;
    const std::size_t budget = std::stoul(argv[1]);
    const PrimeCounting::PhiTables tables(1000000);
    for (int repeat = 0; repeat < 4; ++repeat) {
        const auto start = std::chrono::steady_clock::now();
        const auto result = PrimeCounting::countPhiPrefix(1000000000000LL, 3305, 3369, budget, tables);
        const auto end = std::chrono::steady_clock::now();
        const double ms = std::chrono::duration<double, std::milli>(end-start).count();
        std::cout << repeat << ',' << result << ',' << std::setprecision(12) << ms << '\n';
    }
}
"""
inputs={'before':root/'before/staged-P2.cpp','after':source/'src/P2.cpp'}
for variant,path in inputs.items():
    directory=out/variant
    directory.mkdir()
    (directory/'VariantP2.cpp').write_bytes(path.read_bytes())
    (directory/'timing.cpp').write_text(harness)
    invoke(flags+[directory/'timing.cpp']+objects+['-o',directory/'timing'],out/(variant+'-build.log'))
print('Timing binaries built; starting 20 fresh processes.',flush=True)
rng=random.Random(20261111)
rows=[]
expected=None
for block in range(5):
    budgets=[8,36];rng.shuffle(budgets)
    for budget in budgets:
        variants=['before','after'];rng.shuffle(variants)
        for order,variant in enumerate(variants):
            data=invoke([out/variant/'timing',budget],out/(str(block)+'-'+str(budget)+'-'+variant+'.log'))
            for row in csv.reader(data.splitlines()):
                repeat,result,ms=int(row[0]),int(row[1]),float(row[2])
                if expected is None: expected=result
                if result!=expected: raise RuntimeError('count mismatch')
                rows.append(dict(block=block,budget=budget,variant=variant,order=order,repeat=repeat,result=result,ms=ms))
with (out/'timings.csv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
summaries=[]
for budget in (8,36):
    for state in ('first','warm'):
        values={}
        for variant in ('before','after'):
            values[variant]=[statistics.median(r['ms'] for r in rows if r['budget']==budget and r['variant']==variant and r['block']==block and ((r['repeat']==0)==(state=='first'))) for block in range(5)]
        deltas=[a-b for a,b in zip(values['after'],values['before'])]
        ratios=[100*(a/b-1) for a,b in zip(values['after'],values['before'])]
        summaries.append(dict(budget=budget,state=state,before_median_ms=statistics.median(values['before']),after_median_ms=statistics.median(values['after']),paired_delta_ms=statistics.median(deltas),paired_percent=statistics.median(ratios),all_slower=all(d>0 for d in deltas),all_faster=all(d<0 for d in deltas)))
result=dict(rows=len(rows),processes=20,expected=expected,seed=20261111,geometry=dict(x=1000000000000,L=1000000,first=3305,last=3369,terms=64,jobs=33),note='Five randomized matched blocks; fresh process first call and median of three warm calls. Time only Phi prefix, excluding table setup. Before is maintainer staged P2, after is cleanup. Budget8 is an unchanged-launch control; budget36 changes36 launches to33. Timings on this host are evidence for these cases only, not selector tuning or general scaling. The normal R thread ceiling is10; this private-helper diagnostic deliberately exercises a larger budget.',summary=summaries,commands=commands,source_sha256={v:hashlib.sha256(p.read_bytes()).hexdigest() for v,p in inputs.items()})
(out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(summaries,indent=2),flush=True)
print('PASS:',len(rows),'matching counts',flush=True)
