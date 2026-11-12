#!/usr/bin/env python3
"""Sanitizer validation for temporary initial-chunk kernels; never time here."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

HERE=Path(__file__).resolve().parent
SANITIZERS={'asan-ubsan':'address,undefined','tsan':'thread'}

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def invoke(command,log):
    result=subprocess.run(list(map(str,command)),capture_output=True,text=True,timeout=300)
    log.write_text(result.stdout+result.stderr)
    if result.returncode or result.stderr:
        raise RuntimeError(f'failed {command}\n{result.stdout}\n{result.stderr}')
    return result.stdout

def main():
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args()
    root=a.root.resolve();source=root/'source';generated=root/'build';out=root/'initial-validation'
    if out.exists():raise RuntimeError('refusing existing validation directory')
    out.mkdir()
    units=('Phi.cpp','PhiParallel.cpp','P3.cpp','PiMeissel.cpp','PiLehmer.cpp','Eratosthenes.cpp')
    commands=[];checks=[]
    for name,sanitizer in SANITIZERS.items():
        directory=out/name;directory.mkdir()
        flags=['clang++','-std=c++11','-O1','-g','-pthread','-Wall','-Wextra','-Wpedantic',
            '-fno-omit-frame-pointer','-fno-sanitize-recover=all','-fsanitize='+sanitizer,
            '-I'+str(source/'inst/include'),'-I'+str(generated),'-I'+str(HERE.parent)]
        common=[]
        for path,defs in [(source/'src'/n,[]) for n in units]+[(source/'src/P2.cpp',['-DP2=P2Reference']),(HERE.parent/'check.cpp',[])]:
            target=directory/('P2Reference.o' if defs else path.stem+'.o')
            cmd=flags+defs+['-c',path,'-o',target];invoke(cmd,directory/(target.name+'.log'));commands.append(list(map(str,cmd)));common.append(target)
        for variant in ('baseline','head1','head2','head4','head8','head32','hybrid16'):
            if variant=='baseline':path=source/'src/P2.cpp';defs=[]
            elif variant=='hybrid16':path=HERE.parent/'P2Hybrid.cpp';defs=['-DP2_HYBRID_MULTIPLIER=16']
            else:path=generated/'P2Initial.cpp';defs=['-DP2_INITIAL_CHUNK='+variant[4:]]
            target=directory/('check-'+variant)
            cmd=flags+defs+[path]+common+['-o',target]
            invoke(cmd,directory/('build-'+variant+'.log'));commands.append(list(map(str,cmd)))
            result=invoke([target],directory/('run-'+variant+'.log'));commands.append([str(target)])
            match=re.fullmatch(r'small_checks=(\d+)\nPASS checks=(\d+)\n',result)
            if not match:raise RuntimeError('unexpected common check output '+result)
            checks.append(dict(sanitizer=name,variant=variant,checks=int(match[2]),small_checks=int(match[1])))
            print(name,variant,result.strip(),flush=True)
        target=directory/'check-queue'
        cmd=flags+[HERE/'queue-check.cpp','-o',target]
        invoke(cmd,directory/'build-queue.log');commands.append(list(map(str,cmd)))
        result=invoke([target],directory/'run-queue.log');commands.append([str(target)])
        match=re.fullmatch(r'PASS queue_checks=(\d+)\n',result)
        if not match:raise RuntimeError('unexpected queue check output '+result)
        checks.append(dict(sanitizer=name,variant='queue',checks=int(match[1])))
        print(name,result.strip(),flush=True)
    # Generator rejects a source mutation and refuses output overwrite without writes.
    import importlib.util
    spec=importlib.util.spec_from_file_location('initial_generate',HERE/'generate.py');gen=importlib.util.module_from_spec(spec);spec.loader.exec_module(gen)
    rejection_checks=0
    for text in ((source/'src/P2.cpp').read_text()+'\n',(source/'src/P2.cpp').read_text().replace('CHUNK_SIZE = 32','CHUNK_SIZE = 16')):
        try:gen.generate(text)
        except RuntimeError:rejection_checks+=1
        else:raise RuntimeError('source drift accepted')
    before={n:sha(generated/n) for n in ('P2Initial.cpp','InitialChunks.h','generation.json')}
    reject=subprocess.run(['python3',str(HERE/'generate.py'),str(source),str(generated)],capture_output=True,text=True)
    if not reject.returncode or 'refusing to overwrite' not in reject.stderr:raise RuntimeError('overwrite accepted')
    if before!={n:sha(generated/n) for n in before}:raise RuntimeError('overwrite modified existing outputs')
    rejection_checks+=1
    paths=[HERE/n for n in ('InitialChunks.h','generate.py','queue-check.cpp','validate.py')]+[HERE.parent/'check.cpp',HERE.parent/'P2Hybrid.cpp',generated/'P2Initial.cpp']+[source/'src'/n for n in units+('P2.cpp',)]
    summary=dict(checks=checks,generator_rejections=rejection_checks,commands=commands,
        input_sha256={str(path):sha(path) for path in paths},
        scope='Common independent/differential checks unchanged; pure ticket queue boundary/overflow/concurrency/no-barrier checks added. No timings collected.')
    (out/'validation.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('saved',out/'validation.json',flush=True)

if __name__=='__main__':main()
