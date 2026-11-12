#!/usr/bin/env python3
"""Plot retained measurements using base R; no counting kernels are run."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import statistics
import subprocess

p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args()
with gzip.open(a.root/'query-profile.jsonl.gz','rt') as f:
    rows=[json.loads(s) for s in f]
rows=[r for r in rows if r['x']==80000000000 and r['cutoff']=='zero' and r['threads']==2 and r['repeat']>0]
duration=lambda r:sum(q['end_us']-q['begin_us'] for q in r['queries'])
blocks={b:[r for r in rows if r['block']==b] for b in sorted({r['block'] for r in rows})}
medians={b:statistics.median(duration(r) for r in v) for b,v in blocks.items()}
middle=statistics.median(medians.values());b=min(medians,key=lambda b:abs(medians[b]-middle))
row=min(blocks[b],key=lambda r:abs(duration(r)-medians[b]))
with (a.root/'plot-query.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=('index','prime','phi_ms'));w.writeheader()
    for q in row['queries']:w.writerow(dict(index=q['index'],prime=q['prime'],phi_ms=(q['end_us']-q['begin_us'])/1000))
(a.root/'plot-selection.json').write_text(json.dumps(dict(
    rule='Closest block median to median of warmed first32 sums, then closest warmed call to that block median.',
    case=row['case'],block=row['block'],repeat=row['repeat'],first32_us=duration(row)),indent=2)+'\n')
subprocess.run(['Rscript',str(Path(__file__).with_name('plot.R')),str(a.root)],check=True)
