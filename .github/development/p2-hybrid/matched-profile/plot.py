#!/usr/bin/env python3
"""Select typical warmed traces by elapsed-time medians, then plot with base R."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import statistics
import subprocess

p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args()
with gzip.open(a.root/'timeline.jsonl.gz','rt') as f: rows=[json.loads(s) for s in f]
records=[];selected=[]
for cutoff in ('zero','meissel'):
    for mode in ('baseline','hybrid'):
        pool=[r for r in rows if r['cutoff']==cutoff and r['mode']==mode and r['repeat']>0]
        blocks={b:[r for r in pool if r['block']==b] for b in sorted({r['block'] for r in pool})}
        med={b:statistics.median(r['elapsed_us'] for r in v) for b,v in blocks.items()}
        target=statistics.median(med.values())
        block=min(med,key=lambda b:abs(med[b]-target))
        row=min(blocks[block],key=lambda r:abs(r['elapsed_us']-med[block]))
        selected.append({k:row[k] for k in ('cutoff','mode','block','repeat','elapsed_us')})
        for w in row['workers']:
            if w['group']==1:
                events=[(w['begin_us'],w['end_us'],'stream')]
            else:
                events=[(s,e,'prefix' if last<=1725 else 'tail' if first>=1725 else 'boundary')
                        for first,last,s,e in w['chunks']]
            for begin,end,kind in events:
                records.append(dict(cutoff=cutoff,mode=mode,group=w['group'],slot=w['slot'],
                    begin_ms=begin/1000,end_ms=end/1000,kind=kind,total_ms=row['elapsed_us']/1000))
with (a.root/'plot-traces.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
(a.root/'plot-selection.json').write_text(json.dumps(dict(
    rule='Closest block median to median of warmed process medians, then closest warmed call to that block median; ties retain original order.',
    selected=selected),indent=2)+'\n')
subprocess.run(['Rscript',str(Path(__file__).with_name('plot.R')),str(a.root)],check=True)
