#!/usr/bin/env python3
"""Exact-copy instrumentation; never edit production or earlier experiments."""
from pathlib import Path
import argparse

def once(text, old, new):
    if text.count(old) != 1:
        raise RuntimeError('expected one source anchor: '+old[:100])
    return text.replace(old, new)

def instrument(text, hybrid):
    text = '#include "trace.h"\n' + text
    text = once(text, 'const std::int64_t quotient = x / primes[i];',
                'MP_QUERY(i);\n            const std::int64_t quotient = x / primes[i];')
    text = once(text, 'batch.launch([&, slot] {',
                'batch.launch([&, slot] {\n                MP_WORKER(slot);')
    text = once(text, 'subtotal += countP2Range(x, begin,',
                'MP_CHUNK(begin, std::min(begin + CHUNK_SIZE, last));\n                    subtotal += countP2Range(x, begin,')
    if not hybrid:
        text = once(text, 'if (workers == 1) return sum + countP2Range(x, first, last, tables);',
                    'MP_PHASE(0);\n        if (workers == 1) { MP_WORKER(0); MP_CHUNK(first,last); return sum + countP2Range(x, first, last, tables); }')
        return text.replace('std::int64_t P2(', 'std::int64_t P2Production(')
    text = once(text, 'if (workers == 1) return countP2Range(x, first, last, tables);',
                'MP_PHASE(Matched::current->group);\n        if (workers == 1) { MP_WORKER(0); MP_CHUNK(first,last); return countP2Range(x, first, last, tables); }')
    text = once(text, 'void fill() {', 'void fill() {\n            MP_FILL();')
    text = once(text, 'PrimeCountStream stream(start, finish, piAt(start, tables), tables);', '''
#ifdef MATCHED_TRACE
        const double seedBegin = Matched::now();
#endif
#ifdef MATCHED_COUNTERS
        MatchedPhi::Scope seedStats(Matched::worker->stats[2]);
#endif
        const auto seed = piAt(start, tables);
#ifdef MATCHED_TRACE
        Matched::worker->seed += Matched::now() - seedBegin;
        const double setupBegin = Matched::now();
#endif
        PrimeCountStream stream(start, finish, seed, tables);
#ifdef MATCHED_TRACE
        Matched::worker->setup += Matched::now() - setupBegin;
#endif
''')
    text = once(text, 'if (first == last) return 0;\n        const auto& primes = tables.primes();',
                'if (first == last) return 0;\n        MP_PHASE(1);\n        const auto& primes = tables.primes();')
    text = once(text, 'if (workers == 1) {\n#ifdef P2_HYBRID_MEMORY',
                'if (workers == 1) {\n            MP_WORKER(0);\n#ifdef P2_HYBRID_MEMORY')
    text = once(text, 'batch.launch([&, slot, begin, last, lower] {',
                'batch.launch([&, slot, begin, last, lower] {\n                    MP_WORKER(slot);')
    return text

def main():
    p=argparse.ArgumentParser(); p.add_argument('source',type=Path); p.add_argument('hybrid',type=Path); p.add_argument('output',type=Path)
    a=p.parse_args(); a.output.mkdir(exist_ok=True)
    (a.output/'P2Production.cpp').write_text(instrument((a.source/'src/P2.cpp').read_text(),False))
    (a.output/'P2Instrumented.cpp').write_text(instrument(a.hybrid.read_text(),True))

if __name__=='__main__': main()
