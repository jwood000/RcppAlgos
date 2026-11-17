# Frozen P2 test references

These small fixtures preserve the approved behavior without retaining the old
experimental variants and generated source trees.

- `P2Reference.cpp`: the earlier Phi-only P2 implementation, used as a
  differential reference for larger inputs. It shares Phi with production and
  is **not** an independent mathematical oracle.
- `P2SelectorCandidate.h`: the accepted streaming-selection policy, used to
  detect changes at numeric and worker-budget boundaries.
- `InitialChunks.h`: the original unsigned-ticket queue, used alongside an
  independently constructed expected range list to check the signed queue.

The common checker also uses trial division and an independently sieved pi table;
the private stream checker uses its own sieve. Those checks provide independent
coverage beyond agreement with these fixtures. All three files came from the
studies preserved at commit `36303eb268288a2449e348ba6d5b911ffb2fac51`.
