#include "Partitions/PartitionsMultiset.h"
#include "Partitions/NextPartition.h"

// The current approach is not ideal. We must iterate
int CountPartsMultiset(const std::vector<int> &Reps,
                       const std::vector<int> &pz) {

    std::vector<int> z(pz.cbegin(), pz.cend());
    std::vector<int> rpsCnt(Reps.cbegin(), Reps.cend());

    const int lastCol  = pz.size() - 1;
    const int lastElem = Reps.size() - 1;

    int p = 0;
    int e = 0;
    int b = 0;

    // If we have made it here, we know a solution exists
    // (i.e. part.solnExists = true). The way keepGoing works
    // is that it terminates when it can no longer generate
    // new partitions. The current partition still counts
    // hence why we start count at 1.
    int count = 1;
    PrepareMultisetPart(rpsCnt, z, b, p, e, lastCol, lastElem);

    for (; keepGoing(rpsCnt, lastElem, z, e, b);
         NextMultisetGenPart(rpsCnt, z, e, b, p, lastCol, lastElem)) {
        ++count;
    }

    return count;
}
