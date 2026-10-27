#include "Partitions/PartitionsCountMultiset.h"
#include "Partitions/MultisetCountEngine.h"

double CountPartsMultiset(
    int n, int m, const std::vector<int>& allowed,
    const std::vector<int>& Reps, int strtLen
) {
    MultisetCountEngine<false> counter;
    return counter.Count(n, m, allowed, Reps);
}

double CountCompsMultiset(
    int n, int m, const std::vector<int>& allowed,
    const std::vector<int>& Reps, int strtLen
) {
    MultisetCountEngine<true> counter;
    return counter.Count(n, m, allowed, Reps);
}

double CountCompsMultisetZero(
    int n, int m, const std::vector<int>& allowed,
    const std::vector<int>& Reps, int strtLen
) {

    if (strtLen == 0) {
        // This means that z contains only zeros
        return 1;
    } else {
        MultisetCountEngine<true> counter;
        return counter.Count(n, m, allowed, Reps, strtLen);
    }
}

double CountCompsMultisetWeak(
    int n, int m, const std::vector<int>& allowed,
    const std::vector<int>& Reps, int strtLen
) {

    if (strtLen == 0) {
        // This means that z contains only zeros
        return 1;
    } else {
        MultisetCountEngine<true> counter;
        return counter.Count(n, m, allowed, Reps);
    }
}
