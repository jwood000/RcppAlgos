#include "Partitions/NextPartition.h"

static bool keepGoing(const std::vector<int> &rpsCnt,
                      const std::vector<int> &z, int edge, int boundary) {

    if (edge >= 0) {
        const int myDiff = z[boundary] - z[edge];

        if (myDiff < 2) {
            return false;
        } else if (myDiff == 2) {
            return (rpsCnt[z[edge] + 1] > 1);
        } else {
            return (rpsCnt[z[edge] + 1] && rpsCnt[z[boundary] - 1]);
        }
    } else {
        return false;
    }
}

static std::vector<int> rleCpp(const std::vector<int> &x, int first_idx) {

    std::vector<int> lengths;
    int prev = x[first_idx];
    std::size_t i = 0;
    lengths.push_back(1);

    for(auto it = x.cbegin() + first_idx + 1; it != x.cend(); ++it) {
        if (prev == *it) {
            ++lengths[i];
        } else {
            lengths.push_back(1);
            prev = *it;
            ++i;
        }
    }

    return lengths;
}

static double NumPermsWithRep(const std::vector<int> &v, bool includeZero) {

    mpz_class result = 1;

    int first_idx = includeZero ? 0 : std::distance(
        v.cbegin(),
        std::find_if(v.cbegin(), v.cend(), [](int i) {return i != 0;})
    );

    // If all entries are zero or v is empty. This shouldn't happen,
    // but here for safety.
    if (first_idx == static_cast<int>(v.size())) return 1.0;

    std::vector<int> myLens = rleCpp(v, first_idx);
    std::sort(myLens.begin(), myLens.end(), std::greater<int>());

    const int myMax = myLens[0];
    const int numUni = myLens.size();

    for (int i = v.size() - first_idx; i > myMax; --i) {
        result *= i;
    }

    if (numUni > 1) {
        mpz_class div(1);

        for (int i = 1; i < numUni; ++i) {
            div *= mpz_class::factorial(myLens[i]);
        }

        mpz_divexact(result.get_mpz_t(), result.get_mpz_t(), div.get_mpz_t());
    }

    return result.get_d();
}

double CountPartsMultisetBrute(
    const std::vector<int> &Reps, const std::vector<int> &pz,
    bool IsComp, bool IsWeak
) {

    std::vector<int> z(pz.cbegin(), pz.cend());
    std::vector<int> rpsCnt(Reps.cbegin(), Reps.cend());

    const int lastCol  = pz.size() - 1;
    const int lastElem = Reps.size() - 1;

    int p = 0;
    int e = 0;
    int b = 0;

    // If we have made it here, we know a solution exists
    // (i.e. part.solnExists = true). The function keepGoing works by
    // terminating when it can no longer generate new partitions. The
    // current partition still counts hence why we start count at 1.
    PrepareMultisetPart(rpsCnt, z, b, p, e, lastCol, lastElem);
    double count = IsComp ? 0 : 1;

    // NOTE: z is always a partition produced by:
    //
    //     PrepareMultisetPart/NextMultisetGenPart
    //
    // i.e. non decreasing (lex-order generator). No need to sort before
    // calling NumPermsWithRep().
    if (IsComp && IsWeak) {
        for (; keepGoing(rpsCnt, z, e, b);
            NextMultisetGenPart(rpsCnt, z, e, b, p, lastCol, lastElem)) {
            count += NumPermsWithRep(z, true);
        }

        count += NumPermsWithRep(z, true);
    } else if (IsComp) {
        for (; keepGoing(rpsCnt, z, e, b);
            NextMultisetGenPart(rpsCnt, z, e, b, p, lastCol, lastElem)) {
            count += NumPermsWithRep(z, false);
        }

        count += NumPermsWithRep(z, false);
    } else {
        for (; keepGoing(rpsCnt, z, e, b);
            NextMultisetGenPart(rpsCnt, z, e, b, p, lastCol, lastElem)) {
            ++count;
        }
    }

    return count;
}
