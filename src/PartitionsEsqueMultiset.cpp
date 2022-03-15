#include "Constraints/PartitionsEsqueMultiset.h"

template <typename T>
int PartitionsEsqueMultiset<T>::GetLowerBound(
        const std::vector<T> &v, std::vector<int> &z,
        const funcPtr<T> fun, const reducePtr<T> reduce,
        const partialPtr<T> partial, T currPartial,
        int n, int m, int strt
    ) {

    const int lastCol = m - 1;
    const int lenMinusM = freqs.size() - m;

    std::vector<T> vPass(m);

    for (int i = freqs.size() - 1, j = 0; i >= lenMinusM; --i, ++j) {
        vPass[j] = v[freqs[i]];
    }

    T partVal = fun(vPass, m - 1);

    if (strt == 0) {
        const T testMax = partial(partVal, vPass.back(), m);

        if (testMax < tarMin) {
            return 0;
        }
    }

    int zExpCurrPos = freqs.size() - m;
    std::vector<int> repsCounter(Reps.cbegin(), Reps.cend());

    if (strt) {
        for (int i = 0; i < strt; ++i) {
            vPass[i] = v[z[i]];
            partVal = partial(partVal, vPass[i], m);
            --repsCounter[z[i]];
            ++zExpCurrPos;
            reduce(m, partVal, v[freqs[zExpCurrPos]]);
        }

        currPartial = fun(vPass, strt);

        if (z[strt - 1] != freqs.back()) {
            const auto it = std::find(freqs.begin(), freqs.end(), z[strt - 1] + 1);

            // Find the first index in freqs that equals z[strt - 1] + 1
            // We want to get the next index after z[strt - 1], so we must
            // take into account repsCounter, which keeps track of how many
            // of each index is left.
            const int myInd = std::distance(freqs.begin(), it);
            const int freqsStrt = myInd - repsCounter[z[strt - 1]];

            for (int i = strt, j = freqsStrt; i < m; ++i, ++j) {
                vPass[i] = v[freqs[j]];
            }
        } else {
            for (int i = strt; i < m; ++i) {
                vPass[i] = v[freqs.back()];
            }
        }
    } else {
        for (int i = 0; i < m; ++i) {
            vPass[i] = v[freqs[i]];
        }
    }

    const T testMin = fun(vPass, m);

    if (testMin > tarMax) {
        return 0;
    }

    int idx = freqs[freqs.size() - m + strt];
    int lowBnd = 0;

    if (strt) {
        lowBnd = repsCounter[z[strt - 1]] ? z[strt - 1] : z[strt - 1] + 1;
    }

    for (int i = strt; i < lastCol; ++i) {
        if (this->LowerBound(v, tarMin, partVal, idx, lowBnd)) {
            if (idx > lowBnd && repsCounter[idx - 1]) {
                const int numIterLeft = m - i;
                const auto it = std::find(freqs.begin(), freqs.end(), idx + 1);
                const int myInd = std::distance(freqs.begin(), it);
                const int freqsStrt = myInd - repsCounter[idx];

                for (int j = 0, k = freqsStrt; j < numIterLeft; ++j, ++k) {
                    vPass[j] = v[freqs[k]];
                }

                const T minRemaining = fun(vPass, numIterLeft);
                const T currMin = partial(minRemaining, currPartial, m);

                if (currMin > tarMin) {
                    --idx;
                }
            }
        }

        z[i] = idx;
        partVal = partial(partVal, v[idx], m);
        currPartial = partial(currPartial, v[idx], m);

        --repsCounter[idx];

        if (repsCounter[idx] == 0) {
            ++idx;
        }

        ++zExpCurrPos;
        lowBnd = idx;
        idx = freqs[zExpCurrPos];
        reduce(m, partVal, v[idx]);
    }

    this->LowerBoundLast(v, tarMin, partVal, idx, lowBnd);
    z[lastCol] = idx;
    return 1;
}

template <typename T>
void PartitionsEsqueMultiset<T>::NextSection(
        const std::vector<T> &v, const std::vector<T> &targetVals,
        std::vector<T> &testVec, std::vector<int> &z,
        const funcPtr<T> f, const compPtr<T> comp,
        int m, int m1, int m2
    ) {

    for (int i = m2; i >= 0 && !this->check_0; --i) {
        if (z[i] != freqs[pentExtreme + i]) {
            ++z[i];
            testVec[i] = v[z[i]];

            GetLowerBound(v, z, f, reduce, this->partial,
                          currPartial, this->n, m, i + 1);

            for (int j = i + 1, k = zIndex[z[i]] + 1; j < m; ++j, ++k) {
                testVec[j] = v[freqs[k]];
            }

            T testVal = f(testVec, m);
            this->check_0 = comp(testVal, targetVals);
        }
    }
}

template <typename T>
void PartitionsEsqueMultiset<T>::Prepare(const std::string &currComp,
                                         std::vector<T> &v) {

    this->SetComparison(currComp);

    // Since PartitionsEsque only occurs for "==" or "IsBetween" (See
    // ConstraintStructure  in ConstraintsUtils.cpp) we will not need
    // to sort in reverse order as we do in ConstraintsMultiset.cpp
    for (int i = 0; i < (this->n - 1); ++i) {
        for (int j = i + 1; j < this->n; ++j) {
            if (v[i] > v[j]) {
                std::swap(v[i], v[j]);
                std::swap(Reps[i], Reps[j]);
            }
        }
    }

    zIndex.clear();
    freqs.clear();

    for (int i = 0, k = 0; i < this->n; ++i) {
        zIndex.push_back(k);

        for (int j = 0; j < Reps[i]; ++j, ++k) {
            freqs.push_back(i);
        }
    }

    this->check_1 = GetLowerBound(
        v, this->z, this->fun, reduce, this->partial,
        currPartial, this->n, this->m, this->count
    );
}

template <typename T>
PartitionsEsqueMultiset<T>::PartitionsEsqueMultiset(
    const std::vector<std::string> &comparison,
    const std::string &myFun, const std::string &myFunTest,
    int n_, int m_, bool IsComb_, bool xtraCol_,
    const std::vector<T> &targetVals,
    std::vector<int> &Reps_
) : ConstraintsClass<T>(comparison, myFun, myFunTest,
                        n_, m_, IsComb_, xtraCol_),
    tarMin(*std::min_element(targetVals.cbegin(), targetVals.cend())),
    tarMax(*std::max_element(targetVals.cbegin(), targetVals.cend())),
    currPartial(myFun == "prod" ? 1 : 0), reduce(GetReducePtr<T>(myFun)),
    freqsSize(std::accumulate(Reps_.cbegin(), Reps_.cend(), 0)),
    pentExtreme(freqsSize - m_),
    Reps(Reps_) {}

template class PartitionsEsqueMultiset<int>;
template class PartitionsEsqueMultiset<double>;
