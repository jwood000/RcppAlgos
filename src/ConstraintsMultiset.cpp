#include "Constraints/ConstraintsMultiset.h"

template <typename T>
void ConstraintsMultiset<T>::NextSection(
        const std::vector<T> &v, const std::vector<T> &targetVals,
        std::vector<T> &testVec, std::vector<int> &z,
        const funcPtr<T> f, const compPtr<T> comp,
        int m, int m1, int m2
    ) {

    for (int i = m2; i >= 0 && !this->check_0; --i) {
        if (z[i] != freqs[pentExtreme + i]) {
            ++z[i];
            testVec[i] = v[z[i]];

            for (int j = i + 1, k = zIndex[z[i]] + 1; j < m; ++j, ++k) {
                z[j] = freqs[k];
                testVec[j] = v[z[j]];
            }

            T testVal = f(testVec, m);
            this->check_0 = comp(testVal, targetVals);
        }
    }
}

template <typename T>
void ConstraintsMultiset<T>::Prepare(const std::string &currComp,
                                     std::vector<T> &v) {

    this->SetComparison(currComp);

    if (currComp == ">" || currComp == ">=") {
        for (int i = 0; i < (this->n - 1); ++i) {
            for (int j = i + 1; j < this->n; ++j) {
                if (v[i] < v[j]) {
                    std::swap(v[i], v[j]);
                    std::swap(Reps[i], Reps[j]);
                }
            }
        }
    } else {
        for (int i = 0; i < (this->n - 1); ++i) {
            for (int j = i + 1; j < this->n; ++j) {
                if (v[i] > v[j]) {
                    std::swap(v[i], v[j]);
                    std::swap(Reps[i], Reps[j]);
                }
            }
        }
    }

    this->z.clear();
    zIndex.clear();
    freqs.clear();

    for (int i = 0, k = 0; i < this->n; ++i) {
        zIndex.push_back(k);

        for (int j = 0; j < Reps[i]; ++j, ++k) {
            freqs.push_back(i);
        }
    }

    this->z.assign(freqs.cbegin(), freqs.cbegin() + this->m);
}

template <typename T>
ConstraintsMultiset<T>::ConstraintsMultiset(
    const std::vector<std::string> &comparison,
    const std::string &myFun, const std::string &myFunTest,
    int n_, int m_, bool IsComb_, bool xtraCol_,
    std::vector<int> &Reps_
) : ConstraintsClass<T>(comparison, myFun, myFunTest,
                        n_, m_, IsComb_, xtraCol_),
    freqsSize(std::accumulate(Reps_.cbegin(), Reps_.cend(), 0)),
    pentExtreme(freqsSize - m_),
    Reps(Reps_) {}

template class ConstraintsMultiset<int>;
template class ConstraintsMultiset<double>;
