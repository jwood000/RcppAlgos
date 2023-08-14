#include "Constraints/PartitionsEsqueDistinct.h"

template <typename T>
int PartitionsEsqueDistinct<T>::GetLowerBound(
    const std::vector<T> &v, std::vector<int> &z,
    const funcPtr<T> fun, const reducePtr<T> reduce,
    const partialPtr<T> partial, T currPartial,
    int n, int m, int strt
) {

    const int lastCol = m - 1;
    std::vector<T> vPass(m);
    vPass.assign(v.crbegin(), v.crbegin() + m);
    T partVal = fun(vPass, m - 1);

    if (strt == 0) {
        const T testMax = partial(partVal, vPass.back(), m);

        if (testMax < tarMin) {
            return 0;
        }
    }

    int currPos = n - m;

    if (strt) {
        for (int i = 0; i < strt; ++i) {
            vPass[i] = v[z[i]];
            partVal = partial(partVal, vPass[i], m);
            ++currPos;
            reduce(m, partVal, v[currPos]);
        }

        currPartial = fun(vPass, strt);

        for (int i = strt, j = 1; i < m; ++i, ++j) {
            vPass[i] = v[z[strt - 1] + j];
        }
    } else {
        vPass.assign(v.cbegin(), v.cbegin() + m);
    }

    const T testMin = fun(vPass, m);

    if (testMin > tarMax) {
        return 0;
    }

    int idx = n - m + strt;
    int lowBnd = (strt) ? z[strt - 1] + 1 : 0;

    for (int i = strt; i < lastCol; ++i) {
        if (this->LowerBound(v, tarMin, partVal, idx, lowBnd)) {
            if (idx > lowBnd) {
                const int numIterLeft = m - i;

                for (int j = 0, k = idx; j < numIterLeft; ++j, ++k) {
                    vPass[j] = v[k];
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

        ++idx;
        ++currPos;

        lowBnd = idx;
        idx = currPos;
        reduce(m, partVal, v[currPos]);
    }

    this->LowerBoundLast(v, tarMin, partVal, idx, lowBnd);
    z[lastCol] = idx;
    return 1;
}

template <typename T>
void PartitionsEsqueDistinct<T>::NextSection(
    const std::vector<T> &v, const std::vector<T> &targetVals,
    std::vector<T> &testVec, std::vector<int> &z,
    const funcPtr<T> f, const compPtr<T> comp,
    int m, int m1, int m2
) {

    for (int i = m2; i >= 0 && !this->check_0; --i) {
        if (z[i] != (nMinusM + i)) {
            ++z[i];
            testVec[i] = v[z[i]];

            GetLowerBound(v, z, f, reduce, this->partial,
                          currPartial, this->n, m, i + 1);

            for (int k = (i + 1); k < m; ++k) {
                testVec[k] = v[z[k]];
            }

            T testVal = f(testVec, m);
            this->check_0 = comp(testVal, targetVals);
        }
    }
}

template <typename T>
void PartitionsEsqueDistinct<T>::Prepare(const std::string &currComp,
                                         std::vector<T> &v) {

    this->SetComparison(currComp);
    std::sort(v.begin(), v.end());
    std::iota(this->z.begin(), this->z.end(), 0);

    this->check_1 = GetLowerBound(
        v, this->z, this->fun, reduce,
        this->partial, currPartial, this->n, this->m
    );
}

template <typename T>
PartitionsEsqueDistinct<T>::PartitionsEsqueDistinct(
    const std::vector<std::string> &comparison,
    const std::string &myFun, const std::string &myFunTest,
    int n_, int m_, bool IsComb_, bool xtraCol_,
    const std::vector<T> &targetVals
) : ConstraintsClass<T>(comparison, myFun, myFunTest,
                        n_, m_, IsComb_, xtraCol_),
    tarMin(*std::min_element(targetVals.cbegin(), targetVals.cend())),
    tarMax(*std::max_element(targetVals.cbegin(), targetVals.cend())),
    currPartial(myFun == "prod" ? 1 : 0), reduce(GetReducePtr<T>(myFun)),
    nMinusM(n_ - m_) {}

template class PartitionsEsqueDistinct<int>;
template class PartitionsEsqueDistinct<double>;
