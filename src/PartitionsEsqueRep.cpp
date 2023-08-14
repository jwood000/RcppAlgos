#include "Constraints/PartitionsEsqueRep.h"

template <typename T>
int PartitionsEsqueRep<T>::GetLowerBound(
    const std::vector<T> &v, std::vector<int> &z,
    const funcPtr<T> fun, const reducePtr<T> reduce,
    const partialPtr<T> partial, T currPartial,
    int n, int m, int strt
) {

    const int lastElem = n - 1;
    const int lastCol  = m - 1;

    std::vector<T> vPass(m);
    std::fill(vPass.begin(), vPass.end(), v.back());
    T partVal = fun(vPass, m - 1);

    if (strt == 0) {
        const T testMax = partial(partVal, vPass.back(), m);

        if (testMax < tarMin) {
            return 0;
        }
    }

    if (strt) {
        for (int i = 0; i < strt; ++i) {
            vPass[i] = v[z[i]];
            partVal = partial(partVal, vPass[i], m);
            reduce(m, partVal, v[lastElem]);
        }

        currPartial = fun(vPass, strt);

        for (int i = strt; i < m; ++i) {
            vPass[i] = v[z[strt - 1]];
        }
    } else {
        std::fill(vPass.begin(), vPass.end(), v[0]);
    }

    const T testMin = fun(vPass, m);

    if (testMin > tarMax) {
        return 0;
    }

    int idx = lastElem;
    int lowBnd = (strt) ? z[strt - 1] : 0;

    for (int i = strt; i < lastCol; ++i) {
        if (this->LowerBound(v, tarMin, partVal, idx, lowBnd)) {
            if (idx > lowBnd) {
                const int numIterLeft = m - i;

                for (int j = 0; j < numIterLeft; ++j) {
                    vPass[j] = v[idx];
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

        lowBnd = idx;
        idx = lastElem;
        reduce(m, partVal, v[lastElem]);
    }

    this->LowerBoundLast(v, tarMin, partVal, idx, lowBnd);
    z[lastCol] = idx;
    return 1;
}

template <typename T>
void PartitionsEsqueRep<T>::NextSection(
    const std::vector<T> &v, const std::vector<T> &targetVals,
    std::vector<T> &testVec, std::vector<int> &z,
    const funcPtr<T> f, const compPtr<T> comp,
    int m, int m1, int m2
) {

    for (int i = m2; i >= 0 && !this->check_0; --i) {
        if (z[i] != this->maxZ) {
            ++z[i];
            testVec[i] = v[z[i]];

            GetLowerBound(v, z, f, reduce, this->partial,
                          currPartial, this->n, m, i + 1);

            for (int k = i + 1; k < m; ++k) {
                testVec[k] = v[z[k]];
            }

            T testVal = f(testVec, m);
            this->check_0 = comp(testVal, targetVals);
        }
    }
}

template <typename T>
void PartitionsEsqueRep<T>::Prepare(const std::string &currComp,
                                    std::vector<T> &v) {

    this->SetComparison(currComp);
    std::sort(v.begin(), v.end());
    this->z.assign(this->m, 0);

    this->check_1 = GetLowerBound(
        v, this->z, this->fun, reduce,
        this->partial, currPartial, this->n, this->m
    );
}

template <typename T>
PartitionsEsqueRep<T>::PartitionsEsqueRep(
    const std::vector<std::string> &comparison,
    const std::string &myFun, const std::string &myFunTest,
    int n_, int m_, bool IsComb_, bool xtraCol_,
    const std::vector<T> &targetVals
) : ConstraintsClass<T>(comparison, myFun, myFunTest,
                        n_, m_, IsComb_, xtraCol_),
    tarMin(*std::min_element(targetVals.cbegin(), targetVals.cend())),
    tarMax(*std::max_element(targetVals.cbegin(), targetVals.cend())),
    currPartial(myFun == "prod" ? 1 : 0), reduce(GetReducePtr<T>(myFun)) {}

template class PartitionsEsqueRep<int>;
template class PartitionsEsqueRep<double>;
