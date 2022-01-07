#include "Constraints/ConstraintsDistinct.h"

template <typename T>
void ConstraintsDistinct<T>::NextSection(
        const std::vector<T> &v, const std::vector<T> &targetVals,
        std::vector<T> &testVec, std::vector<int> &z,
        const funcPtr<T> f, const compPtr<T> comp,
        int m, int m1, int m2, bool check_0, bool &check_1
    ) {

    if (check_1) {
        bool noChange = true;

        for (int i = m2; i >= 0 && !check_0; --i) {
            if (z[i] != (nMinusM + i)) {
                ++z[i];
                testVec[i] = v[z[i]];

                for (int k = i + 1; k < m; ++k) {
                    z[k] = z[k - 1] + 1;
                    testVec[k] = v[z[k]];
                }

                T testVal = f(testVec, m);
                check_0 = comp(testVal, targetVals);
                noChange = false;
            }
        }

        check_1 = (!noChange && check_0);
    }
}

template <typename T>
void ConstraintsDistinct<T>::Prepare(const std::string &currComp,
                                     std::vector<T> &v) {

    this->SetComparison(currComp);

    if (currComp == ">" || currComp == ">=") {
        std::sort(v.begin(), v.end(), std::greater<T>());
    } else {
        std::sort(v.begin(), v.end());
    }

    std::iota(this->z.begin(), this->z.end(), 0);
}

template <typename T>
ConstraintsDistinct<T>::ConstraintsDistinct(
    const std::vector<std::string> &comparison,
    const std::string &myFun, int n_, int m_,
    bool IsComb_, bool xtraCol_
) : ConstraintsClass<T>(comparison, myFun, n_, m_, IsComb_, xtraCol_),
    nMinusM(n_ - m_) {}

template class ConstraintsDistinct<int>;
template class ConstraintsDistinct<double>;
