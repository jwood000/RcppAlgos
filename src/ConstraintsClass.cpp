#include "Constraints/PartitionsEsqueMultiset.h"
#include "Constraints/PartitionsEsqueDistinct.h"
#include "Constraints/PartitionsEsqueRep.h"
#include "Constraints/ConstraintsMultiset.h"
#include "Constraints/ConstraintsDistinct.h"
#include "Constraints/ConstraintsRep.h"
#include <memory>
#include <chrono>

// Used for checking whether user has interrupted computation
constexpr auto timeout = std::chrono::milliseconds(1000);

FunType GetFunType(const std::string &myFun) {

    if (myFun == "min") {
        return FunType::Min;
    } else if (myFun == "max") {
        return FunType::Max;
    } else if (myFun == "sum") {
        return FunType::Sum;
    } else if (myFun == "prod") {
        return FunType::Prod;
    } else {
        return FunType::Mean;
    }
}

template <typename T>
double ConstraintsClass<T>::GetBound(double tarMin, double partVal) {

    switch (ftype) {
        case FunType::Sum : {
            return tarMin - partVal;
        } case FunType::Prod : {
            return tarMin / partVal;
        } default : {
            return (tarMin * m) - (partVal * (m - 1));
        }
    }
}

template <typename T>
bool ConstraintsClass<T>::LowerBound(
        const std::vector<T> &v, T tarMin,
        T partVal, int &idx, int low
    ) {

    const double bound = GetBound(tarMin, partVal);

    if (v[idx] <= bound) {
        return false;
    } else if (v[low] < bound) {
        auto lower = std::find_if(
            v.cbegin() + low, v.cbegin() + idx,
            [=](T v_i) {return v_i >= bound;}
        );

        idx = std::distance(v.cbegin(), lower);
        return v[idx] > bound;
    } else {
        idx = low;
        return false;
    }
}

template <typename T>
void ConstraintsClass<T>::LowerBoundLast(
        const std::vector<T> &v, T tarMin,
        T partVal, int &idx, int low
    ) {

    const double bound = GetBound(tarMin, partVal);

    if (v[idx] > bound && v[low] < bound) {
        while (idx > low && v[idx] > bound) {
            --idx;
        }
    } else {
        idx = low;
    }
}

template <typename T>
void ConstraintsClass<T>::SetComparison(const std::string &currComp) {
    compOne = GetCompPtr<T>(currComp);
    compTwo = compOne;

    const auto itComp = std::find(compSpecial.cbegin(),
                                  compSpecial.cend(), currComp);

    if (itComp != compSpecial.end()) {
        const int myIndex = std::distance(compSpecial.cbegin(), itComp);
        compTwo = GetCompPtr<T>(compHelper[myIndex]);
    }

    testVec.assign(m, 0);
    check_0 = true;
    check_1 = true;
}

template <typename T>
void ConstraintsClass<T>::PopulateVec(
        const std::vector<T> &v, std::vector<T> &cnstrntVec, int limit
    ) {

    if (IsComb) {
        for (int k = 0; k < m; ++k) {
            cnstrntVec.push_back(v[z[k]]);
        }

        ++count;
    } else {
        do {
            for (int k = 0; k < m; ++k) {
                cnstrntVec.push_back(v[z[k]]);
            }

            ++count;
        } while (count < limit && std::next_permutation(z.begin(), z.end()));
    }
}

template <typename T>
void ConstraintsClass<T>::FilterProspects(
        const std::vector<T> &v, const std::vector<T> &targetVals,
        std::vector<T> &cnstrntVec, std::vector<T> &resVec, int limit
    ) {

    for (int i = 0; i < m; ++i) {
        testVec[i] = v[z[i]];
    }

    const T partialVal = fun(testVec, m1);
    T testVal = partial(partialVal, testVec.back(), m);
    check_0 = compTwo(testVal, targetVals);

    while (check_0 && check_1) {
        if (compOne(testVal, targetVals)) {
            const int myStart = count;
            PopulateVec(v, cnstrntVec, limit);

            for (int i = myStart; xtraCol && i < count; ++i) {
                if (ftesttype == FunType::Mean) {
                    resVec.push_back(testVal / m);
                } else {
                    resVec.push_back(testVal);
                }
            }

            check_1 = count < limit;
        }

        check_0 = z[m1] != maxZ;

        if (check_0) {
            ++z[m1];
            testVec[m1] = v[z[m1]];
            testVal = partial(partialVal, testVec.back(), m);
            check_0 = compTwo(testVal, targetVals);
        }
    }
}

template <typename T>
void ConstraintsClass<T>::GetSolutions(
        const std::vector<T> &v, const std::vector<T> &targetVals,
        std::vector<T> &cnstrntVec, std::vector<T> &resVec, int limit
    ) {

    check_1 = count < limit;

    if (m == 1) {
        int ind = 0;
        T testVal = v[ind];
        check_0 = compTwo(testVal, targetVals);

        while (check_0 && check_1) {
            if (compOne(testVal, targetVals)) {
                for (int k = 0; k < m; ++k) {
                    cnstrntVec.push_back(v[ind]);
                }

                ++count;
                check_1 = (count < limit);
                if (xtraCol) {resVec.push_back(testVal);}
            }

            check_0 = ind != maxZ;

            if (check_0) {
                ++ind;
                testVal = v[ind];
                check_0 = compTwo(testVal, targetVals);
            }
        }
    } else {
        auto check_point_1 = std::chrono::steady_clock::now();

        while (check_0 && check_1) {
            FilterProspects(v, targetVals, cnstrntVec, resVec, limit);
            NextSection(v, targetVals, testVec, z,
                        fun, compTwo, m, m1, m2);

            const auto check_point_2 = std::chrono::steady_clock::now();

            if (check_point_2 - check_point_1 > timeout) {
                cpp11::check_user_interrupt();
                check_point_1 = std::chrono::steady_clock::now();
            }
        }
    }
}

template <typename T>
ConstraintsClass<T>::ConstraintsClass(
    const std::vector<std::string> &comparison,
    const std::string &myFun, const std::string &myFunTest,
    int n_, int m_, bool IsComb_, bool xtraCol_
) : maxZ(n_ - 1), n(n_), m(m_), m1(m - 1),
    m2(m - 2), IsComb(IsComb_), xtraCol(xtraCol_),
    ftype(GetFunType(myFun)), ftesttype(GetFunType(myFunTest)),
    fun(GetFuncPtr<T>(myFun)), partial(GetPartialPtr<T>(myFun)) {

    z.assign(m, 0);
    testVec.assign(m, 0);
    count = 0;
}

template <typename T>
std::unique_ptr<ConstraintsClass<T>> MakeConstraints(
        const std::vector<std::string> &comparison, const std::string &myFun,
        const std::string &myFunTest, std::vector<int> &Reps,
        const std::vector<T> &targetVals, ConstraintType ctype, int n,
        int m, bool IsComb, bool xtraCol, bool IsMult, bool IsRep
    ) {

    if (ctype == ConstraintType::PartitionEsque) {
        if (IsMult) {
            return std::make_unique<PartitionsEsqueMultiset<T>>(
                comparison, myFun, myFunTest, n, m,
                IsComb, xtraCol, targetVals, Reps
            );
        } else if (IsRep) {
            return std::make_unique<PartitionsEsqueRep<T>>(
                comparison, myFun, myFunTest, n, m,
                IsComb, xtraCol, targetVals
            );
        } else {
            return std::make_unique<PartitionsEsqueDistinct<T>>(
                comparison, myFun, myFunTest, n, m,
                IsComb, xtraCol, targetVals
            );
        }
    } else if (IsMult) {
        return std::make_unique<ConstraintsMultiset<T>>(
            comparison, myFun, myFunTest, n, m, IsComb, xtraCol, Reps
        );
    } else if (IsRep) {
        return std::make_unique<ConstraintsRep<T>>(
            comparison, myFun, myFunTest, n, m, IsComb, xtraCol
        );
    } else {
        return std::make_unique<ConstraintsDistinct<T>>(
            comparison, myFun, myFunTest, n, m, IsComb, xtraCol
        );
    }
}

template class ConstraintsClass<int>;
template class ConstraintsClass<double>;

template std::unique_ptr<ConstraintsClass<int>> MakeConstraints(
        const std::vector<std::string>&, const std::string&,
        const std::string&, std::vector<int>&, const std::vector<int>&,
        ConstraintType, int, int, bool, bool, bool, bool
    );

template std::unique_ptr<ConstraintsClass<double>> MakeConstraints(
        const std::vector<std::string>&, const std::string&,
        const std::string&, std::vector<int>&, const std::vector<double>&,
        ConstraintType, int, int, bool, bool, bool, bool
    );
