#ifndef CONSTRAINTS_DISTINCT_H
#define CONSTRAINTS_DISTINCT_H

#include "Constraints/ConstraintsClass.h"

template <typename T>
class ConstraintsDistinct : public ConstraintsClass<T> {
private:
    const int nMinusM;

    void NextSection(
        const std::vector<T> &v, const std::vector<T> &targetVals,
        std::vector<T> &testVec, std::vector<int> &z,
        const funcPtr<T> f, const compPtr<T> comp,
        int m, int m1, int m2, bool check_0, bool &check_1
    );

public:
    ConstraintsDistinct(
        const std::vector<std::string> &comparison,
        const std::string &myFun, int n_, int m_,
        bool IsComb_, bool xtraCol_
    );

    void Prepare(const std::string &currComp, std::vector<T> &v);
};

#endif
