#pragma once

#include "Constraints/ConstraintsClass.h"

template <typename T>
class PartitionsEsqueRep : public ConstraintsClass<T> {
private:
    const T tarMin;
    const T tarMax;
    const T currPartial;
    const reducePtr<T> reduce;

    int GetLowerBound(const std::vector<T> &v, std::vector<int> &z,
                      const funcPtr<T> fun, const reducePtr<T> reduce,
                      const partialPtr<T> partial, T currPartial,
                      int n, int m, int strt = 0);

    void NextSection(
        const std::vector<T> &v, const std::vector<T> &targetVals,
        std::vector<T> &testVec, std::vector<int> &z,
        const funcPtr<T> f, const compPtr<T> comp,
        int m, int m1, int m2
    );

public:
    PartitionsEsqueRep(
        const std::vector<std::string> &comparison,
        const std::string &myFun, const std::string &myFunTest,
        int n_, int m_, bool IsComb_, bool xtraCol_,
        const std::vector<T> &targetVals
    );

    void Prepare(const std::string &currComp, std::vector<T> &v);
};
