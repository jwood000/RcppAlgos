#pragma once

#include "Constraints/ConstraintsClass.h"

template <typename T>
class ConstraintsMultiset : public ConstraintsClass<T> {
private:
    const int freqsSize;
    const int pentExtreme;

    std::vector<int> Reps;
    std::vector<int> freqs;
    std::vector<int> zIndex;

    void NextSection(
        const std::vector<T> &v, const std::vector<T> &targetVals,
        std::vector<T> &testVec, std::vector<int> &z,
        const funcPtr<T> f, const compPtr<T> comp,
        int m, int m1, int m2
    );

public:
    ConstraintsMultiset(
        const std::vector<std::string> &comparison,
        const std::string &myFun, const std::string &myFunTest,
        int n_, int m_, bool IsComb_, bool xtraCol_,
        std::vector<int> &Reps_
    );

    void Prepare(const std::string &currComp, std::vector<T> &v);
};
