#pragma once

#include "Constraints/ConstraintsUtils.h"

template <typename T>
class ConstraintsClass {
protected:
    int maxZ;
    int count;

    const int n;
    const int m;
    const int m1;
    const int m2;
    const bool IsComb;
    const bool xtraCol;

    const FunType ftype;
    const FunType ftesttype;

    compPtr<T> compOne;
    compPtr<T> compTwo;

    const funcPtr<T> fun;
    const partialPtr<T> partial;

    bool check_0;
    bool check_1;

    std::vector<int> z;
    std::vector<T> testVec;

    double GetBound(double tarMin, double partVal);

    bool LowerBound(
        const std::vector<T> &v, T tarMin,
        T partVal, int &idx, int low
    );

    void LowerBoundLast(
        const std::vector<T> &v, T tarMin,
        T partVal, int &idx, int low
    );

    void PopulateVec(const std::vector<T> &v,
                     std::vector<T> &cnstrntVec, int limit);

    void FilterProspects(
        const std::vector<T> &v, const std::vector<T> &targetVals,
        std::vector<T> &cnstrntVec, std::vector<T> &resVec, int limit
    );

    virtual void NextSection(
        const std::vector<T> &v, const std::vector<T> &targetVals,
        std::vector<T> &testVec, std::vector<int> &z,
        const funcPtr<T> f, const compPtr<T> comp,
        int m, int m1, int m2
    ) = 0;

public:

    virtual ~ConstraintsClass() = default;
    ConstraintsClass(
        const std::vector<std::string> &comparison,
        const std::string &myFun, const std::string &myFunTest,
        int n_, int m_, bool IsComb_, bool xtraCol_
    );

    virtual void Prepare(const std::string &currComp, std::vector<T> &v) = 0;
    void SetComparison(const std::string &currComp);
    void GetSolutions(
        const std::vector<T> &v, const std::vector<T> &targetVals,
        std::vector<T> &cnstrntVec, std::vector<T> &resVec, int limit
    );

    int GetCount() const {return count;}
    void Reset() {count = 0;}
};

template <typename T>
std::unique_ptr<ConstraintsClass<T>> MakeConstraints(
    const std::vector<std::string> &comparison, const std::string &myFun,
    const std::string &myFunTest, std::vector<int> &Reps,
    const std::vector<T> &targetVals, ConstraintType ctype, int n,
    int m, bool IsComb, bool xtraCol, bool IsMult, bool IsRep
);
