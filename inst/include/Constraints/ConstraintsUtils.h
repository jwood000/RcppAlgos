#ifndef CONSTRAINTS_UTILS_H
#define CONSTRAINTS_UTILS_H

#include "Constraints/ConstraintsTypes.h"
#include "Constraints/UserConstraintFuns.h"
#include "Partitions/PartitionsUtils.h"
#include "SetUpUtils.h"
#include <limits>
#include <map>

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

constexpr double dblIntMax = std::numeric_limits<int>::max();

// This fixes inequalites where the symbols are swapped. That is: "=>" amd "=<"
static const std::map<std::string, std::string> compForms = {
    {"<",   "<"}, {">",   ">"},
    {"<=", "<="}, {">=", ">="},
    {"==", "=="}, {"=<", "<="},
    {"=>", ">="}
};

static const std::array<std::string, 5> compSpecial = {
    {"==", ">,<", ">=,<", ">,<=", ">=,<="}
};

static const std::array<std::string, 5> compHelper = {
    {"<=", "<", "<", "<=", "<="}
};

// in R console: print(sqrt(.Machine$double.eps), digits = 16)
// [1] 0.00000001490116119384766
// Which is also 2^(-26)
constexpr double defaultTolerance = std::numeric_limits<float>::epsilon() / 8.0;

template <typename T>
void PopulateVec(const std::vector<T> &v,
                 std::vector<T> &cnstrntVec,
                 std::vector<int> &z, int &count,
                 int m, int nRows, bool IsComb);

template <typename T>
void SectionOne(const std::vector<T> &v, std::vector<T> &testVec,
                std::vector<int> &z, const std::vector<T> &targetVals,
                std::vector<T> &cnstrntVec, std::vector<T> &resVec,
                bool &check_0, bool &check_1, int &count,
                partialPtr<T> partial, funcPtr<T> fun,
                compPtr<T> compOne, compPtr<T> compTwo, int m, int m1,
                int nRows, int maxZ, bool IsComb, bool xtraCol);

bool CheckSpecialCase(bool bLower, const std::string &mainFun,
                      const std::vector<double> &vNum);

bool CheckIsInteger(const std::string &funPass, int n,
                    int m, const std::vector<double> &vNum,
                    const std::vector<double> &targetVals,
                    const funcPtr<double> myFunDbl,
                    bool checkLim = false);

void ConstraintSetup(const std::vector<double> &vNum,
                     const std::vector<int> &Reps,
                     std::vector<double> &targetVals,
                     std::vector<int> &targetIntVals,
                     const funcPtr<double> funDbl, PartDesign &part,
                     ConstraintType &ctype, int n, int m,
                     std::vector<std::string> &compFunVec,
                     const std::string &mainFun, VecType &myType,
                     SEXP Rtarget, SEXP RcompFun, SEXP Rtolerance,
                     SEXP Rlow, bool IsComb, bool bCalcMultiset);

#endif
