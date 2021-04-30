#ifndef CONSTRAINTS_UTILS_H
#define CONSTRAINTS_UTILS_H

#include "Constraints/UserConstraintFuns.h"
#include "Partitions/PartitionsUtils.h"
#include "SetUpUtils.h"
#include <limits>
#include <map>

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

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
                     int n, int m, std::vector<std::string> &compFunVec,
                     const std::string &mainFun, VecType &myType,
                     SEXP Rtarget, SEXP RcompFun, SEXP Rtolerance,
                     SEXP Rlow, bool IsConstrained, bool bCalcMultiset);

#endif
