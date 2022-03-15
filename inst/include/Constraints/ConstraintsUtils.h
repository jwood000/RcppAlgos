#ifndef CONSTRAINTS_UTILS_H
#define CONSTRAINTS_UTILS_H

#include "Constraints/ConstraintsTypes.h"
#include "Constraints/UserConstraintFuns.h"
#include "Partitions/PartitionsUtils.h"
#include "SetUpUtils.h"
#include <limits>
#include <map>

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
void AddResultToParts(T* mat, std::int64_t result,
                      std::size_t numResult,
                      std::size_t width);

template <typename T>
void VectorToMatrix(const std::vector<T> &cnstrntVec,
                    const std::vector<T> &resVec, T* mat,
                    std::int64_t result, std::size_t numResult,
                    std::size_t width, std::size_t upperBound,
                    bool xtraCol, bool IsPart);

bool CheckSpecialCase(bool bLower, const std::string &mainFun,
                      const std::vector<double> &vNum);

bool CheckIsInteger(const std::string &funPass, int n,
                    int m, const std::vector<double> &vNum,
                    const std::vector<double> &targetVals,
                    const funcPtr<double> myFunDbl, bool checkLim,
                    bool IsRep, bool IsMult, bool IsPart);

void ConstraintSetup(const std::vector<double> &vNum,
                     const std::vector<int> &Reps,
                     std::vector<double> &targetVals,
                     std::vector<int> &vInt, std::vector<int> &targetIntVals,
                     const funcPtr<double> funDbl, PartDesign &part,
                     ConstraintType &ctype, int lenV, int m,
                     std::vector<std::string> &compFunVec,
                     const std::string &mainFun, const std::string &funTest,
                     VecType &myType, SEXP Rtarget, SEXP RcompFun,
                     SEXP Rtolerance, SEXP Rlow,
                     bool IsComb, bool bCalcMulti = false);

#endif
