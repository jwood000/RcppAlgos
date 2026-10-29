#pragma once

#include "Constraints/ConstraintsTypes.h"
#include "Constraints/UserConstraintFuns.h"
#include "Partitions/PartitionsUtils.h"
#include "SetUpUtils.h"
#include <limits>
#include <map>

// This fixes inequalities where the symbols are swapped. That is: "=>" and "=<"
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
constexpr double defaultTolerance =
    std::numeric_limits<float>::epsilon() / 8.0;

inline void ValidatePartitionDesign(const PartDesign &part) {
    // Generic constraint fallbacks do not preserve composition semantics.
    if (part.isComp &&
        (part.ptype == PartitionType::CoarseGrained ||
         part.ptype == PartitionType::NotPartition)) {
        cpp11::stop("Generation is not supported for this composition.");
    }

    if (part.isPart && part.ptype == PartitionType::NotMapped) {
        cpp11::stop("Currently, there is no algorithm for this case.");
    }
}

void InitialSetupPartDesign(
    PartDesign &part, SEXP RIsWeak, SEXP RIsComposition,
    bool IsRep, bool IsMult, bool mIsNull, bool IsComb
);

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

bool CheckSpecialCase(const std::vector<double> &vNum,
                      const std::string &mainFun, PartitionType ptype,
                      bool isPerm, bool bLower);

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
                     SEXP Rtolerance, SEXP Rlow);
