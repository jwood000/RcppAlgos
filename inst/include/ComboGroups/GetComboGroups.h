#pragma once

#include "cpp11/function.hpp"
#include "cpp11/strings.hpp"
#include "cpp11/list.hpp"

#include "ComboGroups/ComboGroupsTemplate.h"
#include "CppConvert.h"
#include "SetUpUtils.h"
#include "RMatrix.h"
#include <functional>
#include <numeric>
#include <thread>

CmbGrpClsFuncs GetClassFuncs(
    std::unique_ptr<ComboGroupsTemplate> const &CmbGrp
);

SEXP GetComboGroups(
    SEXP Rv, nextGrpFunc nextCmbGrp, nthFuncDbl nthCmbGrp,
    nthFuncGmp nthCmbGrpGmp, finalTouchFunc FinalTouch,
    const std::vector<double> &vNum, const std::vector<int> &vInt,
    const std::vector<int> &startZ, const VecType &myType,
    const std::vector<double> &mySample, const std::vector<mpz_class> &myVec,
    mpz_class lowerMpz, double lower, int n, int numResults, int nThreads,
    bool IsArray, bool IsNamed, bool Parallel, bool IsSample, bool IsGmp
);
