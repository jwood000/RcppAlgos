#include "ClassUtils/ConstraintsClass.h"
#include "ClassUtils/ComboApplyClass.h"
#include "ClassUtils/ExposeClass.h"
#include "ClassUtils/ComboClass.h"
#include "CheckReturn.h"

static void Finalizer(SEXP ext) {

    if (NULL == R_ExternalPtrAddr(ext)) {
        return;
    }

    class Combo* ptr = (class Combo*) R_ExternalPtrAddr(ext);
    Free(ptr);
    R_ClearExternalPtr(ext);
}

// RVals contains: v, vNum, vInt, m, RcompRows, maxThreads, & nThreads
// RboolVec contains: IsFac, IsComb, IsMult, IsRep, IsGmp, & IsFull
// freqInfo contains: myReps & freqs
SEXP ComboNew(SEXP RVals, SEXP RboolVec, SEXP freqInfo) {

    const std::vector<int> bVec  = CleanConvert::GetNumVec<int>(RboolVec);
    const std::vector<int> Rreps = CleanConvert::GetNumVec<int>(
        VECTOR_ELT(freqInfo, 0)
    );
    
    const std::vector<int> Rfreqs = CleanConvert::GetNumVec<int>(
        VECTOR_ELT(freqInfo, 1)
    );
    
    const int Rm = Rf_asInteger(VECTOR_ELT(RVals, 3));
    const int maxThreads = Rf_asInteger(VECTOR_ELT(RVals, 5));
    
    const std::vector<double> RvNum = CleanConvert::GetNumVec<double>(
        VECTOR_ELT(RVals, 1)
    );
    const std::vector<int> RvInt = CleanConvert::GetNumVec<int>(
        VECTOR_ELT(RVals, 2)
    );
    
    VecType myType;
    SetType(myType, VECTOR_ELT(RVals, 0));
    
    class Combo* ptr = new Combo(VECTOR_ELT(RVals, 0), Rm,
                                 VECTOR_ELT(RVals, 4), bVec, Rreps, Rfreqs,
                                 RvInt, RvNum, myType, maxThreads,
                                 VECTOR_ELT(RVals, 6));
    SEXP ext = PROTECT(R_MakeExternalPtr(ptr, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(ext, Finalizer, TRUE);
    
    UNPROTECT(1);
    return ext;
}

SEXP ComboApplyNew(SEXP RVals, SEXP RboolVec, SEXP freqInfo,
                   SEXP RstdFun, SEXP Rrho, SEXP R_RFunVal) {
    
    const std::vector<int> bVec  = CleanConvert::GetNumVec<int>(RboolVec);
    const std::vector<int> Rreps = CleanConvert::GetNumVec<int>(
        VECTOR_ELT(freqInfo, 0)
    );
    
    const std::vector<int> Rfreqs = CleanConvert::GetNumVec<int>(
        VECTOR_ELT(freqInfo, 1)
    );
    
    const int Rm = Rf_asInteger(VECTOR_ELT(RVals, 3));
    const int maxThreads = Rf_asInteger(VECTOR_ELT(RVals, 5));
    
    const std::vector<double> RvNum = CleanConvert::GetNumVec<double>(
        VECTOR_ELT(RVals, 1)
    );
    const std::vector<int> RvInt = CleanConvert::GetNumVec<int>(
        VECTOR_ELT(RVals, 2)
    );
    
    VecType myType;
    SetType(myType, VECTOR_ELT(RVals, 0));
    
    class ComboApply* ptr = new ComboApply(VECTOR_ELT(RVals, 0), Rm,
                                           VECTOR_ELT(RVals, 4), bVec, Rreps,
                                           Rfreqs, RvInt, RvNum, myType,
                                           maxThreads, VECTOR_ELT(RVals, 6),
                                           RstdFun, Rrho, R_RFunVal);
    SEXP ext = PROTECT(R_MakeExternalPtr(ptr, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(ext, Finalizer, TRUE);
    
    UNPROTECT(1);
    return ext;
}

void StartOverGlue(SEXP ext) {
    class Combo* ptr = (class Combo*) R_ExternalPtrAddr(ext);
    ptr->startOver();
}

SEXP NextCombGlue(SEXP ext) {
    class Combo* ptr = (class Combo*) R_ExternalPtrAddr(ext);
    return ptr->nextComb();
}

SEXP NextNumCombGlue(SEXP ext, SEXP Rnum) {
    class Combo* ptr = (class Combo*) R_ExternalPtrAddr(ext);
    return ptr->nextNumCombs(Rnum);
}

SEXP NextGatherGlue(SEXP ext) {
    class Combo* ptr = (class Combo*) R_ExternalPtrAddr(ext);
    return ptr->nextGather();
}

SEXP PrevCombGlue(SEXP ext) {
    class Combo* ptr = (class Combo*) R_ExternalPtrAddr(ext);
    return ptr->prevComb();
}

SEXP PrevNumCombGlue(SEXP ext, SEXP Rnum) {
    class Combo* ptr = (class Combo*) R_ExternalPtrAddr(ext);
    return ptr->prevNumCombs(Rnum);
}

SEXP PrevGatherGlue(SEXP ext) {
    class Combo* ptr = (class Combo*) R_ExternalPtrAddr(ext);
    return ptr->prevGather();
}

SEXP CurrCombGlue(SEXP ext) {
    class Combo* ptr = (class Combo*) R_ExternalPtrAddr(ext);
    return ptr->currComb();
}

SEXP RandomAccessGlue(SEXP ext, SEXP RIndexVec) {
    class Combo* ptr = (class Combo*) R_ExternalPtrAddr(ext);
    return ptr->randomAccess(RIndexVec);
}

SEXP SourceVectorGlue(SEXP ext) {
    class Combo* ptr = (class Combo*) R_ExternalPtrAddr(ext);
    return ptr->sourceVector();
}

SEXP FrontGlue(SEXP ext) {
    class Combo* ptr = (class Combo*) R_ExternalPtrAddr(ext);
    return ptr->front();
}

SEXP BackGlue(SEXP ext) {
    class Combo* ptr = (class Combo*) R_ExternalPtrAddr(ext);
    return ptr->back();
}

SEXP SummaryGlue(SEXP ext) {
    class Combo* ptr = (class Combo*) R_ExternalPtrAddr(ext);
    return ptr->summary();
}

SEXP ConstraintsNew(SEXP RVals, SEXP RboolVec, SEXP freqInfo, SEXP RmainFun,
                    SEXP RcompFun, SEXP Rlimits, SEXP RKeepRes, SEXP Rtarget,
                    SEXP Rtolerance, SEXP RmIsNull) {
    
    const std::vector<int> bVec  = CleanConvert::GetNumVec<int>(RboolVec);
    const std::vector<int> Rreps = CleanConvert::GetNumVec<int>(
        VECTOR_ELT(freqInfo, 0)
    );
    
    const std::vector<int> Rfreqs = CleanConvert::GetNumVec<int>(
        VECTOR_ELT(freqInfo, 1)
    );
    
    const int Rm = Rf_asInteger(VECTOR_ELT(RVals, 3));
    const int maxThreads = Rf_asInteger(VECTOR_ELT(RVals, 5));
    
    std::vector<double> RvNum = CleanConvert::GetNumVec<double>(
        VECTOR_ELT(RVals, 1)
    );
    std::vector<int> RvInt = CleanConvert::GetNumVec<int>(
        VECTOR_ELT(RVals, 2)
    );
    
    VecType myType;
    SetType(myType, VECTOR_ELT(RVals, 0));
    
    const bool IsConstrained = CheckConstrnd(RmainFun, RcompFun, Rtarget);
    const int n = RvNum.size();
    
    if (!Rf_isString(RmainFun) || Rf_length(RmainFun) != 1) {
        Rf_error("contraintFun must be one of the following:"
                     " 'prod', 'sum', 'mean', 'max', or 'min'");
    }
    
    const std::string mainFun(CHAR(STRING_ELT(RmainFun, 0)));
    const auto funIt = std::find(mainFunSet.begin(), mainFunSet.end(), mainFun);
    
    if (funIt == mainFunSet.end()) {
        Rf_error("contraintFun must be one of the following:"
                     " 'prod', 'sum', 'mean', 'max', or 'min'");
    }
    
    std::vector<int> tarIntVals;
    const funcPtr<double> funDbl = GetFuncPtr<double>(mainFun);
    const bool IsComb = static_cast<bool>(bVec[1]);
    const bool IsStdGmp = static_cast<bool>(bVec[4]);

    std::vector<std::string> compVec;
    std::vector<double> tarVals;

    ConstraintType ctype = ConstraintType::NoConstraint;
    PartDesign part;

    part.isRep = static_cast<bool>(bVec[3]);
    part.isMult = static_cast<bool>(bVec[2]);
    part.mIsNull = static_cast<bool>(Rf_asLogical(RmIsNull));

    if (IsConstrained) {
        ConstraintSetup(RvNum, Rreps, tarVals, RvInt, tarIntVals,
                        funDbl, part, ctype, n, Rm, compVec, mainFun, myType,
                        Rtarget, RcompFun, Rtolerance, R_NilValue, IsComb, false);
    }

    const double computedRows = (part.count > 0) ? part.count :
        (IsStdGmp ? 0 : Rf_asReal(VECTOR_ELT(RVals, 4)));

    const bool IsGmp = (computedRows > Significand53);
    mpz_t computedRowsMpz;
    mpz_init(computedRowsMpz);

    // if (IsGmp && part.isPart) {
    //     mpz_set(computedRowsMpz, part.bigCount);
    // } else if (IsGmp) {
    //     GetComputedRowMpz(computedRowsMpz, IsMult,
    //                       IsComb, IsRep, n, m, Rm, freqs, myReps);
    // }

    // This variable is used in determining the number of results. If the
    // output is constrained and the ConstraintType is "General" or
    // "PartitionEsque", it means we really don't know how many results
    // we have. The computedRows above is a strict upper bound but not
    // necessarily the least upper bound. In these cases, we don't want
    // to unnecessarily throw an error when computedRows exceeds 2^31 - 1.
    const bool numUnknown = ctype == ConstraintType::PartitionEsque ||
        ctype == ConstraintType::SpecialCnstrnt ||
        ctype == ConstraintType::General        ||
        part.numUnknown;
    // 
    // std::vector<int> startZ(m);
    // const int cap     = n - part.includeZero;
    // const int strtLen = std::count_if(part.startZ.cbegin(),
    //                                   part.startZ.cend(),
    //                                   [](int i){return i > 0;});
    // 
    // if (ctype < ConstraintType::PartMapping) {
    //     SetStartZ(myReps, freqs, startZ, IsComb, n, m,
    //               lower, lowerMpz[0], IsRep, IsMult, IsGmp);
    // } else {
    //     if (bLower) {
    //         nthPartsPtr nthPartFun = GetNthPartsFunc(part.ptype, IsGmp);
    //         startZ = nthPartFun(part.mapTar, part.width,
    //                             cap, strtLen, lower, lowerMpz[0]);
    //         
    //         if (ctype == ConstraintType::PartStandard && !part.includeZero) {
    //             for (auto &z_i: startZ) {
    //                 ++z_i;
    //             }
    //         }
    //     } else {
    //         startZ = part.startZ;
    //     }
    // }
    // 
    // // This is used when we are unable to calculate the number of results
    // // upfront (E.g. comboGeneral(rnorm(10), 5, constraintFun = "sum,
    // //                            comparisonFun = "<=", limitConstraints = 1))
    // double userNum = 0;
    // const bool bSetNum = !numUnknown ||
    //     ctype == ConstraintType::SpecialCnstrnt;
    // 
    // SetNumResults(IsGmp, bLower, bUpper, bSetNum, upperMpz.get(),
    //               lowerMpz.get(), lower, upper, computedRows,
    //               computedRowsMpz, nRows, userNum);
    // 
    // 
    // class ComboApply* ptr = new ComboApply(VECTOR_ELT(RVals, 0), Rm,
    //                                        VECTOR_ELT(RVals, 4), bVec, Rreps,
    //                                        Rfreqs, RvInt, RvNum, myType,
    //                                        maxThreads, VECTOR_ELT(RVals, 6),
    //                                        RstdFun, Rrho, R_RFunVal);
    // SEXP ext = PROTECT(R_MakeExternalPtr(ptr, R_NilValue, R_NilValue));
    // R_RegisterCFinalizerEx(ext, Finalizer, TRUE);
    // 
    // UNPROTECT(1);
    // return ext;
}
