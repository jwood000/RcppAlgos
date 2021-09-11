#include "ClassUtils/ComboApplyClass.h"
#include "ClassUtils/ExposeClass.h"
#include "ClassUtils/ComboClass.h"

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

    const std::vector<int> bVec = CleanConvert::GetNumVec<int>(RboolVec);
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
    
    const std::vector<int> bVec = CleanConvert::GetNumVec<int>(RboolVec);
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
