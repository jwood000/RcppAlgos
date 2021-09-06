#include "ClassUtils/ExposeCombo.h"
#include "ClassUtils/ComboClass.h"

static void Finalizer(SEXP ext) {

    if (NULL == R_ExternalPtrAddr(ext)) {
        return;
    }

    Rprintf("finalizing\n");
    class Combo* ptr = (class Combo*) R_ExternalPtrAddr(ext);
    Free(ptr);
    R_ClearExternalPtr(ext);
}

// RVals contains: v, vNum, vInt, m, RcompRows, maxThreads, & nThreads
// RboolVec contains: IsFac, IsComb, IsMult, IsRep, IsGmp, & IsFull
// freqInfo contains: myReps & freqs
SEXP ComboNew(SEXP RVals, SEXP RboolVec, SEXP freqInfo) {
    
    std::vector<int> bVec;
    std::vector<int> Rreps;
    std::vector<int> Rfreqs;

    CleanConvert::convertVector(RboolVec, bVec, VecType::Integer,
                                "bVec", true, true, true);
    CleanConvert::convertVector(VECTOR_ELT(freqInfo, 0), Rreps,
                                VecType::Integer, "myReps");
    CleanConvert::convertVector(VECTOR_ELT(freqInfo, 1), Rfreqs,
                                VecType::Integer, "freqs");

    int maxThreads = 1;
    int Rm = Rf_asInteger(VECTOR_ELT(RVals, 3));
    
    CleanConvert::convertPrimitive(VECTOR_ELT(RVals, 5), maxThreads,
                                   VecType::Integer, "maxThreads");
    
    std::vector<int> RvInt;
    std::vector<double> RvNum;
    CleanConvert::convertVector(VECTOR_ELT(RVals, 1), RvInt,
                                VecType::Integer, "myReps");
    CleanConvert::convertVector(VECTOR_ELT(RVals, 2), RvNum,
                                VecType::Numeric, "vNum");
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
