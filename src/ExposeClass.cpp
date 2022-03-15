#include "Constraints/CnstrntsSpecialClass.h"
#include "Constraints/CnstrntsToRClass.h"
#include "Partitions/PartitionsClass.h"
#include "ClassUtils/ComboApplyClass.h"
#include "ClassUtils/ComboResClass.h"
#include "ClassUtils/ExposeClass.h"
#include "ClassUtils/ComboClass.h"
#include "CheckReturn.h"

static void Finalizer(SEXP ext) {

    if (NULL == R_ExternalPtrAddr(ext)) {
        return;
    }

    class Combo* ptr = (class Combo*) R_ExternalPtrAddr(ext);
    R_ClearExternalPtr(ext);
    if (ptr) delete ptr;
}

// RVals contains: v, vNum, vInt, m, RcompRows, maxThreads, & nThreads
// RboolVec contains: IsFac, IsComb, IsMult, IsRep, IsGmp, & IsFull
// freqInfo contains: myReps & freqs
SEXP CombClassNew(SEXP RVals, SEXP RboolVec, SEXP freqInfo, SEXP Rparallel,
                  SEXP RstdFun, SEXP Rrho, SEXP R_RFunVal, SEXP RmainFun,
                  SEXP RcompFun, SEXP Rtarget, SEXP RKeepRes,
                  SEXP Rtolerance, SEXP RmIsNull, SEXP RretVal) {

    const int ReturnValue = Rf_asInteger(RretVal);
    const std::vector<int> bVec   = CleanConvert::GetNumVec<int>(RboolVec);
    const std::vector<int> myReps = CleanConvert::GetNumVec<int>(
        VECTOR_ELT(freqInfo, 0)
    );

    const std::vector<int> freqs = CleanConvert::GetNumVec<int>(
        VECTOR_ELT(freqInfo, 1)
    );

    const int m = Rf_asInteger(VECTOR_ELT(RVals, 3));
    const int maxThreads = Rf_asInteger(VECTOR_ELT(RVals, 5));

    const std::vector<double> vNum = CleanConvert::GetNumVec<double>(
        VECTOR_ELT(RVals, 1)
    );
    std::vector<int> vInt = CleanConvert::GetNumVec<int>(
        VECTOR_ELT(RVals, 2)
    );

    VecType myType;
    SetType(myType, VECTOR_ELT(RVals, 0));
    const bool Parallel = CleanConvert::convertFlag(Rparallel, "Parallel");

    if (ReturnValue == 1) {
        class Combo* ptr = new Combo(
            VECTOR_ELT(RVals, 0), m, VECTOR_ELT(RVals, 4), bVec, myReps,
            freqs, vInt, vNum, myType, maxThreads, VECTOR_ELT(RVals, 6),
            Parallel
        );

        SEXP ext = PROTECT(R_MakeExternalPtr(ptr, R_NilValue, R_NilValue));
        R_RegisterCFinalizerEx(ext, Finalizer, TRUE);

        UNPROTECT(1);
        return ext;
    } else if (ReturnValue == 2) {
        class ComboApply* ptr = new ComboApply(
            VECTOR_ELT(RVals, 0), m, VECTOR_ELT(RVals, 4), bVec, myReps,
            freqs, vInt, vNum, myType, maxThreads, VECTOR_ELT(RVals, 6),
            Parallel, RstdFun, Rrho, R_RFunVal
        );

        SEXP ext = PROTECT(R_MakeExternalPtr(ptr, R_NilValue, R_NilValue));
        R_RegisterCFinalizerEx(ext, Finalizer, TRUE);

        UNPROTECT(1);
        return ext;
    } else {
        const bool KeepRes = CleanConvert::convertFlag(RKeepRes, "keepResults");
        const bool IsConstrained = CheckConstrnd(RmainFun, RcompFun, Rtarget);
        const int n = vNum.size();

        if (!Rf_isString(RmainFun) || Rf_length(RmainFun) != 1) {
            Rf_error("contraintFun must be one of the following:"
                         " 'prod', 'sum', 'mean', 'max', or 'min'");
        }

        const std::string funTest(CHAR(STRING_ELT(RmainFun, 0)));
        const auto funIt = std::find(mainFunSet.begin(), mainFunSet.end(), funTest);

        if (funIt == mainFunSet.end()) {
            Rf_error("contraintFun must be one of the following:"
                         " 'prod', 'sum', 'mean', 'max', or 'min'");
        }

        const std::string mainFun = funTest == "mean" ? "sum" : funTest;
        funcPtr<double> funDbl = GetFuncPtr<double>(mainFun);

        const bool IsComb   = static_cast<bool>(bVec[1]);
        const bool IsMult   = static_cast<bool>(bVec[2]);
        const bool IsRep    = static_cast<bool>(bVec[3]);
        const bool IsStdGmp = static_cast<bool>(bVec[4]);

        ConstraintType ctype = ConstraintType::NoConstraint;
        PartDesign part;

        part.isRep   = IsRep;
        part.isMult  = IsMult;
        part.mIsNull = static_cast<bool>(Rf_asLogical(RmIsNull));

        std::vector<std::string> compVec;
        std::vector<double> tarVals;
        std::vector<int> tarIntVals;

        if (IsConstrained) {
            ConstraintSetup(vNum, myReps, tarVals, vInt, tarIntVals,
                            funDbl, part, ctype, n, m, compVec, mainFun,
                            funTest, myType, Rtarget, RcompFun,
                            Rtolerance, R_NilValue, IsComb, true);
        }

        auto computedRowsMpz = FromCpp14::make_unique<mpz_t[]>(1);
        mpz_init(computedRowsMpz[0]);

        if (IsStdGmp) {
            createMPZArray(VECTOR_ELT(RVals, 4), computedRowsMpz.get(), 1,
                           "computedRowsMpz");
        }

        const double computedRows = (part.count > 0) ? part.count :
            (IsStdGmp ? mpz_get_d(computedRowsMpz[0]) :
                 Rf_asReal(VECTOR_ELT(RVals, 4)));
        const bool IsGmp = (computedRows > Significand53);

        if (IsGmp && part.isPart) {
            mpz_set(computedRowsMpz[0], part.bigCount);
        }

        // See comments in ConstraintsMain.cpp
        const bool numUnknown = ctype == ConstraintType::PartitionEsque ||
                                ctype == ConstraintType::SpecialCnstrnt ||
                                ctype == ConstraintType::General        ||
                                part.numUnknown;

        std::vector<int> startZ(m);
        const int cap     = n - part.includeZero;
        const int strtLen = std::count_if(part.startZ.cbegin(),
                                          part.startZ.cend(),
                                          [](int i){return i > 0;});

        if (ctype < ConstraintType::PartMapping) {
            mpz_t zero;
            mpz_init(zero);
            SetStartZ(myReps, freqs, startZ, IsComb, n, m,
                      0, zero, IsRep, IsMult, IsGmp);
            mpz_clear(zero);
        } else {
            startZ = part.startZ;
        }

        if (ctype == ConstraintType::NoConstraint) {
            funDbl = GetFuncPtr<double>(funTest);

            if (myType == VecType::Integer) {
                if (!CheckIsInteger(funTest, n, m, vNum, vNum, funDbl,
                                    false, IsRep, IsMult, false)) {
                    myType = VecType::Numeric;
                }
            }

            // The true literal below is for KeepRes. We don't set the variable
            // to true because it is declared as a const, which we want to keep
            class ComboRes* ptr = new ComboRes(
                VECTOR_ELT(RVals, 0), m, VECTOR_ELT(RVals, 4), bVec, myReps,
                freqs, vInt, vNum, myType, maxThreads, VECTOR_ELT(RVals, 6),
                Parallel, part, compVec, tarVals, tarIntVals, startZ, funTest,
                funTest, funDbl, ctype, strtLen, cap, true, numUnknown,
                computedRows, computedRowsMpz[0]
            );

            SEXP ext = PROTECT(R_MakeExternalPtr(ptr, R_NilValue, R_NilValue));
            R_RegisterCFinalizerEx(ext, Finalizer, TRUE);

            UNPROTECT(1);
            return ext;
        } else if (ctype == ConstraintType::General ||
                   ctype == ConstraintType::PartitionEsque) {

            class CnstrntsToR* ptr = new CnstrntsToR(
                VECTOR_ELT(RVals, 0), m, VECTOR_ELT(RVals, 4), bVec, myReps,
                freqs, vInt, vNum, myType, maxThreads, VECTOR_ELT(RVals, 6),
                Parallel, part, compVec, tarVals, tarIntVals, startZ, mainFun,
                funTest, funDbl, ctype, strtLen, cap, KeepRes, numUnknown,
                computedRows, computedRowsMpz[0]
            );

            SEXP ext = PROTECT(R_MakeExternalPtr(ptr, R_NilValue, R_NilValue));
            R_RegisterCFinalizerEx(ext, Finalizer, TRUE);

            UNPROTECT(1);
            return ext;
        } else if (ctype == ConstraintType::SpecialCnstrnt) {
            // We must use a single thread to ensure the proper constraint
            // algorithm is called here: ConstraintsSpecial.cpp. Note that
            // the multithreaded algo is constrained by the maximum number
            // of results (i.e. count is incremented in the main body and
            // is checked in the do while. In the single threaded case,
            // count is incremented when the condition is met and isn't
            // relied on in the do while).
            class CnstrntsSpecial* ptr = new CnstrntsSpecial(
                VECTOR_ELT(RVals, 0), m, VECTOR_ELT(RVals, 4), bVec, myReps,
                freqs, vInt, vNum, myType, maxThreads, Rf_ScalarInteger(1),
                Parallel, part, compVec, tarVals, tarIntVals, startZ, mainFun,
                funTest, funDbl, ctype, strtLen, cap, KeepRes, numUnknown,
                computedRows, computedRowsMpz[0]
            );

            SEXP ext = PROTECT(R_MakeExternalPtr(ptr, R_NilValue, R_NilValue));
            R_RegisterCFinalizerEx(ext, Finalizer, TRUE);

            UNPROTECT(1);
            return ext;
        } else {
            class Partitions* ptr = new Partitions(
                VECTOR_ELT(RVals, 0), m, VECTOR_ELT(RVals, 4), bVec, myReps,
                freqs, vInt, vNum, myType, maxThreads, VECTOR_ELT(RVals, 6),
                Parallel, part, compVec, tarVals, tarIntVals, startZ, mainFun,
                funTest, funDbl, ctype, strtLen, cap, KeepRes, numUnknown,
                computedRows, computedRowsMpz[0]
            );

            SEXP ext = PROTECT(R_MakeExternalPtr(ptr, R_NilValue, R_NilValue));
            R_RegisterCFinalizerEx(ext, Finalizer, TRUE);

            UNPROTECT(1);
            return ext;
        }
    }
}

SEXP StartOverGlue(SEXP ext) {
    class Combo* ptr = (class Combo*) R_ExternalPtrAddr(ext);
    ptr->startOver();
    return R_NilValue;
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
