#include "Constraints/CnstrntsSpecialClass.h"
#include "ComboGroups/ComboGroupsClass.h"
#include "Constraints/CnstrntsToRClass.h"
#include "Partitions/PartitionsClass.h"
#include "ClassUtils/ComboApplyClass.h"
#include "Cartesian/CartesianClass.h"
#include "ClassUtils/ComboResClass.h"
#include "ClassUtils/ComboClass.h"
#include "CheckReturn.h"

static void Finalizer(SEXP ext) {

    if (NULL == R_ExternalPtrAddr(ext)) {
        return;
    }

    class Iterator* ptr = (class Iterator*) R_ExternalPtrAddr(ext);
    R_ClearExternalPtr(ext);
    if (ptr) delete ptr;
}

[[cpp11::register]]
SEXP CartClassNew(SEXP Rv_RList, SEXP RNumThreads,
                  SEXP RmaxThreads, SEXP RForce_DF) {

    cpp11::list RList(Rv_RList);
    const int nCols = RList.size();;
    std::vector<std::vector<int>> myVec(nCols);
    std::vector<int> typeCheck(N_TYPES, 0);

    std::vector<int> IsFactor(nCols);
    std::vector<int> lenGrps(nCols);
    CartesianInitialPrep(RList, IsFactor, lenGrps, nCols);

    const int sumLength = std::accumulate(
        lenGrps.begin(), lenGrps.end(), 0
    );

    cpp11::writable::strings charVec(sumLength);
    std::vector<Rcomplex> cmplxVec(sumLength);
    std::vector<Rbyte> rawVec(sumLength);
    std::vector<double> dblVec(sumLength);
    std::vector<int> intVec(sumLength);
    std::vector<int> boolVec(sumLength);

    VecType myType = VecType::Integer;
    bool Force_DF = CppConvert::convertFlag(RForce_DF, "Return_DF");
    bool IsDF = true;

    ProductPrepare(
        RList, IsFactor, lenGrps, myVec, charVec, cmplxVec, rawVec,
        dblVec, intVec, boolVec, typeCheck, myType, nCols, IsDF
    );

    IsDF = IsDF || Force_DF;
    int maxThreads = 1;
    CppConvert::convertPrimitive(RmaxThreads, maxThreads,
                                 VecType::Integer, "maxThreads");

    const double computedRows = CartesianCount(lenGrps);
    const bool IsGmp = computedRows > Significand53;

    mpz_class computedRowsMpz;

    if (IsGmp) {
        CartesianCountGmp(computedRowsMpz, lenGrps);
    }

    const int maxLen = *std::max_element(lenGrps.begin(), lenGrps.end());
    std::vector<int> mat_idx(maxLen * nCols);

    // transposing myVec so that each row represents
    // each element of the given list
    for (int i = 0; i < nCols; ++i) {
        for (int j = 0; j < lenGrps[i]; ++j) {
            mat_idx[i + j * nCols] = myVec[i][j];
        }
    }

    // Transform lenGrps to be used in nextProduct
    for (auto &v_i: lenGrps) {
        v_i = nCols * (v_i - 1);
    }

    cpp11::sexp compRows = CppConvert::GetCount(
        IsGmp, computedRowsMpz, computedRows
    );

    class CartesianClass* ptr = new CartesianClass(
        Rv_RList, compRows, maxThreads, RNumThreads, false, IsGmp,
        mat_idx, typeCheck, IsFactor, intVec, dblVec, boolVec,
        cmplxVec, rawVec, charVec, lenGrps, IsDF, nCols, myType
    );

    cpp11::sexp ext = R_MakeExternalPtr(ptr, R_NilValue, R_NilValue);
    R_RegisterCFinalizerEx(ext, Finalizer, TRUE);

    return ext;
}

// RVals contains: v, vNum, vInt, m, RcompRows, maxThreads, nThreads,
//                 nGrps, grpSizes, & retType
//
// RboolVec contains: IsFac, IsComb, IsMult, IsRep,
//                    IsGmp, IsFull, IsComp, & IsWeak
//
// freqInfo contains: myReps & freqs

[[cpp11::register]]
SEXP CombClassNew(SEXP RVals, SEXP RboolVec, SEXP freqInfo, SEXP Rparallel,
                  SEXP RstdFun, SEXP Rrho, SEXP R_RFunVal, SEXP RmainFun,
                  SEXP RcompFun, SEXP Rtarget, SEXP RKeepRes,
                  SEXP Rtolerance, SEXP RmIsNull, SEXP RretVal) {

    const int ReturnValue = Rf_asInteger(RretVal);
    const std::vector<int> bVec =
        CppConvert::GetVec<int>(RboolVec);
    const std::vector<int> myReps =
        CppConvert::GetVec<int>(VECTOR_ELT(freqInfo, 0));
    const std::vector<int> freqs =
        CppConvert::GetVec<int>(VECTOR_ELT(freqInfo, 1));

    const int m = Rf_asInteger(VECTOR_ELT(RVals, 3));
    const int maxThreads = Rf_asInteger(VECTOR_ELT(RVals, 5));

    const std::vector<double> vNum =
        CppConvert::GetVec<double>(VECTOR_ELT(RVals, 1));
    std::vector<int> vInt =
        CppConvert::GetVec<int>(VECTOR_ELT(RVals, 2));

    VecType myType;
    SetType(myType, VECTOR_ELT(RVals, 0));
    const bool Parallel = CppConvert::convertFlag(Rparallel, "Parallel");

    if (ReturnValue == 1) {
        class Combo* ptr = new Combo(
            VECTOR_ELT(RVals, 0), m, VECTOR_ELT(RVals, 4), bVec, myReps,
            freqs, vInt, vNum, myType, maxThreads, VECTOR_ELT(RVals, 6),
            Parallel
        );

        cpp11::sexp ext = R_MakeExternalPtr(ptr, R_NilValue, R_NilValue);
        R_RegisterCFinalizerEx(ext, Finalizer, TRUE);

        return ext;
    } else if (ReturnValue == 2) {
        class ComboApply* ptr = new ComboApply(
            VECTOR_ELT(RVals, 0), m, VECTOR_ELT(RVals, 4), bVec, myReps,
            freqs, vInt, vNum, myType, maxThreads, VECTOR_ELT(RVals, 6),
            Parallel, RstdFun, Rrho, R_RFunVal
        );

        cpp11::sexp ext = R_MakeExternalPtr(ptr, R_NilValue, R_NilValue);
        R_RegisterCFinalizerEx(ext, Finalizer, TRUE);

        return ext;
    } else if (ReturnValue == 4) {
        class ComboGroupsClass* ptr = new ComboGroupsClass(
            VECTOR_ELT(RVals, 0), m, VECTOR_ELT(RVals, 4), bVec,
            myReps, freqs, vInt, vNum, myType, maxThreads,
            VECTOR_ELT(RVals, 6), Parallel, VECTOR_ELT(RVals, 7),
            VECTOR_ELT(RVals, 8), VECTOR_ELT(RVals, 9)
        );

        cpp11::sexp ext = R_MakeExternalPtr(ptr, R_NilValue, R_NilValue);
        R_RegisterCFinalizerEx(ext, Finalizer, TRUE);

        return ext;
    } else {
        const bool KeepRes = CppConvert::convertFlag(RKeepRes, "keepResults");
        const bool IsConstrained = CheckConstrnd(RmainFun, RcompFun, Rtarget);
        const int n = vNum.size();

        if (!Rf_isString(RmainFun) || Rf_length(RmainFun) != 1) {
            cpp11::stop("contraintFun must be one of the following:"
                         " 'prod', 'sum', 'mean', 'max', or 'min'");
        }

        const std::string funTest(CHAR(STRING_ELT(RmainFun, 0)));
        const auto funIt = std::find(mainFunSet.begin(), mainFunSet.end(), funTest);

        if (funIt == mainFunSet.end()) {
            cpp11::stop("contraintFun must be one of the following:"
                         " 'prod', 'sum', 'mean', 'max', or 'min'");
        }

        const std::string mainFun = funTest == "mean" ? "sum" : funTest;
        funcPtr<double> funDbl = GetFuncPtr<double>(mainFun);

        const bool IsComp   = static_cast<bool>(bVec[6]);
        const bool IsComb   = static_cast<bool>(bVec[1]) && !IsComp;
        const bool IsMult   = static_cast<bool>(bVec[2]);
        const bool IsRep    = static_cast<bool>(bVec[3]);
        const bool IsStdGmp = static_cast<bool>(bVec[4]);

        ConstraintType ctype = ConstraintType::NoConstraint;
        PartDesign part;

        part.isRep   = IsRep;
        part.isMult  = IsMult;
        part.mIsNull = static_cast<bool>(Rf_asLogical(RmIsNull));
        part.isComp  = IsComp;
        part.isComb  = IsComb;
        part.isWeak  = static_cast<bool>(bVec[7]);

        std::vector<std::string> compVec;
        std::vector<double> tarVals;
        std::vector<int> tarIntVals;

        if (IsConstrained) {
            ConstraintSetup(vNum, myReps, tarVals, vInt, tarIntVals,
                            funDbl, part, ctype, n, m, compVec, mainFun,
                            funTest, myType, Rtarget, RcompFun,
                            Rtolerance, R_NilValue, true);
        }

        mpz_class computedRowsMpz;

        if (IsStdGmp) {
            CppConvert::convertMpzClass(
                VECTOR_ELT(RVals, 4), computedRowsMpz, "computedRowsMpz"
            );
        }

        const bool usePartCount = part.isPart &&
                                  !part.isGmp &&
                                  !part.numUnknown;

        const double computedRows = usePartCount ? part.count :
            (IsStdGmp ? computedRowsMpz.get_d() :
                 Rf_asReal(VECTOR_ELT(RVals, 4)));
        const bool IsGmp = (computedRows > Significand53);

        if (IsGmp && part.isPart) {
            computedRowsMpz = part.bigCount;
        }

        // See comments in ConstraintsMain.cpp
        const bool numUnknown = ctype == ConstraintType::PartitionEsque ||
                                ctype == ConstraintType::SpecialCnstrnt ||
                                ctype == ConstraintType::General        ||
                                (part.isPart && part.numUnknown);

        std::vector<int> startZ(m);
        const int cap     = n - part.includeZero;
        const int strtLen = std::count_if(part.startZ.cbegin(),
                                          part.startZ.cend(),
                                          [](int i){return i > 0;});

        if (ctype < ConstraintType::PartMapping) {
            mpz_class zero(0);
            SetStartZ(myReps, freqs, startZ, IsComb, n, m,
                      0, zero, IsRep, IsMult, IsGmp);
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
                computedRows, computedRowsMpz
            );

            cpp11::sexp ext = R_MakeExternalPtr(ptr, R_NilValue, R_NilValue);
            R_RegisterCFinalizerEx(ext, Finalizer, TRUE);

            return ext;
        } else if (ctype == ConstraintType::General ||
                   ctype == ConstraintType::PartitionEsque) {

            class CnstrntsToR* ptr = new CnstrntsToR(
                VECTOR_ELT(RVals, 0), m, VECTOR_ELT(RVals, 4), bVec, myReps,
                freqs, vInt, vNum, myType, maxThreads, VECTOR_ELT(RVals, 6),
                Parallel, part, compVec, tarVals, tarIntVals, startZ, mainFun,
                funTest, funDbl, ctype, strtLen, cap, KeepRes, numUnknown,
                computedRows, computedRowsMpz
            );

            cpp11::sexp ext = R_MakeExternalPtr(ptr, R_NilValue, R_NilValue);
            R_RegisterCFinalizerEx(ext, Finalizer, TRUE);

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
                computedRows, computedRowsMpz
            );

            cpp11::sexp ext = R_MakeExternalPtr(ptr, R_NilValue, R_NilValue);
            R_RegisterCFinalizerEx(ext, Finalizer, TRUE);

            return ext;
        } else {
            class Partitions* ptr = new Partitions(
                VECTOR_ELT(RVals, 0), m, VECTOR_ELT(RVals, 4), bVec, myReps,
                freqs, vInt, vNum, myType, maxThreads, VECTOR_ELT(RVals, 6),
                Parallel, part, compVec, tarVals, tarIntVals, startZ, mainFun,
                funTest, funDbl, ctype, strtLen, cap, KeepRes, numUnknown,
                computedRows, computedRowsMpz
            );

            cpp11::sexp ext = R_MakeExternalPtr(ptr, R_NilValue, R_NilValue);
            R_RegisterCFinalizerEx(ext, Finalizer, TRUE);

            return ext;
        }
    }
}

[[cpp11::register]]
SEXP StartOverGlue(SEXP ext) {
    class Iterator* ptr = (class Iterator*) R_ExternalPtrAddr(ext);
    ptr->startOver();
    return R_NilValue;
}

[[cpp11::register]]
SEXP NextIterGlue(SEXP ext) {
    class Iterator* ptr = (class Iterator*) R_ExternalPtrAddr(ext);
    return ptr->nextIter();
}

[[cpp11::register]]
SEXP NextNumIterGlue(SEXP ext, SEXP Rnum) {
    class Iterator* ptr = (class Iterator*) R_ExternalPtrAddr(ext);
    return ptr->nextNumIters(Rnum);
}

[[cpp11::register]]
SEXP NextGatherGlue(SEXP ext) {
    class Iterator* ptr = (class Iterator*) R_ExternalPtrAddr(ext);
    return ptr->nextGather();
}

[[cpp11::register]]
SEXP PrevIterGlue(SEXP ext) {
    class Iterator* ptr = (class Iterator*) R_ExternalPtrAddr(ext);
    return ptr->prevIter();
}

[[cpp11::register]]
SEXP PrevNumIterGlue(SEXP ext, SEXP Rnum) {
    class Iterator* ptr = (class Iterator*) R_ExternalPtrAddr(ext);
    return ptr->prevNumIters(Rnum);
}

[[cpp11::register]]
SEXP PrevGatherGlue(SEXP ext) {
    class Iterator* ptr = (class Iterator*) R_ExternalPtrAddr(ext);
    return ptr->prevGather();
}

[[cpp11::register]]
SEXP CurrIterGlue(SEXP ext) {
    class Iterator* ptr = (class Iterator*) R_ExternalPtrAddr(ext);
    return ptr->currIter();
}

[[cpp11::register]]
SEXP RandomAccessGlue(SEXP ext, SEXP RIndexVec) {
    class Iterator* ptr = (class Iterator*) R_ExternalPtrAddr(ext);
    return ptr->randomAccess(RIndexVec);
}

[[cpp11::register]]
SEXP SourceVectorGlue(SEXP ext) {
    class Iterator* ptr = (class Iterator*) R_ExternalPtrAddr(ext);
    return ptr->sourceVector();
}

[[cpp11::register]]
SEXP FrontGlue(SEXP ext) {
    class Iterator* ptr = (class Iterator*) R_ExternalPtrAddr(ext);
    return ptr->front();
}

[[cpp11::register]]
SEXP BackGlue(SEXP ext) {
    class Iterator* ptr = (class Iterator*) R_ExternalPtrAddr(ext);
    return ptr->back();
}

[[cpp11::register]]
SEXP SummaryGlue(SEXP ext) {
    class Iterator* ptr = (class Iterator*) R_ExternalPtrAddr(ext);
    return ptr->summary();
}
