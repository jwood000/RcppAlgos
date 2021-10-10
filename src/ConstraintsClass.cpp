#include "Constraints/ConstraintsClass.h"

template <typename T>
void Constraints<T>::PrepareConstraints() {
    
}

template <typename T>
Constraints<T>::Constraints(
    SEXP Rv, int Rm, SEXP RcompRows, const std::vector<int> &bVec,
    const std::vector<int> &Rreps, const std::vector<int> &Rfreqs,
    const std::vector<int> &RvInt, const std::vector<double> &RvNum,
    VecType typePass, int RmaxThreads, SEXP RnumThreads, bool Rparallel,
    const PartDesign &Rpart, const std::vector<std::string> &RcompVec,
    std::vector<double> &RtarVals, std::vector<int> &RtarIntVals,
    std::vector<int> &RstartZ, const std::string &RmainFun,
    funcPtr<double> RfunDbl, ConstraintType Rctype, int RstrtLen,
    int Rcap, bool RKeepRes, bool RnumUnknown, double RcnstrtRows,
    mpz_t RcnstrtRowsMpz
) : ComboRes(Rv, Rm, RcompRows, bVec, Rreps, Rfreqs, RvInt, RvNum, typePass,
             RmaxThreads, RnumThreads, Rparallel, Rpart, RcompVec, RtarVals,
             RtarIntVals, RstartZ, RmainFun, RfunDbl, Rctype, RstrtLen, Rcap,
             RKeepRes, RnumUnknown, RcnstrtRows, RcnstrtRowsMpz),
    nextCnstrnt(GetCnstrtPtr<T>(IsMult, IsRep), m2(m - 2), nMinusM(n - m),
    freqsSize(std::accumulate(Rreps.cbegin(), Rreps.cend(), 0)),
    pentExtreme(freqsSize - m), fun(GetFuncPtr<T>(RmainFun))) {
    
    check_0 = true;
    check_1 = true;
}

template <typename T>
void Constraints<T>::startOver() {
    Combo::startOver();
}

template <typename T>
SEXP Constraints<T>::nextComb() {
    
    if (CheckEqSi(IsGmp, mpzIndex, dblIndex, 0)) {
        increment(IsGmp, mpzIndex, dblIndex);
        return VecReturn();
    } else if (CheckIndLT(IsGmp, mpzIndex, dblIndex,
                          cnstrtCountMpz, cnstrtCount)) {
        increment(IsGmp, mpzIndex, dblIndex);
        nextCnstrnt(v, targetVals, freqs, zIndex, testVec, z, fun, comp, m,
                    m1, m2, nMinusM, maxZ, pentExtreme, check_0, check_1);
        return VecReturn();
    } else if (CheckEqInd(IsGmp, mpzIndex, dblIndex,
                          cnstrtCountMpz, cnstrtCount)) {
        const std::string message = "No more results. To see the last "
                                    "result, use the prevIter method(s)\n\n";
        Rprintf(message.c_str());
        increment(IsGmp, mpzIndex, dblIndex);
        return Rf_ScalarLogical(false);
    } else {
        return Rf_ScalarLogical(false);
    }
}

template <typename T>
SEXP Constraints<T>::nextNumCombs(SEXP RNum) {
    
    int num;
    CleanConvert::convertPrimitive(RNum, num, VecType::Integer,
                                   "The number of results");
    
    if (CheckIndLT(IsGmp, mpzIndex, dblIndex,
                   cnstrtCountMpz, cnstrtCount)) {
        
        int nRows = 0;
        int numIncrement = 0;
        
        if (IsGmp) {
            mpz_sub(mpzTemp, cnstrtCountMpz, mpzIndex);
            nRows = mpz_cmp_si(mpzTemp, num) < 0 ? mpz_get_si(mpzTemp) : num;
            numIncrement = mpz_cmp_si(mpzTemp, num) < 0 ? (nRows + 1) : nRows;
        } else {
            dblTemp = cnstrtCount - dblIndex;
            nRows = num > dblTemp ? dblTemp : num;
            numIncrement = num > dblTemp ? (nRows + 1) : nRows;
        }
        
        if (CheckGrTSi(IsGmp, mpzIndex, dblIndex, 0)) {
            nextCnstrnt(v, targetVals, freqs, zIndex, testVec,
                        z, fun, comp, m, m1, m2, nMinusM, maxZ,
                        pentExtreme, check_0, check_1);
        }
        
        SEXP res = PROTECT(MatrixReturn(nRows));
        increment(IsGmp, mpzIndex, dblIndex, numIncrement);
        zUpdateIndex(vNum, vInt, z, sexpVec, res, width, nRows);
        UNPROTECT(1);
        return res;
    } else if (CheckEqInd(IsGmp, mpzIndex, dblIndex,
                          cnstrtCountMpz, cnstrtCount)) {
        const std::string message = "No more results. To see the last result"
                                    ", use the prevIter method(s)\n\n";
        Rprintf(message.c_str());
        increment(IsGmp, mpzIndex, dblIndex);
        return Rf_ScalarLogical(false);
    } else {
        return Rf_ScalarLogical(false);
    }
}

template <typename T>
SEXP Constraints<T>::nextGather() {
    
    if (IsGmp) {
        mpz_sub(mpzTemp, cnstrtCountMpz, mpzIndex);
        
        if (mpz_cmp_si(mpzTemp, std::numeric_limits<int>::max()) > 0) {
            Rf_error("The number of requested rows is greater than ",
                     std::to_string(std::numeric_limits<int>::max()).c_str());
        }
    } else {
        dblTemp = cnstrtCount - dblIndex;
        
        if (dblTemp > std::numeric_limits<int>::max()) {
            Rf_error("The number of requested rows is greater than ",
                     std::to_string(std::numeric_limits<int>::max()).c_str());
        }
    }
    
    const int nRows = (IsGmp) ? mpz_get_si(mpzTemp) : dblTemp;
    
    if (nRows > 0) {
        if (CheckGrTSi(IsGmp, mpzIndex, dblIndex, 0)) {
            nextCnstrnt(freqs, z, n1, m1);
        }
        
        SEXP res = PROTECT(MatrixReturn(nRows));
        
        if (IsGmp) {
            mpz_add_ui(mpzIndex, cnstrtCountMpz, 1u);
        } else {
            dblIndex = cnstrtCount + 1;
        }
        
        zUpdateIndex(vNum, vInt, z, sexpVec, res, width, nRows);
        UNPROTECT(1);
        return res;
    } else {
        return Rf_ScalarLogical(false);
    }
}

template <typename T>
SEXP Constraints<T>::currComb() {
    
    if (CheckIndGrT(IsGmp, mpzIndex, dblIndex,
                    cnstrtCountMpz, cnstrtCount)) {
        const std::string message = "No more results. To see the last "
                                    "result, use the prevIter method(s)\n\n";
        Rprintf(message.c_str());
        return Rf_ScalarLogical(0);
    } else if (CheckGrTSi(IsGmp, mpzIndex, dblIndex, 0)) {
        return VecReturn();
    } else {
        const std::string message = "Iterator Initialized. To see the first "
                                    "result, use the nextIter method(s)\n\n";
        Rprintf(message.c_str());
        return Rf_ScalarLogical(0);
    }
}

template <typename T>
SEXP Constraints<T>::summary() {
    SEXP res = PROTECT(Combo::summary());
    std::string desc(R_CHAR(STRING_ELT(VECTOR_ELT(res, 0), 0)));
    desc += " with " + mainFun + " applied to each result";
    SET_VECTOR_ELT(res, 0, Rf_mkString(desc.c_str()));
    UNPROTECT(1);
    return res;
}
