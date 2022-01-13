#include "Constraints/CnstrntsToRClass.h"

template <int sexpType, typename T>
SEXP CnstrtVecReturn(const std::vector<T> &v) {

    SEXP res = PROTECT(Rf_allocVector(sexpType, v.size()));

    if (sexpType == INTSXP) {
        int* ptrOut = INTEGER(res);

        for (std::size_t j = 0; j < v.size(); ++j) {
            ptrOut[j] = v[j];
        }
    } else {
        double* ptrOut = REAL(res);

        for (std::size_t j = 0; j < v.size(); ++j) {
            ptrOut[j] = v[j];
        }
    }

    UNPROTECT(1);
    return res;
}

template <typename T>
void GetNSolutions(const std::vector<std::string> &compVec,
                   std::unique_ptr<ConstraintsClass<T>> &Cnstrt,
                   std::vector<T> &cnstrntVec, std::vector<T> &resVec,
                   std::vector<T> &v, std::vector<T> &tar, int nRows) {

    int limit = Cnstrt->GetCount() + nRows;
    Cnstrt->GetSolutions(v, tar, cnstrntVec, resVec, limit);

    if (Cnstrt->GetCount() < limit && compVec.size() == 2 && tar.size() == 2) {
        tar.erase(tar.begin());
        Cnstrt->Prepare(compVec.back(), v);
        Cnstrt->GetSolutions(v, tar, cnstrntVec, resVec, limit);
    }
}

template <int sexpType, typename T>
SEXP GetNextCnstrt(const std::vector<std::string> &compVec,
                   std::unique_ptr<ConstraintsClass<T>> &Cnstrt,
                   std::vector<T> &v, std::vector<T> &tar,
                   std::vector<T> &curr, bool Keep, bool &keepGoing) {

    std::vector<T> resVec;
    std::vector<T> cnstrntVec;
    GetNSolutions(compVec, Cnstrt, cnstrntVec, resVec, v, tar, 1);

    if (cnstrntVec.size()) {
        if (Keep) cnstrntVec.push_back(resVec.front());
        curr = cnstrntVec;
        return CnstrtVecReturn<sexpType>(cnstrntVec);
    } else {
        keepGoing = false;
        const std::string message = "No more results.\n\n";
        Rprintf(message.c_str());
        return Rf_ScalarLogical(false);
    }
}

SEXP CnstrntsToR::GetNext() {

    if (RTYPE == INTSXP) {
        return GetNextCnstrt<INTSXP>(compVec, CnstrtInt, vInt, tarIntVals,
                                     currIntVec, KeepRes, keepGoing);
    } else {
        return GetNextCnstrt<REALSXP>(compVec, CnstrtDbl, vNum, tarVals,
                                      currDblVec, KeepRes, keepGoing);
    }
}

template <typename T>
void SetCurrVec(const std::vector<T> &cnstrntVec,
                const std::vector<T> &resVec,
                std::vector<T> &curr, std::size_t m, bool Keep) {

    std::vector<T> newCurr(cnstrntVec.end() - std::min(m, cnstrntVec.size()),
                           cnstrntVec.end());
    if (Keep) newCurr.push_back(resVec.back());
    curr = newCurr;
}

SEXP CnstrntsToR::GetNextN(int n) {

    if (RTYPE == INTSXP) {
        std::vector<int> resVec;
        std::vector<int> cnstrntVec;
        GetNSolutions(compVec, CnstrtInt, cnstrntVec, resVec, vInt, tarIntVals, n);

        if (cnstrntVec.size()) {
            SetCurrVec(cnstrntVec, resVec, currIntVec, width, KeepRes);
            const int vecLen = cnstrntVec.size();
            const int numResult = vecLen / m;

            SEXP res = PROTECT(Rf_allocMatrix(INTSXP, numResult, nCols));
            int* matInt = INTEGER(res);

            VectorToMatrix(cnstrntVec, resVec, matInt, 0, numResult,
                           width, upperBoundInt, KeepRes, false);
            UNPROTECT(1);
            return res;
        }
    } else {
        std::vector<double> resVec;
        std::vector<double> cnstrntVec;
        GetNSolutions(compVec, CnstrtDbl, cnstrntVec, resVec, vNum, tarVals, n);

        if (cnstrntVec.size()) {
            SetCurrVec(cnstrntVec, resVec, currDblVec, width, KeepRes);
            const int vecLen = cnstrntVec.size();
            const int numResult = vecLen / m;

            SEXP res = PROTECT(Rf_allocMatrix(REALSXP, numResult, nCols));
            double* matNum = REAL(res);

            VectorToMatrix(cnstrntVec, resVec, matNum, 0, numResult,
                           width, upperBoundDbl, KeepRes, false);
            UNPROTECT(1);
            return res;
        }
    }

    keepGoing = false;
    const std::string message = "No more results.\n\n";
    Rprintf(message.c_str());
    return Rf_ScalarLogical(false);
}

CnstrntsToR::CnstrntsToR(
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
    maxRows(std::min(dblIntMax, RcnstrtRows)),
    origTarIntVals(RtarIntVals), origTarVals(RtarVals) {

    if (RTYPE == INTSXP) {
        CnstrtInt = MakeConstraints(compVec, mainFun, myReps,
                                    tarIntVals, ctype, n, m, IsComb,
                                    KeepRes, IsMult, IsRep);
        CnstrtInt->Prepare(compVec.front(), vInt);
    } else {
        CnstrtDbl = MakeConstraints(compVec, mainFun, myReps, tarVals, ctype,
                                    n, m, IsComb, KeepRes, IsMult, IsRep);
        CnstrtDbl->Prepare(compVec.front(), vNum);
    }

    std::vector<double> dblVec;
    double vecMax = std::floor(dblVec.max_size() / m);
    upperBoundDbl = std::min(vecMax, dblIntMax);

    std::vector<int> intVec;
    vecMax = std::floor(intVec.max_size() / m);
    upperBoundInt = std::min(vecMax, dblIntMax);
}

void CnstrntsToR::startOver() {
    keepGoing = true;

    if (RTYPE == INTSXP) {
        tarIntVals = origTarIntVals;
        CnstrtInt->Reset();
        CnstrtInt->Prepare(compVec.front(), vInt);
    } else {
        tarVals = origTarVals;
        CnstrtDbl->Reset();
        CnstrtDbl->Prepare(compVec.front(), vNum);
    }
}

SEXP CnstrntsToR::nextComb() {

    if (keepGoing) {
        return GetNext();
    } else {
        return Rf_ScalarLogical(false);
    }
}

SEXP CnstrntsToR::nextNumCombs(SEXP RNum) {

    int num;
    CleanConvert::convertPrimitive(RNum, num, VecType::Integer,
                                   "The number of results");

    if (keepGoing) {
        return GetNextN(num);
    } else {
        return Rf_ScalarLogical(false);
    }
}

SEXP CnstrntsToR::nextGather() {

    if (keepGoing) {
        const int num = (RTYPE == INTSXP) ? maxRows - CnstrtInt->GetCount() :
                                            maxRows - CnstrtDbl->GetCount();
        return GetNextN(num);
    } else {
        return Rf_ScalarLogical(false);
    }
}

SEXP CnstrntsToR::currComb() {

    if (!keepGoing) {
        return Rf_ScalarLogical(0);
    } else if (RTYPE == INTSXP && CnstrtInt->GetCount()) {
        return CnstrtVecReturn<INTSXP>(currIntVec);
    } else if (RTYPE == REALSXP && CnstrtDbl->GetCount()) {
        return CnstrtVecReturn<REALSXP>(currDblVec);
    } else {
        const std::string message = "Iterator Initialized. To see the first "
                                    "result, use the nextIter method(s)\n\n";
        Rprintf(message.c_str());
        return Rf_ScalarLogical(0);
    }
}

SEXP CnstrntsToR::summary() {
    SEXP parent = PROTECT(Combo::summary());
    std::string desc(R_CHAR(STRING_ELT(VECTOR_ELT(parent, 0), 0)));

    std::string val1 = (RTYPE == INTSXP) ?
                       std::to_string(tarIntVals.front()) :
                       std::to_string(tarVals.front());

    desc += " where the " + mainFun + " is ";

    if (origTarVals.size() == 2) {
        std::string val2 = (RTYPE == INTSXP) ?
                            std::to_string(tarIntVals.back()) :
                            std::to_string(tarVals.back());

        if (compVec.size() == 1) {
            desc += "between (" + compVec.front() +
                    ") " + val1 + " and " + val2;
        } else {
            desc += compVec.front() + " " + val1 +
                    " or " + compVec.back() + " " + val2;
        }
    } else {
        desc += compVec.front() + " " + val1;
    }

    const int idx = (RTYPE == INTSXP) ? CnstrtInt->GetCount() :
                                        CnstrtDbl->GetCount();
    const char *names[] = {"description", "currentIndex",
                           "totalResultsWithoutConstraints", ""};

    SEXP res = PROTECT(Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(res, 0, Rf_mkString(desc.c_str()));
    SET_VECTOR_ELT(res, 1, Rf_ScalarInteger(idx));
    SET_VECTOR_ELT(res, 2, VECTOR_ELT(parent, 2));

    UNPROTECT(2);
    return res;
}
