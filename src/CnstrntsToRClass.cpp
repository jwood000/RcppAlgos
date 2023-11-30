#include "Constraints/CnstrntsToRClass.h"
#include <sstream>
#include <iomanip>
#include <limits>

template <typename T>
void SetCurrVec(const std::vector<T> &cnstrntVec,
                const std::vector<T> &resVec,
                std::vector<T> &curr, std::size_t m, bool Keep) {

    std::vector<T> newCurr(cnstrntVec.end() - std::min(m, cnstrntVec.size()),
                           cnstrntVec.end());
    if (Keep) newCurr.push_back(resVec.back());
    curr = newCurr;
}

template <int sexpType, typename T>
SEXP CnstrtVecReturn(const std::vector<T> &v) {

    cpp11::sexp res = Rf_allocVector(sexpType, v.size());

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
        Rprintf("%s", message.c_str());
        return R_NilValue;
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

SEXP CnstrntsToR::GetNextN(int n) {

    if (RTYPE == INTSXP) {
        std::vector<int> resVec;
        std::vector<int> cnstrntVec;

        GetNSolutions(
            compVec, CnstrtInt, cnstrntVec, resVec, vInt, tarIntVals, n
        );

        if (cnstrntVec.size()) {
            SetCurrVec(cnstrntVec, resVec, currIntVec, width, KeepRes);
            const std::size_t vecLen = cnstrntVec.size();
            const std::size_t numResult = vecLen / m;

            cpp11::sexp res = Rf_allocMatrix(INTSXP, numResult, nCols);
            int* matInt = INTEGER(res);

            VectorToMatrix(cnstrntVec, resVec, matInt, 0, numResult,
                           width, upperBoundInt, KeepRes, false);
            return res;
        }
    } else {
        std::vector<double> resVec;
        std::vector<double> cnstrntVec;
        GetNSolutions(compVec, CnstrtDbl, cnstrntVec, resVec, vNum, tarVals, n);

        if (cnstrntVec.size()) {
            SetCurrVec(cnstrntVec, resVec, currDblVec, width, KeepRes);
            const std::size_t vecLen = cnstrntVec.size();
            const std::size_t numResult = vecLen / m;

            cpp11::sexp res = Rf_allocMatrix(REALSXP, numResult, nCols);
            double* matNum = REAL(res);

            VectorToMatrix(cnstrntVec, resVec, matNum, 0, numResult,
                           width, upperBoundDbl, KeepRes, false);
            return res;
        }
    }

    keepGoing = false;
    const std::string message = "No more results.\n\n";
    Rprintf("%s", message.c_str());
    return R_NilValue;
}

CnstrntsToR::CnstrntsToR(
    SEXP Rv, int Rm, SEXP RcompRows, const std::vector<int> &bVec,
    const std::vector<int> &Rreps, const std::vector<int> &Rfreqs,
    const std::vector<int> &RvInt, const std::vector<double> &RvNum,
    VecType typePass, int RmaxThreads, SEXP RnumThreads, bool Rparallel,
    const PartDesign &Rpart, const std::vector<std::string> &RcompVec,
    std::vector<double> &RtarVals, std::vector<int> &RtarIntVals,
    std::vector<int> &RstartZ, const std::string &RmainFun,
    const std::string &RFunTest, funcPtr<double> RfunDbl,
    ConstraintType Rctype, int RstrtLen, int Rcap, bool RKeepRes,
    bool RnumUnknown, double RcnstrtRows, const mpz_class &RcnstrtRowsMpz
) : ComboRes(Rv, Rm, RcompRows, bVec, Rreps, Rfreqs, RvInt, RvNum, typePass,
             RmaxThreads, RnumThreads, Rparallel, Rpart, RcompVec, RtarVals,
             RtarIntVals, RstartZ, RmainFun, RFunTest, RfunDbl, Rctype,
             RstrtLen, Rcap, RKeepRes, RnumUnknown, RcnstrtRows,
             RcnstrtRowsMpz),
    maxRows(std::min(dblIntMax, RcnstrtRows)),
    origTarIntVals(RtarIntVals), origTarVals(RtarVals) {

    if (RTYPE == INTSXP) {
        CnstrtInt = MakeConstraints(compVec, mainFun, funTest, myReps,
                                    tarIntVals, ctype, n, m, IsComb,
                                    KeepRes, IsMult, IsRep);
        CnstrtInt->Prepare(compVec.front(), vInt);
    } else {
        CnstrtDbl = MakeConstraints(compVec, mainFun, funTest, myReps,
                                    tarVals, ctype, n, m, IsComb,
                                    KeepRes, IsMult, IsRep);
        CnstrtDbl->Prepare(compVec.front(), vNum);
    }

    std::vector<double> dblVec;
    double vecMax = std::floor(dblVec.max_size() / m);
    upperBoundDbl = std::min(vecMax, dblIntMax);

    std::vector<int> intVec;
    vecMax = std::floor(intVec.max_size() / m);
    upperBoundInt = std::min(vecMax, dblIntMax);
    prevIterAvailable = false;
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
        return R_NilValue;
    }
}

SEXP CnstrntsToR::nextNumCombs(SEXP RNum) {

    int num;
    CppConvert::convertPrimitive(RNum, num, VecType::Integer,
                                   "The number of results");

    if (keepGoing) {
        return GetNextN(num);
    } else {
        return R_NilValue;
    }
}

SEXP CnstrntsToR::nextGather() {

    if (keepGoing) {
        const int num = (RTYPE == INTSXP) ? maxRows - CnstrtInt->GetCount() :
                                            maxRows - CnstrtDbl->GetCount();
        return GetNextN(num);
    } else {
        return R_NilValue;
    }
}

SEXP CnstrntsToR::currComb() {

    if (!keepGoing) {
        return R_NilValue;
    } else if (RTYPE == INTSXP && CnstrtInt->GetCount()) {
        return CnstrtVecReturn<INTSXP>(currIntVec);
    } else if (RTYPE == REALSXP && CnstrtDbl->GetCount()) {
        return CnstrtVecReturn<REALSXP>(currDblVec);
    } else {
        return ToSeeFirst(false);
    }
}

SEXP CnstrntsToR::summary() {
    cpp11::sexp res = Combo::summary();
    std::string desc(R_CHAR(STRING_ELT(VECTOR_ELT(res, 0), 0)));

    constexpr auto max_digits10 = std::numeric_limits<double>::max_digits10;
    const double testVal1 = (funTest == "mean") ?
                            tarVals.front() / static_cast<double>(m) :
                            tarVals.front();

    const double testVal2 = (funTest == "mean") ?
                            tarVals.back() / static_cast<double>(m) :
                            tarVals.back();

    const double dblVal1 = (RTYPE == INTSXP) ? tarIntVals.front() :
                           testVal1;

    std::stringstream ss1;
    ss1 << std::setprecision(max_digits10) << dblVal1;
    std::string val1;
    ss1 >> val1;

    desc += " where the " + funTest + " is ";

    if (origTarVals.size() == 2) {
        const double dblVal2 = (RTYPE == INTSXP) ? tarIntVals.back() :
                               testVal2;

        std::stringstream ss2;
        ss2 << std::setprecision(max_digits10) << dblVal2;
        std::string val2;
        ss2 >> val2;

        const bool is_equal = (RTYPE == INTSXP) ?
                              tarIntVals.front() == tarIntVals.back() :
                              tarVals.front() == tarVals.back();

        if (is_equal) {
            desc += "equal to " + val1;
        } else if (compVec.size() == 1) {
            // N.B. The smallest value corresponds to back(). That is why
            // val2 comes before val1. From UserConstraintFuns.cpp:
            //
            // template <typename T>
            // bool greaterEqlLessEql(T x, const std::vector<T> &y) {return x <= y[0] && x >= y[1];}
            //
            // Here we see that the first element is the largest
            desc += "between " + val2 + " and " + val1;
        } else {
            desc += compVec.front() + " " + val1 +
                    " or " + compVec.back() + " " + val2;
        }
    } else {
        desc += compVec.front() + " " + val1;
    }

    const int idx = (RTYPE == INTSXP) ? CnstrtInt->GetCount() :
                                        CnstrtDbl->GetCount();

    SET_VECTOR_ELT(res, 0, Rf_mkString(desc.c_str()));
    SET_VECTOR_ELT(res, 1, Rf_ScalarInteger(idx));
    SET_VECTOR_ELT(res, 2, Rf_ScalarReal(R_NaReal));
    SET_VECTOR_ELT(res, 3, Rf_ScalarReal(R_NaReal));
    return res;
}
