#include "Constraints/CnstrntsSpecialClass.h"

CnstrntsSpecial::CnstrntsSpecial(
    SEXP Rv, int Rm, SEXP RcompRows, const std::vector<int> &bVec,
    const std::vector<int> &Rreps, const std::vector<int> &Rfreqs,
    const std::vector<int> &RvInt, const std::vector<double> &RvNum,
    VecType typePass, int RmaxThreads, SEXP RnumThreads, bool Rparallel,
    const PartDesign &Rpart, const std::vector<std::string> &RcompVec,
    std::vector<double> &RtarVals, std::vector<int> &RtarIntVals,
    std::vector<int> &RstartZ, const std::string &RmainFun,
    const std::string &RFunTest, funcPtr<double> RfunDbl,
    ConstraintType Rctype, int RstrtLen, int Rcap, bool RKeepRes,
    bool RnumUnknown, double RcnstrtRows, mpz_t RcnstrtRowsMpz
) : ComboRes(Rv, Rm, RcompRows, bVec, Rreps, Rfreqs, RvInt, RvNum, typePass,
             RmaxThreads, RnumThreads, Rparallel, Rpart, RcompVec, RtarVals,
             RtarIntVals, RstartZ, RmainFun, RFunTest, RfunDbl, Rctype,
             RstrtLen, Rcap, RKeepRes, RnumUnknown, RcnstrtRows,
             RcnstrtRowsMpz) {
    count = 0;
    keepGoing = true;
    prevIterAvailable = false;
}

void CnstrntsSpecial::startOver() {
    count = 0;
    keepGoing = true;
    Combo::startOver();
}

SEXP CnstrntsSpecial::nextComb() {
    if (keepGoing) {
        SEXP res = PROTECT(ComboRes::nextNumCombs(Rf_ScalarInteger(1)));

        if (Rf_isNull(res)) {
            keepGoing = false;
            UNPROTECT(1);
            return res;
        } else {
            if (Rf_nrows(res)) {
                count = dblIndex;
                Rf_setAttrib(res, R_DimSymbol, R_NilValue);
                UNPROTECT(1);
                return res;
            } else {
                UNPROTECT(1);
                keepGoing = false;
                return ToSeeLast();
            }
        }
    } else {
        keepGoing = false;
        return R_NilValue;
    }
}

SEXP CnstrntsSpecial::nextNumCombs(SEXP RNum) {

    if (keepGoing) {
        SEXP res = PROTECT(ComboRes::nextNumCombs(RNum));

        if (Rf_isNull(res)) {
            keepGoing = false;
            UNPROTECT(1);
            return res;
        } else {
            int num;
            CleanConvert::convertPrimitive(RNum, num, VecType::Integer,
                                           "The number of results");

            if (Rf_nrows(res)) {
                const int returned_nrows = Rf_nrows(res);
                keepGoing = num == returned_nrows;
                count = dblIndex - (num - returned_nrows);
                UNPROTECT(1);
                return res;
            } else {
                UNPROTECT(1);
                keepGoing = false;
                return ToSeeLast();
            }
        }
    } else {
        keepGoing = false;
        return R_NilValue;
    }
}

SEXP CnstrntsSpecial::nextGather() {

    if (keepGoing) {
        SEXP res = PROTECT(ComboRes::nextGather());

        if (Rf_isNull(res)) {
            keepGoing = false;
            UNPROTECT(1);
            return res;
        } else if (Rf_nrows(res)) {
            count += Rf_nrows(res);
            keepGoing = false;
            UNPROTECT(1);
            return res;
        } else {
            UNPROTECT(1);
            keepGoing = false;
            return ToSeeLast();
        }
    } else {
        keepGoing = false;
        return R_NilValue;
    }
}

SEXP CnstrntsSpecial::currComb() {
    return ComboRes::currComb();
}

SEXP CnstrntsSpecial::summary() {
    SEXP res = PROTECT(Combo::summary());
    std::string desc(R_CHAR(STRING_ELT(VECTOR_ELT(res, 0), 0)));

    std::string val1 = (RTYPE == INTSXP) ?
                       std::to_string(tarIntVals.front()) :
                       std::to_string(tarVals.front());

    desc += " where the " + mainFun + " is ";

    if (tarVals.size() == 2) {
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

    SET_VECTOR_ELT(res, 0, Rf_mkString(desc.c_str()));
    SET_VECTOR_ELT(res, 1, Rf_ScalarInteger(count));
    SET_VECTOR_ELT(res, 2, Rf_ScalarReal(R_NaReal));
    SET_VECTOR_ELT(res, 3, Rf_ScalarReal(R_NaReal));

    UNPROTECT(1);
    return res;
}
