#include "ClassUtils/Iterator.h"

SEXP Iterator::ToSeeLast(bool AdjustIdx) {

    std::string message = "No more results.";

    if (prevIterAvailable) {
        message += " To see the last result, use the prevIter method(s)\n\n";
    } else {
        message += "\n\n";
    }

    Rprintf("%s", message.c_str());
    if (AdjustIdx) increment(IsGmp, mpzIndex, dblIndex);
    return R_NilValue;
}

SEXP Iterator::ToSeeFirst(bool AdjustIdx) {

    const std::string message = "Iterator Initialized. To see the first"
        " result, use the nextIter method(s)\n\n";

    Rprintf("%s", message.c_str());
    if (AdjustIdx) decrement(IsGmp, mpzIndex, dblIndex);
    return R_NilValue;
}

Iterator::Iterator(SEXP Rv, VecType typePass, SEXP RcompRow, int RmaxThreads,
                   SEXP RnThreads, bool Rparallel, bool IsGmp) :
    n(Rf_length(Rv)), sexpVec(Rv), RTYPE(TYPEOF(Rv)), myType(typePass),
    maxThreads(RmaxThreads), sexpNThreads(RnThreads), Parallel(Rparallel),
    IsGmp(IsGmp), computedRows(IsGmp ? 0 : Rf_asReal(RcompRow)) {

    if (IsGmp) {
        CppConvert::convertMpzClass(RcompRow, computedRowsMpz,
                                    "computedRowsMpz");
    }

    prevIterAvailable = true;
}

SEXP Iterator::sourceVector() const {
    return sexpVec;
}
