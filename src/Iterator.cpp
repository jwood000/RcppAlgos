#include "ClassUtils/Iterator.h"

Iterator::Iterator(SEXP Rv, SEXP RcompRow, int RmaxThreads,
                   SEXP RnThreads, bool Rparallel, bool IsGmp) :
    sexpVec(Rv), n(Rf_length(Rv)), RTYPE(TYPEOF(Rv)),
    maxThreads(RmaxThreads), sexpNThreads(RnThreads), Parallel(Rparallel),
    IsGmp(IsGmp), computedRows(IsGmp ? 0 : Rf_asReal(RcompRow)) {

    if (IsGmp) {
        CppConvert::convertMpzClass(RcompRow, computedRowsMpz,
                                    "computedRowsMpz");
    }
}

SEXP Iterator::sourceVector() const {
    return sexpVec;
}
