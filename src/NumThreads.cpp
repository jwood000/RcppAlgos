#include <thread>
#include "NumThreads.h"

SEXP cpp11GetNumThreads() {
    const int nThreads =  std::thread::hardware_concurrency();
    SEXP sexpRes = Rf_ScalarInteger(nThreads);
    return sexpRes;
}
