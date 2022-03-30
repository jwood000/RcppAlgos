#include "cpp11/R.hpp"
#include "cpp11/sexp.hpp"
#include <thread>

[[cpp11::register]]
SEXP cpp11GetNumThreads() {
    const int nThreads =  std::thread::hardware_concurrency();
    return cpp11::as_sexp(nThreads);
}
