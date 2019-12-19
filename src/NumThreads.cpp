#include <thread>

// [[Rcpp::export]]
int cpp11GetNumThreads() {
    return std::thread::hardware_concurrency();
}
