#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector NumDivisorsSieve (SEXP n) {
    int m;
    double mTest;
    
    switch(TYPEOF(n)) {
        case REALSXP: {
            mTest = as<double>(n);
            break;
        }
        case INTSXP: {
            mTest = as<double>(n);
            break;
        }
        default: {
            stop("n must be of type numeric or integer");
        }
    }
    
    if (mTest > 2147483647) {stop("n must be less than 2^31");}
    if (mTest < 0) {stop("n must be positive");}
    m = mTest;
    
    std::vector<int> numFacs(m, 1);
    int i, j;
    
    for (i = 2; i <= m; i++) {
        for (j = i; j <= m; j+=i) {numFacs[j - 1]++;}
    }
    
    return wrap(numFacs);
}

// [[Rcpp::export]]
List DivisorListRcpp (SEXP n) {
    int m;
    double mTest;
    
    switch(TYPEOF(n)) {
        case REALSXP: {
            mTest = as<double>(n);
            break;
        }
        case INTSXP: {
            mTest = as<double>(n);
            break;
        }
        default: {
            stop("n must be of type numeric or integer");
        }
    }
    
    if (mTest > 2147483647) {stop("n must be less than 2^31");}
    if (mTest < 0) {stop("n must be positive");}
    m = mTest;
    
    std::vector<std::vector<int> > myDivList(m, std::vector<int>(1, 1));
    std::vector<std::vector<int> >::iterator it2d, itEnd;
    itEnd = myDivList.end();
    // Most values will have fewer than 2 times the 
    // maximal number of bits in m (crude analysis).
    // We don't want to consider the few highly
    // composite values when determing our memory
    // reservation as that will allocate way more
    // memory than necessary for the majority of values.
    int i, j, myMalloc = 2*ceil(log2(m));
    
    for (it2d = myDivList.begin(); it2d < itEnd; it2d++) {
        it2d -> reserve(myMalloc);
    }
    
    for (i = 2; i <= m; i++) {
        for (j = i; j <= m; j+=i) {
            myDivList[j - 1].push_back(i);
        }
    }
    
    return wrap(myDivList);
}