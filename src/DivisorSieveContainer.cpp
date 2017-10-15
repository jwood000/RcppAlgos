#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector NumDivisorsSieve (SEXP n) {
    int m;
    
    switch(TYPEOF(n)) {
        case REALSXP: {
            m = as<int>(n);
            break;
        }
        case INTSXP: {
            m = as<int>(n);
            break;
        }
        default: {
            stop("n must be of type numeric or integer");
        }
    }
    
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
    
    switch(TYPEOF(n)) {
        case REALSXP: {
            m = as<int>(n);
            break;
        }
        case INTSXP: {
            m = as<int>(n);
            break;
        }
        default: {
            stop("n must be of type numeric or integer");
        }
    }
    
    std::vector<std::vector<int> > myDivList(m, std::vector<int>(1, 1));
    int i, j;
    
    for (i = 2; i <= m; i++) {
        for (j = i; j <= m; j+=i) {
            myDivList[j - 1].push_back(i);
        }
    }
    
    return wrap(myDivList);
}