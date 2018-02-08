#include <Rcpp.h>
#include <math.h>
#include <stdint.h>
#include "PollardRho.h"
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
    if (mTest <= 0) {stop("n must be positive");}
    m = (int)ceil(mTest);
    
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
    if (mTest <= 0) {stop("n must be positive");}
    m = (int)ceil(mTest);
    
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

template <typename typeRcpp, typename typeStd>
typeRcpp Factorize (typeStd t, std::vector<int64_t>& factors) {
    
    if (t == 1) {
        std::vector<typeStd> trivialReturn;
        if (factors.size() > 0) {trivialReturn.push_back(factors[0]);}
        trivialReturn.push_back(1);
        return wrap(trivialReturn);
    } else {
        std::vector<unsigned long int> lengths;
        std::vector<int64_t>::iterator it, facEnd;
        facEnd = factors.end();
        int64_t prev = factors[0];
        
        unsigned long int i, j, k, n = factors.size(), numUni = 0;
        std::vector<typeStd> uniFacs(n);
        uniFacs[0] = factors[0];
        lengths.reserve(n);
        lengths.push_back(1);
        
        for(it = factors.begin() + 1; it < facEnd; it++) {
            if (prev == *it) {
                lengths[numUni]++;
            } else {
                numUni++;
                prev = *it;
                lengths.push_back(1);
                uniFacs[numUni] = (typeStd)*it;
            }
        }
        
        unsigned long int ind, facSize = 1, numFacs = 1;
        for (i = 0; i <= numUni; i++) {numFacs *= (lengths[i]+1);}
        
        std::vector<typeStd> myFacs(numFacs);
        typeStd temp;
        
        for (i = 0; i <= lengths[0]; ++i) {
            myFacs[i] = (typeStd)std::pow(uniFacs[0], i);
        }
        
        if (numUni > 0) {
            for (j = 1; j <= numUni; j++) {
                facSize *= (lengths[j-1] + 1);
                for (i = 1; i <= lengths[j]; i++) {
                    ind = i*facSize;
                    for (k = 0; k < facSize; k++) {
                        temp = (typeStd)std::pow(uniFacs[j], i);
                        temp *= myFacs[k];
                        myFacs[ind + k] = temp;
                    }
                }
            }
        }
        
        std::sort(myFacs.begin(), myFacs.end());
        return wrap(myFacs);
    }
}


// [[Rcpp::export]]
SEXP getAllDivisorsRcpp (SEXP n) {
    double m;
    int64_t mPass;
    std::vector<int64_t> factors;
    bool isNegative = false;
    
    switch(TYPEOF(n)) {
        case REALSXP: {
            m = as<double>(n);
            break;
        }
        case INTSXP: {
            m = as<double>(n);
            break;
        }
        default: {
            stop("n must be of type numeric or integer");
        }
    }
    
    if (m < 0) {
        m = std::abs(m);
        isNegative = true;
    }
    
    mPass = (int64_t) m;
    
    if (m > 0) {
        getPrimefactors(mPass, factors);
        if (m <= INT_MAX) {
            int mInt = (int) m;
            IntegerVector myIntDivisors = Factorize<IntegerVector>(mInt, factors);
            
            if (isNegative) {
                unsigned int facSize = myIntDivisors.size();
                std::vector<int> tempInt(2 * facSize);
                unsigned int posInd = facSize, negInd = facSize - 1;
                
                for (std::size_t i = 0; i < facSize; i++, posInd++, negInd--) {
                    tempInt[negInd] = -1 * myIntDivisors[i];
                    tempInt[posInd] = myIntDivisors[i];
                }
                
                return wrap(tempInt);
            } else {
                return myIntDivisors;
            }
        } else {
            NumericVector myDblDivisors = Factorize<NumericVector>(m, factors);
            
            if (isNegative) {
                unsigned int facSize = myDblDivisors.size();
                std::vector<double> tempDbl(2 * facSize);
                unsigned int posInd = facSize, negInd = facSize - 1;
                
                for (std::size_t i = 0; i < facSize; i++, posInd++, negInd--) {
                    tempDbl[negInd] = -1 * myDblDivisors[i];
                    tempDbl[posInd] = myDblDivisors[i];
                }
                
                return wrap(tempDbl);
            } else {
                return myDblDivisors;
            }
        }
    } else {
        return wrap(factors);
    }
}