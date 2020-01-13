#ifndef MOTLEY_START_INDEX_H
#define MOTLEY_START_INDEX_H

#include <Rcpp.h>

namespace MotleyPrimes {
    // This function is slightly different than the getStartingIndex
    // in the DivisorsContainer.cpp file. The step passed in this
    // function is a power of prime and requires an additional check
    // (i.e. else if (myPrime < lowerB)).
    template <typename typeInt>
    inline typeInt getStartIndexPowP(typeInt lowerB, typeInt step, typeInt myPrime) {
        
        typeInt retStrt;
        typeInt remTest = lowerB % step;
        
        if (remTest == 0) {
            retStrt = 0;
        } else if (myPrime < lowerB) {
            retStrt = step - remTest;
        } else {
            retStrt = step - lowerB;
        }
        
        return retStrt;
    }
}

#endif
