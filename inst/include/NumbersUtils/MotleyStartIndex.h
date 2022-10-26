#pragma once

#include <cstdint>

namespace MotleyPrimes {
    // This function is slightly different than the getStartingIndex
    // in the DivisorsContainer.cpp file. The step passed in this
    // function is a power of prime and requires an additional check
    // (i.e. else if (myPrime < lowerB)).
    template <typename T>
    inline T getStartIndexPowP(T lowerB, T step, T myPrime) {

        T remTest = lowerB % step;

        if (remTest == 0) {
            return 0;
        } else if (myPrime < lowerB) {
            return (step - remTest);
        } else {
            return (step - lowerB);
        }
    }
}
