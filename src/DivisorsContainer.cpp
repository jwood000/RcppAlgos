#include <Rcpp.h>
#include <cmath>
#include <libdivide.h>
#include <PollardRho.h>

template <typename typeInt>
inline typeInt getStartingIndex (typeInt lowerB, typeInt step) {
    
    if (step >= lowerB)
        return (2 * step - lowerB);
        
    typeInt retStrt, remTest = lowerB % step;
    retStrt = (remTest == 0) ? 0 : (step - remTest);

    return retStrt;
}

template <typename typeInt>
Rcpp::IntegerVector NumDivisorsSieve (typeInt m, typeInt n, bool keepNames) {
    
    typeInt myRange = n;
    myRange += (1 - m);
    std::vector<int> numFacs;
    std::vector<double> myNames;
    
    if (keepNames){
        myNames.resize(myRange);
        double dblM = (double) m;
        for (std::size_t k = 0; dblM <= n; ++dblM, ++k)
            myNames[k] = dblM;
    }
    
    numFacs = std::vector<int> (myRange, 2);
    typeInt j, myNum, sqrtBound = (typeInt) std::sqrt((double) n);
    
    for (typeInt i = 2; i <= sqrtBound; ++i) {
        myNum = (i * sqrtBound) - m;
        j = getStartingIndex(m, i);
        
        for (; j <= myNum; j += i)
            ++numFacs[j];
        
        for (; j < myRange; j += i)
            numFacs[j] += 2;
    }
    
    // Subtract 1 from the first entry as 1 has only itself
    // as a divisor. N.B. myRange was initialized with 2
    if (m < 2)
        --numFacs[0];
    
    Rcpp::IntegerVector myVector = Rcpp::wrap(numFacs);
    if (keepNames)
        myVector.attr("names") = myNames;
    
    return myVector;
}

template <typename typeInt, typename typeReturn>
Rcpp::List DivisorsSieve (typeInt m, typeReturn retN, bool keepNames) {
    
    typeInt n = (typeInt) retN;
    typeInt myRange = n;
    myRange += (1 - m);
    
    std::vector<std::vector<typeReturn> > 
        MyDivList(myRange,std::vector<typeReturn>(1, 1));
    
    typename std::vector<std::vector<typeReturn> >::iterator it2d, itEnd;
    itEnd = MyDivList.end();
    
    typeInt i, j, myNum = m;
    std::vector<typeReturn> myNames;
    
    if (keepNames){
        myNames.resize(myRange);
        typeReturn retM = (typeReturn) m;
        for (std::size_t k = 0; retM <= retN; ++retM, ++k)
            myNames[k] = retM;
    }
    
    Rcpp::IntegerVector myMemory = NumDivisorsSieve(m, (typeInt) n, false);
    Rcpp::IntegerVector::iterator myMalloc;
    
    if (m < 2)
        myMalloc = myMemory.begin() + 1;
    else
        myMalloc = myMemory.begin();
    
    if (m < 2) {
        for (it2d = MyDivList.begin() + 1; it2d < itEnd; ++it2d, ++myMalloc)
            it2d->reserve(*myMalloc);
        
        for (i = 2; i <= n; ++i)
            for (j = i; j <= n; j += i)
                MyDivList[j - 1].push_back((typeReturn) i);
        
    } else {
        int_fast32_t sqrtBound = floor(sqrt((double)n));
        typeInt myStart, testNum;
        
        for (it2d = MyDivList.begin(); it2d < itEnd; ++it2d, ++myNum, ++myMalloc) {
            it2d->reserve(*myMalloc);
            it2d->push_back((typeReturn) myNum);
        }
        
        for (i = sqrtBound; i >= 2; --i) {
            
            myStart = getStartingIndex(m, i);
            typeInt myNum = m + myStart;
            libdivide::divider<typeInt> fastDiv(i);
            
            for (j = myStart; j < myRange; j += i, myNum += i) {
                // Put element in the second position. (see comment below)
                MyDivList[j].insert(MyDivList[j].begin() + 1, (typeReturn) i);
                testNum = myNum / fastDiv;
                
                // Ensure we won't duplicate adding an element. If
                // testNum <= sqrtBound, it will be added in later
                // iterations. Also, we insert this element in the
                // pentultimate position as it will be the second
                // to the largest element at the time of inclusion.
                // E.g. let i = 5, myNum = 100, so the current
                // vectors looks like so: v = 1, 5, 10, 100 (5 was
                // added to the second position above). With i = 5,
                // testNum = 100 / 5 = 20, thus we add it the second
                // to last position to give v = 1 5 10 20 100.
                if (testNum > sqrtBound)
                    MyDivList[j].insert(MyDivList[j].end() - 1, (typeReturn) testNum);
            }
        }
    }
    
    Rcpp::List myList = Rcpp::wrap(MyDivList);
    if (keepNames)
        myList.attr("names") = myNames;
    
    return myList;
}

// [[Rcpp::export]]
SEXP DivisorsGeneral (SEXP Rb1, SEXP Rb2, 
                       SEXP RIsList, SEXP RNamed) {
    double bound1, bound2, myMax, myMin;
    bool isList = false, isNamed = false;
    
    switch(TYPEOF(Rb1)) {
        case REALSXP: {
            bound1 = Rcpp::as<double>(Rb1);
            break;
        }
        case INTSXP: {
            bound1 = Rcpp::as<double>(Rb1);
            break;
        }
        default: {
            Rcpp::stop("bound1 must be of type numeric or integer");
        }
    }
    
    isList = Rcpp::as<bool>(RIsList);
    isNamed = Rcpp::as<bool>(RNamed);
    
    if (bound1 <= 0 || bound1 > Significand53)
        Rcpp::stop("bound1 must be a positive number less than 2^53");
    
    if (Rf_isNull(Rb2)) {
        myMax = floor(bound1);
        
        if (isList) {
            if (myMax < 2) {
                std::vector<std::vector<int> > trivialRet(1, std::vector<int>(1, 1));
                Rcpp::List z = Rcpp::wrap(trivialRet);
                if (isNamed)
                    z.attr("names") = 1;
                
                return z;
            }
            if (myMax > (INT_MAX - 1)) {
                return DivisorsSieve((int_fast64_t) 1, (double) myMax, isNamed);
            }
            return DivisorsSieve((int_fast32_t) 1,
                                          (int_fast32_t) myMax, isNamed);
        } else {
            if (myMax < 2) {
                Rcpp::IntegerVector v(1, 1);
                if (isNamed)
                    v.attr("names") = 1;
                
                return v;
            }
            if (myMax > (INT_MAX - 1))
                return NumDivisorsSieve((int_fast64_t) 1, (int_fast64_t) myMax, isNamed);
                
            return NumDivisorsSieve((int_fast32_t) 1, (int_fast32_t) myMax, isNamed);
        }
    } else {
        switch(TYPEOF(Rb2)) {
            case REALSXP: {
                bound2 = Rcpp::as<double>(Rb2);
                break;
            }
            case INTSXP: {
                bound2 = Rcpp::as<double>(Rb2);
                break;
            }
            default: {
                Rcpp::stop("bound2 must be of type numeric or integer");
            }
        }
        if (bound2 <= 0 || bound2 > Significand53)
            Rcpp::stop("bound2 must be a positive number less than 2^53");
        
        if (bound1 > bound2) {
            myMax = bound1;
            myMin = bound2;
        } else {
            myMax = bound2;
            myMin = bound1;
        }
        
        myMin = ceil(myMin);
        myMax = floor(myMax);
        
        if (isList) {
            if (myMax < 2) {
                std::vector<std::vector<int> > trivialRet(1, std::vector<int>(1, 1));
                Rcpp::List z = Rcpp::wrap(trivialRet);
                if (isNamed)
                    z.attr("names") = 1;
                
                return z;
            }
            if (myMax > (INT_MAX - 1))
                return DivisorsSieve((int64_t) myMin, (double) myMax, isNamed);
            
            return DivisorsSieve((int32_t) myMin, (int32_t) myMax, isNamed);
        } else {
            if (myMax < 2) {
                Rcpp::IntegerVector v(1, 1);
                if (isNamed)
                    v.attr("names") = 1;
                
                return v;
            }
            if (myMax > (INT_MAX - 1))
                return NumDivisorsSieve((int64_t) myMin, (int64_t) myMax, isNamed);
                
            return NumDivisorsSieve((int32_t) myMin, (int32_t) myMax, isNamed);
        }
    }
}

template <typename typeReturn>
std::vector<typeReturn> Factorize (std::vector<typeReturn>& factors) {
    
    unsigned long int n = factors.size();
    
    if (n == 1) {
        std::vector<typeReturn> primeReturn;
        primeReturn.push_back(1);
        primeReturn.push_back(factors[0]);
        return primeReturn;
    } else {
        std::vector<unsigned long int> lengths;
        typename std::vector<typeReturn>::iterator it, facEnd;
        facEnd = factors.end();
        typeReturn prev = factors[0];
        
        unsigned long int i, j, k, numUni = 0;
        std::vector<typeReturn> uniFacs(n);
        uniFacs[0] = factors[0];
        lengths.reserve(n);
        lengths.push_back(1);
        
        for(it = factors.begin() + 1; it < facEnd; ++it) {
            if (prev == *it) {
                ++lengths[numUni];
            } else {
                ++numUni;
                prev = *it;
                lengths.push_back(1);
                uniFacs[numUni] = *it;
            }
        }
        
        unsigned long int ind, facSize = 1, numFacs = 1;
        for (i = 0; i <= numUni; ++i)
            numFacs *= (lengths[i]+1);
        
        std::vector<typeReturn> myFacs(numFacs);
        typeReturn temp;
        
        for (i = 0; i <= lengths[0]; ++i)
            myFacs[i] = std::pow(uniFacs[0], i);
        
        if (numUni > 0) {
            for (j = 1; j <= numUni; ++j) {
                facSize *= (lengths[j-1] + 1);
                for (i = 1; i <= lengths[j]; ++i) {
                    ind = i*facSize;
                    for (k = 0; k < facSize; ++k) {
                        temp = std::pow(uniFacs[j], i);
                        temp *= myFacs[k];
                        myFacs[ind + k] = temp;
                    }
                }
            }
        }
        
        std::sort(myFacs.begin(), myFacs.end());
        return myFacs;
    }
}

template <typename typeReturn>
Rcpp::List FactorList (std::vector<double> myNums, bool namedList) {
    
    int64_t mPass;
    bool isNegative = false;
    unsigned int myLen = myNums.size();
    
    std::vector<std::vector<typeReturn> > 
            MyDivList(myLen, std::vector<typeReturn>());
    
    for (std::size_t j = 0; j < myLen; ++j) {
        std::vector<typeReturn> myDivisors;
        mPass = (int64_t) myNums[j];
        
        if (mPass < 0) {
            mPass = std::abs(mPass);
            isNegative = true;
        } else {
            isNegative = false;
        }
        
        if (mPass > 1) {
            std::vector<typeReturn> factors;
            getPrimefactors(mPass, factors);
            myDivisors = Factorize<typeReturn>(factors);
            if (isNegative) {
                unsigned int facSize = myDivisors.size();
                std::vector<typeReturn> tempInt(2 * facSize);
                unsigned int posInd = facSize, negInd = facSize - 1;
                
                for (std::size_t i = 0; i < facSize; ++i, ++posInd, --negInd) {
                    tempInt[negInd] = -1 * myDivisors[i];
                    tempInt[posInd] = myDivisors[i];
                }
                
                myDivisors = tempInt;
            }
        } else {
            if (isNegative)
                myDivisors.push_back(-1);
            if (mPass > 0)
                myDivisors.push_back(1);
        }
        
        MyDivList[j] = myDivisors;
    }
    
    Rcpp::List myList = Rcpp::wrap(MyDivList);
    if (namedList)
        myList.attr("names") = myNums;
    
    return myList;
}


// [[Rcpp::export]]
SEXP getAllDivisorsRcpp (SEXP Rv, SEXP RNamed) {
    std::vector<double> myNums;
    bool isNamed = Rcpp::as<bool>(RNamed);
    
    switch(TYPEOF(Rv)) {
        case REALSXP: {
            myNums = Rcpp::as<std::vector<double> >(Rv);
            break;
        }
        case INTSXP: {
            myNums = Rcpp::as<std::vector<double> >(Rv);
            break;
        }
        default: {
            Rcpp::stop("v must be of type numeric or integer");
        }
    }
    
    if (myNums.size() > 1) {
        double myMax = *std::max_element(myNums.begin(), myNums.end());
        double myMin = *std::min_element(myNums.begin(), myNums.end());
        
        if (std::abs(myMin) > myMax)
            myMax = std::abs(myMin);
        
        if (myMax > Significand53)
            Rcpp::stop("the abs value of each element must be less than 2^53");
        
        if (myMax > INT_MAX)
            return FactorList<double>(myNums, isNamed);
        else
            return FactorList<int>(myNums, isNamed);
    } else {
        int64_t mPass = (int64_t) myNums[0];
        bool isNegative = false;
        
        if (mPass < 0) {
            mPass = std::abs(mPass);
            isNegative = true;
        } else {
            isNegative = false;
        }
        
        if (mPass > Significand53)
            Rcpp::stop("the abs value of each element must be less than 2^53");
        
        if (mPass > INT_MAX) {
            std::vector<double> myDivisors;
            std::vector<double> factors;
            getPrimefactors(mPass, factors);
            myDivisors = Factorize<double>(factors);
            
            if (isNegative) {
                unsigned int facSize = myDivisors.size();
                std::vector<double> tempInt(2 * facSize);
                unsigned int posInd = facSize, negInd = facSize - 1;
                
                for (std::size_t i = 0; i < facSize; ++i, ++posInd, --negInd) {
                    tempInt[negInd] = -1 * myDivisors[i];
                    tempInt[posInd] = myDivisors[i];
                }
                
                return Rcpp::wrap(tempInt);
            }
            
            return Rcpp::wrap(myDivisors);
            
        } else if (mPass > 1) {
            std::vector<int> myDivisors;
            std::vector<int> factors;
            getPrimefactors(mPass, factors);
            myDivisors = Factorize<int>(factors);
            
            if (isNegative) {
                unsigned int facSize = myDivisors.size();
                std::vector<int> tempInt(2 * facSize);
                unsigned int posInd = facSize, negInd = facSize - 1;
                
                for (std::size_t i = 0; i < facSize; ++i, ++posInd, --negInd) {
                    tempInt[negInd] = -1 * myDivisors[i];
                    tempInt[posInd] = myDivisors[i];
                }
                
                return Rcpp::wrap(tempInt);
            }
            
            return Rcpp::wrap(myDivisors);
            
        } else {
            std::vector<int> myDivisors;
            if (isNegative)
                myDivisors.push_back(-1);
            if (mPass > 0)
                myDivisors.push_back(1);
            
            return Rcpp::wrap(myDivisors);
        }
    }
}
