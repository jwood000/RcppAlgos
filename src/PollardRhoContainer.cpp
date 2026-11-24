#include "CppConvert.h"

#include "NumbersUtils/PollardRhoDepends.h"
#include "NumbersUtils/PollardRhoUtils.h"
#include "ThreadBatch.h"
#include <cmath>

template <typename T>
void PollardRhoMain(const std::vector<double> &myNums, T myMax,
                    bool bPrimeFacs, bool bAllFacs,
                    std::vector<std::vector<T>> &MyList,
                    int* primeTest, std::size_t myRange,
                    int nThreads, int maxThreads) {

    bool Parallel = false;

    if (nThreads > 1 && myRange > 1 && maxThreads > 1) {
        Parallel = true;
        if (nThreads > maxThreads) {nThreads = maxThreads;}
        if ((myRange / nThreads) < 1) {nThreads = myRange;}
    }

    if (Parallel) {
        const std::size_t chunkSize = myRange / nThreads;
        std::size_t n = chunkSize - 1;
        std::size_t m = 0;
        RcppAlgos::ThreadBatch batch(nThreads);

        for (int j = 0; j < (nThreads - 1); m = n, n += chunkSize, ++j) {
            if (bPrimeFacs) {
                batch.launch([&, m, n] {
                    PrimeFacList<T>(m, n, myNums, MyList);
                });
            } else if (bAllFacs) {
                batch.launch([&, m, n] {
                    FactorList<T>(m, n, myNums, MyList);
                });
            } else {
                batch.launch([&, m, n] {
                    IsPrimeVec(m, n, myNums, primeTest);
                });
            }
        }

        if (bPrimeFacs) {
            batch.launch([&, m] {
                PrimeFacList<T>(m, myRange, myNums, MyList);
            });
        } else if (bAllFacs) {
            batch.launch([&, m] {
                FactorList<T>(m, myRange, myNums, MyList);
            });
        } else {
            batch.launch([&, m] {
                IsPrimeVec(m, myRange, myNums, primeTest);
            });
        }

        batch.finish();
    } else {
        if (bPrimeFacs) {
            PrimeFacList(0u, myRange, myNums, MyList);
        } else if (bAllFacs) {
            FactorList(0u, myRange, myNums, MyList);
        } else {
            IsPrimeVec(0u, myRange, myNums, primeTest);
        }
    }
}

template <typename T>
SEXP PolRhoGlue(const std::vector<double> &myNums, T myMax,
                bool bPrimeFacs, bool bAllFacs, bool keepNames,
                int nThreads, int maxThreads) {

    std::size_t myRange = myNums.size();
    int* tempVec = nullptr;

    if (bPrimeFacs) {
        if (myRange == 1) {
            std::int64_t mPass = static_cast<std::int64_t>(myNums[0]);
            if (mPass == 0) {return Rf_allocVector(INTSXP, 0);}
            std::vector<T> factors;

            if (mPass < 0) {
                mPass = std::abs(mPass);
                factors.push_back(-1);
            }

            GetPrimeFactors(mPass, factors);
            return cpp11::writable::r_vector<T>(factors);
        } else {
            std::vector<std::vector<T>>
                MyPrimeList(myRange, std::vector<T>());

            PollardRhoMain(myNums, myMax, bPrimeFacs, false, MyPrimeList,
                           tempVec, myRange, nThreads, maxThreads);

            cpp11::writable::list myList(myRange);

            for (std::size_t i = 0; i < myRange; ++i) {
                myList[i] = cpp11::writable::r_vector<T>(MyPrimeList[i]);
            }

            if (keepNames) CppConvert::SetNames(myList, myNums);
            return myList;
        }
    } else if (bAllFacs) {
        if (myRange == 1) {
            std::int64_t mPass = static_cast<std::int64_t>(myNums[0]);
            std::vector<T> myDivisors;
            bool isNegative = false;

            if (mPass < 0) {
                mPass = std::abs(mPass);
                isNegative = true;
            }

            if (mPass > 1) {
                std::vector<T> factors;
                GetPrimeFactors(mPass, factors);
                myDivisors = Factorize<T>(factors);

                if (isNegative) {
                    const std::size_t facSize = myDivisors.size();
                    std::vector<T> negPosFacs(2 * facSize);
                    std::size_t posInd = facSize, negInd = facSize - 1;

                    for (std::size_t i = 0; i < facSize; ++i, ++posInd, --negInd) {
                        negPosFacs[negInd] = -1 * myDivisors[i];
                        negPosFacs[posInd] = myDivisors[i];
                    }

                    return cpp11::writable::r_vector<T>(negPosFacs);
                }

                return cpp11::writable::r_vector<T>(myDivisors);
            } else  {
                if (isNegative) {
                    myDivisors.push_back(-1);
                } if (mPass > 0) {
                    myDivisors.push_back(1);
                }

                return cpp11::writable::r_vector<T>(myDivisors);
            }
        }

        std::vector<std::vector<T>>
            MyPrimeList(myRange, std::vector<T>());

        PollardRhoMain(myNums, myMax, bPrimeFacs, true, MyPrimeList,
                       tempVec, myRange, nThreads, maxThreads);

        cpp11::writable::list myList(myRange);

        for (std::size_t i = 0; i < myRange; ++i) {
            myList[i] = cpp11::writable::r_vector<T>(MyPrimeList[i]);
        }

        if (keepNames) CppConvert::SetNames(myList, myNums);
        return myList;
    } else {
        cpp11::sexp isPrimeVec = Rf_allocVector(LGLSXP, myRange);
        int* primeTest = INTEGER(isPrimeVec);
        std::fill_n(primeTest, myRange, 1);
        std::vector<std::vector<T>> tempList;

        PollardRhoMain(myNums, myMax, bPrimeFacs, bAllFacs, tempList,
                       primeTest, myRange, nThreads, maxThreads);

        if (keepNames) CppConvert::SetNames(isPrimeVec, myNums);
        return isPrimeVec;
    }
}

[[cpp11::register]]
SEXP PollardRhoContainer(SEXP Rv, SEXP RNamed,
                         SEXP RbPrimeFacs, SEXP RbAllFacs,
                         SEXP RNumThreads, SEXP RmaxThreads) {

    int nThreads = 1;
    int maxThreads = 1;

    CppConvert::convertPrimitive(RmaxThreads, maxThreads,
                                 VecType::Integer, "maxThreads");
    const bool bPrimeFacs = CppConvert::convertFlag(RbPrimeFacs, "bPrimeFacs");
    const bool bAllFacs = CppConvert::convertFlag(RbAllFacs, "bAllFacs");

    std::vector<double> myNums;
    bool IsNamed = CppConvert::convertFlag(RNamed, "namedList");

    if (bPrimeFacs || bAllFacs) {  // numOnly = true, checkWhole = true, negPoss = true
        CppConvert::convertVector(Rv, myNums, VecType::Numeric,
                                  "v", true, true, true);
    } else {
        CppConvert::convertVector(Rv, myNums, VecType::Numeric, "v");
    }

    if (!Rf_isNull(RNumThreads)) {
        CppConvert::convertPrimitive(
            RNumThreads, nThreads, VecType::Integer, "nThreads"
        );
    }

    if (myNums.empty()) {
        cpp11::sexp result = Rf_allocVector(
            (bPrimeFacs || bAllFacs) ? VECSXP : LGLSXP, 0
        );

        if (IsNamed) {
            CppConvert::SetNames(result, myNums);
        }

        return result;
    }

    double myMax = *std::max_element(myNums.cbegin(), myNums.cend());
    double myMin = *std::min_element(myNums.cbegin(), myNums.cend());

    if (std::abs(myMin) > myMax) {
        myMax = std::abs(myMin);
    }

    if (myMax > std::numeric_limits<int>::max()) {
        return PolRhoGlue(myNums, myMax, bPrimeFacs, bAllFacs,
                          IsNamed, nThreads, maxThreads);
    } else {
        return PolRhoGlue(myNums, static_cast<int>(myMax), bPrimeFacs,
                          bAllFacs, IsNamed, nThreads, maxThreads);
    }
}
