#include "NumbersUtils/libdivide.h"
#include "CleanConvert.h"
#include "SetUpUtils.h"
#include <thread>
#include <cmath>

template <typename T>
inline T getStartingIndex(T lowerB, T step) {

    if (step >= lowerB) {
        return (2 * step - lowerB);
    }

    T remTest = lowerB % step;

    return (remTest == 0) ? 0 : (step - remTest);
}

template <typename T, typename U>
void NumDivisorsSieve(T m, T n, T offsetStrt, U* numFacs) {

    const T myRange = offsetStrt + (n - m) + 1;
    const T sqrtBound = static_cast<T>(std::sqrt(n));

    for (T i = 2; i <= sqrtBound; ++i) {
        const T myNum = offsetStrt + (i * sqrtBound) - m;
        T j = offsetStrt + getStartingIndex(m, i);

        for (; j <= myNum; j += i) {
            ++numFacs[j];
        }

        for (; j < myRange; j += i) {
            numFacs[j] += 2;
        }
    }

    // Subtract 1 from the first entry as 1 has only itself
    // as a divisor. N.B. myRange was initialized with 2
    if (m < 2) --numFacs[0];
}

template <typename T, typename U>
void DivisorsSieve(T m, U retN, T offsetStrt,
                   std::vector<std::vector<U>> &MyDivList) {

    const T n = retN;
    constexpr T zeroOffset = 0;
    const T myRange = (n - m) + 1;

    typename std::vector<std::vector<U>>::iterator it2d;
    typename std::vector<std::vector<U>>::iterator itEnd =
        MyDivList.begin() + offsetStrt + myRange;

    std::vector<int> myMemory(myRange, 2);
    int* ptrMemory = &myMemory.front();
    NumDivisorsSieve(m, n, zeroOffset, ptrMemory);

    if (m < 2) {
        it2d = MyDivList.begin() + 1;
    } else {
        it2d = MyDivList.begin() + offsetStrt;
    }

    if (m < 2) {
        for (std::size_t i = 1; it2d < itEnd; ++it2d, ++i) {
            it2d->reserve(myMemory[i]);
            it2d->push_back(1);
        }

        MyDivList[0].push_back(1);

        for (T i = 2; i <= n; ++i) {
            for (T j = i; j <= n; j += i) {
                MyDivList[j - 1].push_back(static_cast<U>(i));
            }
        }
    } else {
        U numRet = m;
        std::vector<int> begIndex(myRange, 0);

        for (std::size_t i = 0; it2d < itEnd; ++it2d, ++i, ++numRet) {
            it2d->resize(myMemory[i]);
            it2d->back() = numRet;
            it2d->front() = 1;
            --myMemory[i];
        }

        T sqrtBound = static_cast<T>(std::sqrt(n));
        T offsetRange = myRange + offsetStrt;

        for (T i = 2; i <= sqrtBound; ++i) {
            const T myStart = getStartingIndex(m, i);
            const libdivide::divider<T> fastDiv(i);

            for (T j = myStart + offsetStrt, myNum = m + myStart,
                     memInd = myStart; j < offsetRange; j += i, myNum += i, memInd += i) {

                MyDivList[j][++begIndex[memInd]] = static_cast<U>(i);
                const T testNum = myNum / fastDiv;

                // Ensure we won't duplicate adding an element. If
                // testNum <= sqrtBound, it will be added in later
                // iterations. Also, we insert this element in the
                // pentultimate position as it will be the second
                // to the largest element at the time of inclusion.
                // E.g. let i = 5, myNum = 100, so the current
                // vectors looks like so: v = 1, 5, 10, 100 (5 was
                // added to the second position above). With i = 5,
                // testNum = 100 / 5 = 20, thus we add it to the
                // pentultimate position to give v = 1 5 10 20 100.
                if (testNum > sqrtBound) {
                    MyDivList[j][--myMemory[memInd]] = static_cast<U>(testNum);
                }
            }
        }
    }
}

template <typename T, typename U, typename V>
void DivisorMain(T myMin, U myMax, bool bDivSieve,
                 V* DivCountV, std::vector<std::vector<U>> &MyDivList,
                 std::size_t myRange, int nThreads, int maxThreads) {

    bool Parallel = false;
    T offsetStrt = 0;
    const T intMax = static_cast<T>(myMax);

    if (nThreads > 1 && maxThreads > 1  && myRange >= 20000) {
        Parallel = true;

        if (nThreads > maxThreads) {
            nThreads = maxThreads;
        }

        // Ensure that each thread has at least 10000
        if ((myRange / nThreads) < 10000) {
            nThreads = myRange / 10000;
        }
    }

    if (Parallel) {
        std::vector<std::thread> threads;
        T lowerBnd = myMin;
        const T chunkSize = myRange / nThreads;
        T upperBnd = lowerBnd + chunkSize - 1;

        for (int ind = 0; ind < (nThreads - 1); offsetStrt += chunkSize,
                    lowerBnd = (upperBnd + 1), upperBnd += chunkSize, ++ind) {
            if (bDivSieve) {
                threads.emplace_back(std::cref(DivisorsSieve<T, U>),
                                     lowerBnd, static_cast<U>(upperBnd),
                                     offsetStrt, std::ref(MyDivList));
            } else {
                threads.emplace_back(std::cref(NumDivisorsSieve<T, V>),
                                     lowerBnd, upperBnd,
                                     offsetStrt, DivCountV);
            }
        }

        if (bDivSieve) {
            threads.emplace_back(std::cref(DivisorsSieve<T, U>),
                                 lowerBnd, myMax, offsetStrt,
                                 std::ref(MyDivList));
        } else {
            threads.emplace_back(std::cref(NumDivisorsSieve<T, V>),
                                 lowerBnd, intMax, offsetStrt, DivCountV);
        }

        for (auto& thr: threads) {
            thr.join();
        }
    } else {
        if (bDivSieve) {
            DivisorsSieve(myMin, myMax, offsetStrt, MyDivList);
        } else {
            NumDivisorsSieve(myMin, intMax, offsetStrt, DivCountV);
        }
    }
}

SEXP GlueInt(int myMin, int myMax, bool bDivSieve,
             bool keepNames, int nThreads, int maxThreads) {

    std::size_t myRange = (myMax - myMin) + 1;

    if (bDivSieve) {
        std::vector<std::vector<int>> MyDivList(myRange, std::vector<int>());
        int* tempNumDivs = nullptr;

        DivisorMain(myMin, myMax, bDivSieve, tempNumDivs,
                    MyDivList, myRange, nThreads, maxThreads);

        cpp11::sexp myList = Rf_allocVector(VECSXP, myRange);

        for (std::size_t i = 0; i < myRange; ++i) {
            SET_VECTOR_ELT(myList, i, GetIntVec(MyDivList[i]));
        }

        if (keepNames) {
            SetIntNames(myList, myRange, myMin, myMax);
        }

        return myList;
    } else {
        std::vector<std::vector<int>> tempList;
        cpp11::sexp facCountV = Rf_allocVector(INTSXP, myRange);
        int* ptrFacCount  = INTEGER(facCountV);
        std::fill_n(ptrFacCount, myRange, 2);

        DivisorMain(myMin, myMax, bDivSieve, ptrFacCount,
                    tempList, myRange, nThreads, maxThreads);

        if (keepNames) {
            SetIntNames(facCountV, myRange, myMin, myMax);
        }

        return facCountV;
    }
}

SEXP GlueDbl(std::int_fast64_t myMin, double myMax,
             bool bDivSieve, bool keepNames,
             int nThreads, int maxThreads) {

    std::size_t myRange = (myMax - myMin) + 1;

    if (bDivSieve) {
        std::vector<std::vector<double>>
            MyDivList(myRange, std::vector<double>());
        double* tempNumDivs = nullptr;

        DivisorMain(myMin, myMax, bDivSieve, tempNumDivs,
                    MyDivList, myRange, nThreads, maxThreads);

        cpp11::sexp myList = Rf_allocVector(VECSXP, myRange);

        for (std::size_t i = 0; i < myRange; ++i) {
            SET_VECTOR_ELT(myList, i, GetDblVec(MyDivList[i]));
        }

        if (keepNames) {
            SetDblNames(myList, myRange, myMin, myMax);
        }

        return myList;
    } else {
        std::vector<std::vector<double>> tempList;
        cpp11::sexp facCountV = Rf_allocVector(INTSXP, myRange);
        int* ptrFacCount  = INTEGER(facCountV);
        std::fill_n(ptrFacCount, myRange, 2);

        DivisorMain(myMin, myMax, bDivSieve, ptrFacCount,
                    tempList, myRange, nThreads, maxThreads);

        if (keepNames) {
            SetDblNames(facCountV, myRange, myMin, myMax);
        }

        return facCountV;
    }
}

[[cpp11::register]]
SEXP DivNumSieveCpp(SEXP Rb1, SEXP Rb2, SEXP RbDivSieve,
                    SEXP RisNamed, SEXP RNumThreads,
                    SEXP RmaxThreads) {

    double bound1;
    double bound2;

    double myMin;
    double myMax;

    int nThreads = 1;
    int maxThreads = 1;

    CleanConvert::convertPrimitive(RmaxThreads, maxThreads,
                                   VecType::Integer, "maxThreads");
    const bool bDivSieve = CleanConvert::convertFlag(RbDivSieve,
                                                        "bDivSieve");

    const std::string namedObject = (bDivSieve) ? "namedList" : "namedVector";
    bool IsNamed = CleanConvert::convertFlag(RisNamed, namedObject);
    CleanConvert::convertPrimitive(Rb1, bound1, VecType::Numeric, "bound1");

    if (Rf_isNull(Rb2)) {
        bound2 = 1;
    } else {
        CleanConvert::convertPrimitive(Rb2, bound2, VecType::Numeric, "bound2");
    }

    if (bound1 > bound2) {
        myMax = std::floor(bound1);
        myMin = std::ceil(bound2);
    } else {
        myMax = std::floor(bound2);
        myMin = std::ceil(bound1);
    }

    if (myMax < 2) {
        if (bDivSieve) {
            cpp11::sexp res = Rf_allocVector(VECSXP, 1);
            SET_VECTOR_ELT(res, 0, GetIntVec(std::vector<int>(1, 1)));

            if (IsNamed) {
                Rf_setAttrib(res, R_NamesSymbol, Rf_mkString("1"));
            }

            return res;
        } else {
            cpp11::sexp res = Rf_allocVector(INTSXP, 1);
            INTEGER(res)[0] = 1;

            if (IsNamed) {
                Rf_setAttrib(res, R_NamesSymbol, Rf_mkString("1"));
            }

            return res;
        }
    }

    if (!Rf_isNull(RNumThreads)) {
        CleanConvert::convertPrimitive(RNumThreads, nThreads,
                                       VecType::Integer, "nThreads");
    }

    if (myMax > std::numeric_limits<int>::max()) {
        std::int_fast64_t intMin = static_cast<std::int_fast64_t>(myMin);
        return GlueDbl(intMin, myMax, bDivSieve,
                       IsNamed, nThreads, maxThreads);
    } else {
        int intMin = static_cast<int>(myMin);
        int intMax = static_cast<int>(myMax);
        return GlueInt(intMin, intMax, bDivSieve,
                       IsNamed, nThreads, maxThreads);
    }
}
