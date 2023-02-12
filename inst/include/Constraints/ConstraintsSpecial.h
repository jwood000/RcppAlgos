#pragma once

#include <vector>
#include "CppConvert/GmpxxCopy.h"

template <typename T>
void ConstraintsSpecial(const std::vector<T> &v,
                        const std::vector<T> &targetVals,
                        const std::vector<std::string> &compVec,
                        const std::vector<int> &myRep,
                        std::vector<int> freqs,
                        std::vector<T> &cnstrntVec, std::vector<T> &resVec,
                        const std::string &mainFun, std::vector<int> &z,
                        double lower, mpz_class &lowerMpz, int n, int m,
                        int maxRows, int nThreads, bool IsRep, bool xtraCol,
                        bool IsComb, bool IsMult, bool IsGmp);
