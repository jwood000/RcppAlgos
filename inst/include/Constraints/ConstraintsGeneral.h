#pragma once

#include "Constraints/ConstraintsTypes.h"
#include <vector>
#include <string>

template <typename T>
void ConstraintsGeneral(std::vector<T> &v, std::vector<int> &Reps,
                        const std::vector<std::string> &comparison,
                        std::vector<T> &cnstrntVec,
                        std::vector<T> &resVec, std::vector<T> &targetVals,
                        const std::string &myFun, const std::string &myFunTest,
                        double numRows, int n, int m, bool IsRep, bool IsComb,
                        bool IsMult, bool bUserRows, bool xtraCol,
                        ConstraintType ctype);
