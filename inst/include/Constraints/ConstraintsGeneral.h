#ifndef CONNSTRAINTS_GENERAL_H
#define CONNSTRAINTS_GENERAL_H

#include <vector>
#include <string>

template <typename T>
void ConstraintsGeneral(std::vector<T> &v, std::vector<int> &Reps,
                        const std::vector<std::string> &comparison,
                        std::vector<T> &combinatoricsVec,
                        std::vector<T> &resultsVec,
                        std::vector<T> &targetVals,
                        const std::string &myFun, double numRows,
                        int n, int m, bool IsRep, bool IsComb,
                        bool IsMult, bool bUserRows, bool xtraCol);

#endif
