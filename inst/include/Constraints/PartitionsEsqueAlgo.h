#ifndef PARTITION_ESQUE_ALGO_H
#define PARTITION_ESQUE_ALGO_H

#include <vector>
#include <string>

template <typename T>
void PartitionsEsqueAlgo(std::vector<T> &v,
                         const std::vector<T> &targetVals,
                         const std::vector<int> &Reps,
                         const std::string &myFun,
                         const std::string &comparison,
                         std::vector<T> &combinatoricsVec,
                         std::vector<T> &resultsVec, int maxRows,
                         int n, int m, bool isRep, bool IsComb,
                         bool xtraCol, bool IsMult, bool bUserRows);

#endif