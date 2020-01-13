#ifndef PARTITION_ESQUE_ALGO_H
#define PARTITION_ESQUE_ALGO_H

#include <vector>
#include <string>

template <typename typeRcpp, typename typeVector>
typeRcpp PartitionEsqueAlgo(int n, int m, std::vector<typeVector> &v, bool isRep, 
                            std::string myFun, const std::string &comparison,
                            std::vector<typeVector> targetVals, double numRows, 
                            bool isComb, bool xtraCol, std::vector<int> &Reps, 
                            bool isMult, bool bUserRows);

#endif
