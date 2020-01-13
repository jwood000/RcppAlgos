#ifndef CONSTRAINTS_GENERAL_H
#define CONSTRAINTS_GENERAL_H

#include <vector>
#include <string>

template <typename typeRcpp, typename typeVector>
typeRcpp CombinatoricsConstraints(int n, int m, std::vector<typeVector> &v, bool isRep, 
                                  std::string myFun, std::vector<std::string> comparison,
                                  std::vector<typeVector> targetVals, double numRows, 
                                  bool isComb, bool xtraCol, std::vector<int> &Reps, 
                                  bool isMult, bool bUserRows);

#endif
