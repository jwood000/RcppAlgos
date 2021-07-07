#ifndef CONSTRAINTS_SPECIAL_H
#define CONSTRAINTS_SPECIAL_H

#include <vector>

template <typename T>
void ConstraintsSpecial(const std::vector<T> &v,
                        const std::vector<T> &targetVals,
                        const std::vector<std::string> &compVec,
                        const std::vector<int> &freqs,
                        std::vector<T> &cnstrntVec, std::vector<T> &resVec,
                        const std::string &mainFun, std::vector<int> &z,
                        int n, int m, int maxRows, bool IsRep, 
                        bool xtraCol, bool IsComb, bool IsMult);

#endif