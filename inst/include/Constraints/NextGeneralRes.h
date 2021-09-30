#ifndef NEXT_GENERAL_RES_H
#define NEXT_GENERAL_RES_H

#include "Constraints/UserConstraintFuns.h"

template <typename T>
using nextCnstrtPtr = void (*const)(const std::vector<T> &v,
                            const std::vector<T> &targetVals,
                            const std::vector<int> &freqs,
                            const std::vector<int> &zIndex,
                            std::vector<T> &testVec, std::vector<int> &z,
                            const funcPtr<T> fun, const compPtr<T> comp,
                            int m, int m1, int m2, int nMinusM, int maxZ,
                            int pentExtreme, bool check_0, bool &check_1);
template <typename T>
nextCnstrtPtr<T> GetCnstrtPtr(bool IsMult, bool IsRep);

#endif
