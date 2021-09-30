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
void NextCnstrntDistinct(const std::vector<T> &v,
                         const std::vector<T> &targetVals,
                         const std::vector<int> &freqs,
                         const std::vector<int> &zIndex,
                         std::vector<T> &testVec, std::vector<int> &z,
                         const funcPtr<T> fun, const compPtr<T> comp,
                         int m, int m1, int m2, int nMinusM, int maxZ,
                         int pentExtreme, bool check_0, bool &check_1) {

    if (check_1) {
        bool noChange = true;

        for (int i = m2; i >= 0; --i) {
            if (z[i] != (nMinusM + i)) {
                ++z[i];
                testVec[i] = v[z[i]];

                for (int k = i + 1; k < m; ++k) {
                    z[k] = z[k - 1] + 1;
                    testVec[k] = v[z[k]];
                }

                T testVal = fun(testVec, m);
                check_0 = comp(testVal, targetVals);
                noChange = false;

                if (check_0) {
                    break;
                }
            }
        }

        check_1 = (!noChange && check_0);
    }
}

template <typename T>
void NextCnstrntMulti(const std::vector<T> &v,
                      const std::vector<T> &targetVals,
                      const std::vector<int> &freqs,
                      const std::vector<int> &zIndex,
                      std::vector<T> &testVec, std::vector<int> &z,
                      const funcPtr<T> fun, const compPtr<T> comp,
                      int m, int m1, int m2, int nMinusM, int maxZ,
                      int pentExtreme, bool check_0, bool &check_1) {

    if (check_1) {
        bool noChange = true;

        for (int i = m2; i >= 0; --i) {
            if (z[i] != freqs[pentExtreme + i]) {
                ++z[i];
                testVec[i] = v[z[i]];

                for (int j = i + 1, k = zIndex[z[i]] + 1; j <= m1; ++j, ++k) {
                    z[j] = freqs[k];
                    testVec[j] = v[z[j]];
                }

                T testVal = fun(testVec, m);
                check_0 = comp(testVal, targetVals);
                noChange = false;

                if (check_0) {
                    break;
                }
            }
        }

        check_1 = (!noChange && check_0);
    }
}

template <typename T>
void NextCnstrntRep(const std::vector<T> &v,
                    const std::vector<T> &targetVals,
                    const std::vector<int> &freqs,
                    const std::vector<int> &zIndex,
                    std::vector<T> &testVec, std::vector<int> &z,
                    const funcPtr<T> fun, const compPtr<T> comp,
                    int m, int m1, int m2, int nMinusM, int maxZ,
                    int pentExtreme, bool check_0, bool &check_1) {

    if (check_1) {
        bool noChange = true;

        for (int i = m2; i >= 0; --i) {
            if (z[i] != maxZ) {
                ++z[i];
                testVec[i] = v[z[i]];

                for (int k = i + 1; k < m; ++k) {
                    z[k] = z[k - 1];
                    testVec[k] = v[z[k]];
                }

                T testVal = fun(testVec, m);
                check_0 = comp(testVal, targetVals);
                noChange = false;

                if (check_0) {
                    break;
                }
            }
        }

        check_1 = (!noChange && check_0);
    }
}

template <typename T>
nextCnstrtPtr<T> GetCnstrtPtr(bool IsMult, bool IsRep) {
    if (IsMult) {
        return(nextCnstrtPtr<T>(NextCnstrntMulti));
    } else if (IsRep) {
        return(nextCnstrtPtr<T>(NextCnstrntRep));
    } else {
        return(nextCnstrtPtr<T>(NextCnstrntDistinct));
    }
}

template nextCnstrtPtr<int> GetCnstrtPtr(bool, bool);
template nextCnstrtPtr<double> GetCnstrtPtr(bool, bool);
