#include "Constraints/UserConstraintFuns.h"

template <typename T>
void NextCnstrntDistinct(const std::vector<T> &v,
                         const std::vector<T> &targetVals,
                         std::vector<T> &testVec, std::vector<int> &z,
                         const funcPtr<T> fun, const compPtr<T> comp, 
                         int m, int m2, int nMinusM,
                         bool check_0, bool &check_1) {
    
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
                      int m, int m1, int m2, int pentExtreme,
                      bool check_0, bool &check_1) {
    
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
                    std::vector<T> &testVec, std::vector<int> &z,
                    const funcPtr<T> fun, const compPtr<T> comp, 
                    int m, int m2, int maxZ, bool check_0, bool &check_1) {
    
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

template void NextCnstrntDistinct(const std::vector<int>&,
                                  const std::vector<int>&,
                                  std::vector<int>&, std::vector<int>&,
                                  const funcPtr<int>, const compPtr<int>, 
                                  int, int, int, bool, bool&);

template void NextCnstrntDistinct(const std::vector<double>&,
                                  const std::vector<double>&,
                                  std::vector<double>&, std::vector<int>&,
                                  const funcPtr<double>, 
                                  const compPtr<double>, 
                                  int, int, int, bool, bool&);

template void NextCnstrntMulti(const std::vector<int>&,
                               const std::vector<int>&,
                               const std::vector<int>&,
                               const std::vector<int>&,
                               std::vector<int>&, std::vector<int>&,
                               const funcPtr<int>, const compPtr<int>, 
                               int, int, int, int, bool, bool&);

template void NextCnstrntMulti(const std::vector<double>&,
                               const std::vector<double>&,
                               const std::vector<int>&,
                               const std::vector<int>&,
                               std::vector<double>&, std::vector<int>&,
                               const funcPtr<double>, const compPtr<double>, 
                               int, int, int, int, bool, bool&);

template void NextCnstrntRep(const std::vector<int>&,
                             const std::vector<int>&,
                             std::vector<int>&, std::vector<int>&,
                             const funcPtr<int>, const compPtr<int>, 
                             int, int, int, bool, bool&);

template void NextCnstrntRep(const std::vector<double>&,
                             const std::vector<double>&,
                             std::vector<double>&, std::vector<int>&,
                             const funcPtr<double>, const compPtr<double>, 
                             int, int, int, bool, bool&);
    