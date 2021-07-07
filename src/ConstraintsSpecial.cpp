#include "Constraints/UserConstraintFuns.h"
#include "Constraints/ConstraintsUtils.h"
#include "NextCombinatorics.h"

// This is called when we can't easily produce a (loose) monotonic sequence overall,
// and we must generate and test every possible combination/permutation. This occurs
// when we are using "prod" and we have negative numbers involved. We also call this
// when lower is invoked implying that we are testing a specific range.
template <typename T>
void ConstraintsSpecial(const std::vector<T> &v,
                        const std::vector<T> &targetVals,
                        const std::vector<std::string> &compVec,
                        const std::vector<int> &freqs,
                        std::vector<T> &cnstrntVec, std::vector<T> &resVec,
                        const std::string &mainFun, std::vector<int> &z,
                        int n, int m, int maxRows, bool IsRep, 
                        bool xtraCol, bool IsComb, bool IsMult) {

    // Needed to determine if nextFullPerm or nextPerm will be called
    const bool IsFullPerm = (IsComb || IsRep) ? false :
                            (m == n || m == static_cast<int>(freqs.size()));
    
    std::vector<T> testVec(m);
    const funcPtr<T> fun = GetFuncPtr<T>(mainFun);
    const compPtr<T> compOne = GetCompPtr<T>(compVec.front());
    const nextIterPtr nextIter = GetNextIterPtr(IsComb, IsMult,
                                                IsRep, IsFullPerm);
    
    const int n1 = IsComb ? n - 1 : (IsMult ? freqs.size() - 1 : n - 1);
    const int m1 = m - 1;

    if (compVec.size() == 1) {
        for (int i = 0; i < maxRows; ++i, nextIter(freqs, z, n1, m1)) {
            for (int i = 0; i < m; ++i) {
                testVec[i] = v[z[i]];
            }

            const T testVal = fun(testVec, m);

            if (compOne(testVal, targetVals)) {
                cnstrntVec.insert(cnstrntVec.end(),
                                  testVec.begin(), testVec.end());

                if (xtraCol) {
                    resVec.push_back(testVal);
                }
            }
        }
    } else {
        compPtr<T> compTwo = GetCompPtr<T>(compVec.back());
        std::vector<T> targetVals2(1, targetVals.back());

        for (int i = 0; i < maxRows; ++i, nextIter(freqs, z, n1, m1)) {
            for (int i = 0; i < m; ++i) {
                testVec[i] = v[z[i]];
            }
            
            const T testVal = fun(testVec, m);

            if (compOne(testVal, targetVals) ||
                compTwo(testVal, targetVals2)) {
                
                cnstrntVec.insert(cnstrntVec.end(),
                                  testVec.begin(), testVec.end());
                
                if (xtraCol) {
                    resVec.push_back(testVal);
                }
            }
        }
    }
}

template void ConstraintsSpecial(const std::vector<int>&,
                                 const std::vector<int>&,
                                 const std::vector<std::string>&,
                                 const std::vector<int>&, std::vector<int>&,
                                 std::vector<int>&, const std::string&,
                                 std::vector<int>&, int, int,
                                 int, bool, bool, bool, bool);

template void ConstraintsSpecial(const std::vector<double>&,
                                 const std::vector<double>&,
                                 const std::vector<std::string>&,
                                 const std::vector<int>&,
                                 std::vector<double>&, std::vector<double>&,
                                 const std::string&, std::vector<int>&,
                                 int, int, int, bool, bool, bool, bool);
