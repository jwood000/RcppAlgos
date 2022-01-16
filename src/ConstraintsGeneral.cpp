#include "Constraints/ConstraintsClass.h"

// This function applys a constraint function to a vector v with respect
// to a constraint value "target". The main idea is that combinations are
// added successively, until a particular combination exceeds the given
// constraint value for a given constraint function. After this point, we
// can safely skip several combinations knowing that they will exceed the
// given constraint value.

// The PartitionEsque* algorithms are similar to the Constraints* algorithm.
// The main difference is in getting the next section of combinations. With
// the constraints algo, we rely on an iterative approach. The other algo
// makes heavy use of the respective "GetLowerBound" algos. Given a particular
// starting index vector as well as the specific index we are trying to set,
// these algorithms return the next lexicographical combination such that
// when the fun (e.g. "sum") is applied, the value doesn't exceed the minimum
// target value. This algo is particularly effective when we need to find
// combinations of a vector such that the value when the fun is applied is
// between a range (e.g. comparisonFun = "==" and tolerance = 100, or
// comparisonFun = c(">", "<") and limitConstraints = c(60, 62))

template <typename T>
void ConstraintsGeneral(std::vector<T> &v, std::vector<int> &Reps,
                        const std::vector<std::string> &comparison,
                        std::vector<T> &cnstrntVec,
                        std::vector<T> &resVec, std::vector<T> &targetVals,
                        const std::string &myFun, const std::string &myFunTest,
                        double numRows, int n, int m, bool IsRep, bool IsComb,
                        bool IsMult, bool bUserRows, bool xtraCol,
                        ConstraintType ctype) {

    // myFun is one of the following general functions: "prod", "sum",
    // "mean", "min", or "max"; The comparison vector contains up to 2 of the
    // following comparison operator:
    //           "<", "<=", ">", ">=", "==", ">,<", ">=,<", ">,<=", ">=,<="
    const int maxRows = std::min(dblIntMax, numRows);

    if (bUserRows) {
        cnstrntVec.reserve(m * maxRows);
        resVec.reserve(maxRows);
    }

    std::unique_ptr<ConstraintsClass<T>> Cnstrt = MakeConstraints<T>(
        comparison, myFun, myFunTest, Reps, targetVals,
        ctype, n, m, IsComb, xtraCol, IsMult, IsRep
    );

    for (auto comp: comparison) {
        Cnstrt->Prepare(comp, v);
        Cnstrt->GetSolutions(v, targetVals, cnstrntVec, resVec, maxRows);
        targetVals.erase(targetVals.begin());
    }
}

template void ConstraintsGeneral(std::vector<int>&, std::vector<int>&,
                                 const std::vector<std::string>&,
                                 std::vector<int>&, std::vector<int>&,
                                 std::vector<int>&, const std::string&,
                                 const std::string&, double, int, int, bool,
                                 bool, bool, bool, bool, ConstraintType);

template void ConstraintsGeneral(std::vector<double>&, std::vector<int>&,
                                 const std::vector<std::string>&,
                                 std::vector<double>&, std::vector<double>&,
                                 std::vector<double>&, const std::string&,
                                 const std::string&, double, int, int, bool,
                                 bool, bool, bool, bool, ConstraintType);
