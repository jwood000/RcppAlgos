#include "Constraints/ConstraintsClass.h"

// This function applys a constraint function to a vector v with respect
// to a constraint value "target". The main idea is that combinations are
// added successively, until a particular combination exceeds the given
// constraint value for a given constraint function. After this point, we
// can safely skip several combinations knowing that they will exceed the
// given constraint value.

template <typename T>
void ConstraintsGeneral(std::vector<T> &v, std::vector<int> &Reps,
                        const std::vector<std::string> &comparison,
                        std::vector<T> &cnstrntVec,
                        std::vector<T> &resVec, std::vector<T> &targetVals,
                        const std::string &myFun, double numRows,
                        int n, int m, bool IsRep, bool IsComb,
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
        comparison, myFun, Reps, targetVals, ctype,
        n, m, IsComb, xtraCol, IsMult, IsRep
    );

    for (std::size_t nC = 0; nC < comparison.size(); ++nC) {
        Cnstrt->Prepare(comparison[nC], v);
        Cnstrt->GetSolutions(v, targetVals, cnstrntVec, resVec, maxRows);
        targetVals.erase(targetVals.begin());
    }
}

template void ConstraintsGeneral(std::vector<int>&, std::vector<int>&,
                                 const std::vector<std::string>&,
                                 std::vector<int>&, std::vector<int>&,
                                 std::vector<int>&, const std::string&,
                                 double, int, int, bool, bool,
                                 bool, bool, bool, ConstraintType);

template void ConstraintsGeneral(std::vector<double>&, std::vector<int>&,
                                 const std::vector<std::string>&,
                                 std::vector<double>&, std::vector<double>&,
                                 std::vector<double>&, const std::string&,
                                 double, int, int, bool, bool,
                                 bool, bool, bool, ConstraintType);
