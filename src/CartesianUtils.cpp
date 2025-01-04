#include "cpp11/strings.hpp"

#include "SetUpUtils.h"
#include <numeric>

double CartesianCount(const std::vector<int> &lenGrps) {
    return std::accumulate(lenGrps.begin(), lenGrps.end(),
                           1.0, std::multiplies<double>());
}

void CartesianCountGmp(mpz_class &result, const std::vector<int> &lenGrps) {

    result = 1;

    for (auto len: lenGrps) {
        result *= len;
    }
}

void CartesianInitialPrep(
    cpp11::list RList, std::vector<int> &IsFactor,
    std::vector<int> &lenGrps, int nCols
) {

    for (int i = 0; i < nCols; ++i) {
        if (Rf_isFactor(RList[i])) {
            IsFactor[i] = 1;
        } else {
            IsFactor[i] = 0;
        }

        lenGrps[i] = Rf_length(RList[i]);
    }
}

void ProductPrepare(
    cpp11::list RList, const std::vector<int> &IsFactor,
    const std::vector<int> &lenGrps, std::vector<std::vector<int>> &myVec,
    cpp11::writable::strings &charVec, std::vector<Rcomplex> &cmplxVec,
    std::vector<Rbyte> &rawVec, std::vector<double> &dblVec,
    std::vector<int> &intVec, std::vector<int> &boolVec,
    std::vector<int> &typeCheck, VecType &myType, int nCols, bool &IsDF
) {

    for (int i = 0, strt = 0; i < nCols; ++i) {
        switch(TYPEOF(RList[i])) {
            case INTSXP : {
                    if (IsFactor[i]) {
                    typeCheck[tFac] = 1;
                } else {
                    typeCheck[tInt] = 1;
                }

                std::vector<int> temp = CppConvert::GetVec<int>(RList[i]);
                std::copy(temp.begin(), temp.end(), intVec.begin() + strt);
                myType = VecType::Integer;
                break;
            } case LGLSXP : {
                std::vector<int> temp = CppConvert::GetVec<int>(RList[i]);
                std::copy(temp.begin(), temp.end(), boolVec.begin() + strt);
                typeCheck[tLog] = 1;
                myType = VecType::Logical;
                break;
            } case CPLXSXP : {
                std::vector<Rcomplex> temp =
                    CppConvert::GetVec<Rcomplex>(RList[i]);
                std::copy(temp.begin(), temp.end(), cmplxVec.begin() + strt);
                typeCheck[tCpx] = 1;
                myType = VecType::Complex;
                break;
            } case RAWSXP : {
                std::vector<Rbyte> temp = CppConvert::GetVec<Rbyte>(RList[i]);
                std::copy(temp.begin(), temp.end(), rawVec.begin() + strt);
                typeCheck[tRaw] = 1;
                myType = VecType::Raw;
                break;
            } case REALSXP : {
                std::vector<double> temp =
                    CppConvert::GetVec<double>(RList[i]);
                std::copy(temp.begin(), temp.end(), dblVec.begin() + strt);
                typeCheck[tDbl] = 1;
                myType = VecType::Numeric;
                break;
            } case STRSXP : {
                for (int j = 0; j < lenGrps[i]; ++j) {
                    charVec[strt + j] = STRING_ELT(RList[i], j);
                }

                typeCheck[tStr] = 1;
                myType = VecType::Character;
                break;
            }
        }

        std::vector<int> idx(lenGrps[i]);
        std::iota(idx.begin(), idx.end(), strt);

        myVec[i] = idx;
        strt += lenGrps[i];
    }

    int mySum = std::accumulate(typeCheck.cbegin(), typeCheck.cend(), 0);

    // We need to check to see if there is overlap in factor levels
    if (typeCheck[tFac] && mySum == 1) {
        mySum += HomoFactors(IsFactor, RList, nCols);
    }

    IsDF = (mySum > 1) ? true : false;
}

std::vector<int> nthProduct(double dblIdx, const std::vector<int> &lenGrp) {

    double index1 = dblIdx;
    const int m = lenGrp.size();

    std::vector<int> res(m);
    double temp = CartesianCount(lenGrp);

    for (int k = 0; k < m; ++k) {
        temp /= lenGrp[k];
        int j = static_cast<int>(index1 / temp);
        res[k] = j;
        index1 -= (temp * j);
    }

    for (auto &v_i: res) {
        v_i *= m;
    }

    return res;
}

std::vector<int> nthProductGmp(const mpz_class &mpzIdx,
                               const std::vector<int> &lenGrp) {

    mpz_class index1(mpzIdx);
    const int m = lenGrp.size();

    std::vector<int> res(m);
    mpz_class temp;
    mpz_class temp2;
    CartesianCountGmp(temp, lenGrp);

    for (int k = 0; k < m; ++k) {
        mpz_divexact_ui(temp.get_mpz_t(), temp.get_mpz_t(), lenGrp[k]);
        temp2 = index1 / temp;
        int j = temp2.get_si();
        res[k] = j;
        index1 -= (temp * j);
    }

    for (auto &v_i: res) {
        v_i *= m;
    }

    return res;
}

bool nextProduct(const std::vector<int> &lenGrps,
                 std::vector<int> &z, int m) {

    if (z.back() < lenGrps.back()) {
        z.back() += m;
        return true;
    } else {
        z.back() = 0;

        for (int i = m - 2; i >= 0; --i) {
            if (z[i] < lenGrps[i]) {
                z[i] += m;
                return true;
            } else {
                z[i] = 0;
            }
        }
    }

    return false;
}

bool prevProduct(const std::vector<int> &lenGrps,
                 std::vector<int> &z, int m) {

    if (z.back() > 0) {
        z.back() -= m;
        return true;
    } else {
        z.back() = lenGrps.back();

        for (int i = m - 2; i >= 0; --i) {
            if (z[i] > 0) {
                z[i] -= m;
                return true;
            } else {
                z[i] = lenGrps[i];
            }
        }
    }

    return false;
}

void GetStartProd(
    const std::vector<int> &lenNxtPr, std::vector<int> &z,
    mpz_class &lowerMpz, double &lower, int stepSize, bool IsGmp
) {

    if (IsGmp) {
        lowerMpz += stepSize;
        z = nthProductGmp(lowerMpz, lenNxtPr);
    } else {
        lower += stepSize;
        z = nthProduct(lower, lenNxtPr);
    }
}
