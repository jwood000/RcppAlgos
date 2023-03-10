#include "ComboGroup/ComboGroupSame.h"

bool nextCmbGrpSame(std::vector<int> &z, int r, int grpSize,
                    int idx1, int idx2, int curr_bnd, int n) {

    while (idx2 > idx1 && z[idx2] > z[idx1]) {
        --idx2;
    }

    if ((idx2 + 1) < n) {
        // Formally, we had the following conditional:
        //
        // if (z[idx2 + 1] > z[idx1]) {std::swap(z[idx1], z[idx2 + 1]);}
        //
        // This is unnecessary because we know that idx2 always starts out as
        // z.size() - 1 and if the while loop is never entered, idx2 + 1 will
        // be equal to z.size(), meaning:
        //
        //       (idx2 + 1) < n -->> false
        //
        // and this section will not be entered.
        //
        // Otherwise, both of the checks are true at least once. That is,
        // idx2 > idx1 and z[idx2] > z[idx1] for at least one iteration. When
        // the check becomes false, we know that the previous iteration was
        // true, thus z[idx2 + 1] > z[idx1] is true. Thus we can safely
        // proceed with our swap.

        std::swap(z[idx1], z[idx2 + 1]);
        return true;
    }

    const auto zbeg = z.begin();

    while (idx1 > 0) {
        const int tipPnt = z[idx2];

        while (idx1 > curr_bnd && tipPnt < z[idx1]) {
            --idx1;
        }

        if (tipPnt > z[idx1]) { // **Crucial Part**
            int idx3 = idx1 + 1;
            std::sort(zbeg + idx3, z.end());

            // length of left range plus one. The plus one is needed as we are
            // rotating at a pivot just to the right (i.e. plus one).
            const int len_rng = curr_bnd + grpSize - idx3;

            while (z[idx3] < z[idx1]) {
                ++idx3;
            }

            std::swap(z[idx3], z[idx1]);
            std::rotate(zbeg + idx1 + 1,
                        zbeg + idx3 + 1, zbeg + idx3 + len_rng);
            return true;
        } else {
            idx1     -= 2;
            idx2     -= grpSize;
            curr_bnd -= grpSize;
        }
    }

    return false;
}

ComboGroupSame::ComboGroupSame(
    int n_, int numGroups, int i1, int i2, int bnd, int size
) : ComboGroup(n_, numGroups, i1, i2, bnd), grpSize(size) {};

bool ComboGroupSame::nextComboGroup(std::vector<int> &z) {
    return nextCmbGrpSame(z, r, grpSize, idx1, idx2, curr_bnd, n);
}

double ComboGroupSame::numGroupCombs() {

    double result = 1;

    for (double i = n; i > r; --i) {
        result *= i;
    }

    if (result < std::numeric_limits<double>::max()) {
        double myDiv = 1;

        for (double i = 2; i <= grpSize; ++i) {
            myDiv *= i;
        }

        result /= std::pow(myDiv, r);
        return std::round(result);
    } else {
        return std::numeric_limits<double>::infinity();
    }
}

mpz_class ComboGroupSame::numGroupCombsGmp() {

    mpz_class result(1);

    for (int i = n; i > r; --i) {
        result *= i;
    }

    mpz_class myDiv(1);

    for (int i = 2; i <= grpSize; ++i) {
        myDiv *= i;
    }

    mpz_pow_ui(myDiv.get_mpz_t(), myDiv.get_mpz_t(), r);
    mpz_divexact(
        result.get_mpz_t(), result.get_mpz_t(), myDiv.get_mpz_t()
    );

    return result;
}

std::vector<int> ComboGroupSame::nthComboGroup(double myIndex) {

    int s = n - 1;
    const int g = grpSize - 1;
    std::int64_t temp   = static_cast<std::int64_t>(nChooseK(s, g));
    std::int64_t secLen = static_cast<std::int64_t>(computedRows) / temp;

    std::vector<int> res(n, 0);
    std::vector<int> idx_used(n, 0);
    std::vector<int> v(s);
    std::iota(v.begin(), v.end(), 1);

    int myMin = 0;
    std::int64_t ind1 = myIndex;
    std::int64_t ind2 = myIndex;

    mpz_class mpzDefault;

    for (int j = 0; j < (r - 1); ++j) {
        ind2 = ind2 / secLen;
        res[j * grpSize] = myMin;
        idx_used[myMin]  = 1;

        SettleRes(v, res, idx_used, mpzDefault,
                  n, s, g, j * grpSize + 1, ind2);

        myMin = v.front();
        v.erase(v.begin());
        ind1 -= ind2 * secLen;
        ind2 = ind1;

        s -= grpSize;
        temp = nChooseK(s, g);
        secLen /= temp;
    }

    res[(r - 1) * grpSize] = myMin;

    for (int k = (r - 1) * grpSize + 1, i = 0; i < g; ++k, ++i) {
        res[k] = v[i];
    }

    return res;
}

std::vector<int> ComboGroupSame::nthComboGroupGmp(const mpz_class &lowerMpz) {

    mpz_class ind1(lowerMpz);
    mpz_class ind2(lowerMpz);

    int s = n - 1;
    const int g = grpSize - 1;

    mpz_class temp(1);
    mpz_class secLen(1);

    nChooseKGmp(temp, s, g);
    mpz_divexact(
        secLen.get_mpz_t(), computedRowsMpz.get_mpz_t(), temp.get_mpz_t()
    );

    std::vector<int> res(n, 0);
    std::vector<int> idx_used(n, 0);

    std::vector<int>  v(s);
    std::iota(v.begin(), v.end(), 1);
    int myMin = 0;

    for (int j = 0; j < (r - 1); ++j) {
        ind2 /= secLen;
        res[j * grpSize] = myMin;
        idx_used[myMin]  = 1;

        SettleResGmp(v, res, idx_used, ind2,
                     n, s, g, j * grpSize + 1);

        myMin = v.front();
        v.erase(v.begin());
        temp = ind2 * secLen;
        ind1 -= temp;
        ind2 = ind1;

        s -= grpSize;
        nChooseKGmp(temp, s, g);
        mpz_divexact(
            secLen.get_mpz_t(), secLen.get_mpz_t(), temp.get_mpz_t()
        );
    }

    res[(r - 1) * grpSize] = myMin;

    for (int k = (r - 1) * grpSize + 1, i = 0; i < g; ++k, ++i) {
        res[k] = v[i];
    }

    return res;
}

void ComboGroupSame::FinalTouch(
    SEXP res, bool IsArray, int nRows, bool IsNamed,
    const std::vector<double> &mySample,
    const std::vector<mpz_class> &myBigSamp, bool IsSample
) {

    std::vector<std::string> myColNames(r, "Grp");

    for (int j = 0; j < r; ++j) {
        myColNames[j] += std::to_string(j + 1);
    }

    if (IsArray) {
        cpp11::integers dim({nRows, grpSize, r});
        Rf_setAttrib(res, R_DimSymbol, dim);
        cpp11::writable::strings myNames(r);

        for (int i = 0; i < r; ++i) {
            myNames[i] = myColNames[i].c_str();
        }

        SetSampleNames(res, IsGmp, nRows, mySample,
                       myBigSamp, IsNamed, myNames, 2);

        if (!IsNamed) {
            cpp11::writable::list dimNames(3);
            dimNames[2] = myNames;
            Rf_setAttrib(res, R_DimNamesSymbol, dimNames);
        }
    } else {
        cpp11::writable::strings myNames(n);

        for (int i = 0, k = 0; i < r; ++i) {
            for (int j = 0; j < grpSize; ++j, ++k) {
                myNames[k] = myColNames[i].c_str();
            }
        }

        SetSampleNames(res, IsGmp, nRows, mySample,
                       myBigSamp, IsNamed, myNames, 1);

        if (!IsNamed) {
            cpp11::writable::list dimNames(2);
            dimNames[1] = myNames;
            Rf_setAttrib(res, R_DimNamesSymbol, dimNames);
        }
    }
}
