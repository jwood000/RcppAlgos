#include "ComboGroup/ComboGroupUnique.h"

// This algorithm is for the case where the length of each group is different.
// The skeleton of the algorithm is exactly the same as the special case. The
// only difference is that we are aware of what group we are in as well as the
// size of each group, hence the vector vGrpSize.
bool nextCmbGrpUni(std::vector<int> &z,
                   const std::vector<int> &vGrpSize,
                   int idx1, int idx2, int lbound) {

    while (idx2 > idx1 && z[idx2] > z[idx1]) {
        --idx2;
    }

    if ((idx2 + 1) < static_cast<int>(z.size())) {
        std::swap(z[idx1], z[idx2 + 1]);
        return true;
    }

    const auto zbeg = z.begin();

    // Start at penultimate group
    for (int g_idx = vGrpSize.size() - 2; g_idx >= 0; --g_idx, --idx1) {

        const int tipPnt = z[idx2];

        while (idx1 > lbound && tipPnt < z[idx1]) {
            --idx1;
        }

        if (z[idx1] < tipPnt) {
            int idx3 = idx1 + 1;
            std::sort(zbeg + idx3, z.end());

            const int len_rng = lbound + vGrpSize[g_idx] - idx1;

            while (z[idx3] < z[idx1]) {
                ++idx3;
            }

            std::swap(z[idx3], z[idx1]);
            std::rotate(zbeg + idx1 + 1,
                        zbeg + idx3 + 1, zbeg + idx3 + len_rng);
            return true;
        } else if (g_idx > 0) {
            idx2   -= vGrpSize[g_idx + 1];
            lbound -= vGrpSize[g_idx - 1];
        }
    }

    return false;
}

double numCmbGrpUni(const std::vector<int> &grp, int n) {

    double result = std::tgamma(n + 1);

    if (result < std::numeric_limits<double>::max()) {
        double myDiv = 1;

        for (auto g: grp) {
            myDiv *= std::tgamma(g + 1);
        }

        result /= myDiv;
        return std::round(result);
    } else {
        return std::numeric_limits<double>::infinity();
    }
}

mpz_class numCmbGrpUniGmp(const std::vector<int> &grp, int n) {

    mpz_class result(1);
    mpz_fac_ui(result.get_mpz_t(), n);

    mpz_class myDiv(1);
    mpz_class temp(1);

    for (auto g: grp) {
        mpz_fac_ui(temp.get_mpz_t(), g);
        myDiv *= temp;
    }

    mpz_divexact(result.get_mpz_t(), result.get_mpz_t(), myDiv.get_mpz_t());
    return result;
}

ComboGroupUnique::ComboGroupUnique(
    int n_, int numGroups, int i1, int i2,
    int bnd, const std::vector<int> &grp_
) : ComboGroup(n_, numGroups, i1, i2, bnd), grp(grp_) {}

bool ComboGroupUnique::nextComboGroup(std::vector<int> &z) {
    return nextCmbGrpUni(z, grp, idx1, idx2, curr_bnd);
}

double ComboGroupUnique::numGroupCombs() {
    return numCmbGrpUni(grp, n);
}

mpz_class ComboGroupUnique::numGroupCombsGmp() {
    return numCmbGrpUniGmp(grp, n);
}

void removeFirst(std::vector<int> &v, int &a) {
    if (v.size()) {
        a -= v.front();
        v.erase(v.begin());
    }
}

std::vector<int> ComboGroupUnique::nthComboGroup(double myIndex) {

    int p = n;
    int q = n;

    std::vector<int> grpCopy(grp.begin(), grp.end());
    std::int64_t intIdx = myIndex;
    mpz_class mpzDefault;

    std::vector<int> res(n, 0);
    std::vector<int> idx_used(n, 0);
    std::vector<int> v(n);
    std::iota(v.begin(), v.end(), 0);

    for (int i = 0, j = 0, g = grp.front();
         i < (r - 1); ++i, j += g, g = grp[i]) {

        removeFirst(grpCopy, p);
        const std::int64_t secLen = static_cast<std::int64_t>(
            numCmbGrpUni(grpCopy, p)
        );
        const std::int64_t idx = intIdx / secLen;

        SettleRes(v, res, idx_used, mpzDefault, n, q, g, j, idx);
        q = p;
        intIdx -= (idx * secLen);
    }

    for (int k = n - 1, i = 1, s = v.size() - 1;
         i <= grp.back(); --k, --s, ++i) {
        res[k] = v[s];
    }

    return res;
}

std::vector<int> ComboGroupUnique::nthComboGroupGmp(
    const mpz_class &lowerMpz
) {

    int p = n;
    int q = n;

    std::vector<int> grpCopy(grp.begin(), grp.end());
    std::vector<int> res(n, 0);
    std::vector<int> idx_used(n, 0);

    std::vector<int>  v(n);
    std::iota(v.begin(), v.end(), 0);

    mpz_class mpzIndex(lowerMpz);
    mpz_class secLen(1);
    mpz_class idx(1);

    for (int i = 0, j = 0, g = grp.front();
         i < (r - 1); ++i, j += g, g = grp[i]) {

        removeFirst(grpCopy, p);
        secLen = numCmbGrpUniGmp(grpCopy, p);
        idx = mpzIndex / secLen;

        SettleResGmp(v, res, idx_used, idx, n, q, g, j);
        q = p;
        mpzIndex -= (idx * secLen);
    }

    for (int k = n - 1, i = 1, s = v.size() - 1;
         i <= grp.back(); --k, --s, ++i) {
        res[k] = v[s];
    }

    return res;
}

void ComboGroupUnique::FinalTouch(
    SEXP res, bool IsArray, int nRows, bool IsNamed,
    const std::vector<double> &mySample,
    const std::vector<mpz_class> &myBigSamp, bool IsSample
) {
    FinalTouchMisc(res, IsArray, nRows, IsNamed, grp,
                   mySample, myBigSamp, IsSample, IsGmp, r, n);
}
