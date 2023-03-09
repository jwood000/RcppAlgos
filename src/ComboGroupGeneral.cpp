#include "ComboGroup/ComboGroupGeneral.h"
#include <unordered_map>
#include <algorithm>

bool nextCmbGrpGen(std::vector<int> &z, int idx1, int idx2,
                   int curr_bnd, const Group &MyGrp) {

    while (idx2 > idx1 && z[idx2] > z[idx1]) {
        --idx2;
    }

    if ((idx2 + 1) < static_cast<int>(z.size())) {
        std::swap(z[idx1], z[idx2 + 1]);
        return true;
    }

    bool clean = true;

    // Start at penultimate group
    for (int i = MyGrp.get_size() - 2; i >= 0; --i) {

        const int tipPnt = z[idx2];
        const int low    = MyGrp.get_low(curr_bnd, i);

        while (idx1 > low && tipPnt < z[idx1]) {
            --idx1;
        }

        if (z[idx1] < tipPnt) {
            MyGrp.balance(z, idx1, curr_bnd, i);
            return true;
        } else if (clean && MyGrp.require_external(z, i)) {
            if (MyGrp.flip_external(z, idx1, i)) return true;
            clean = false;
        } else if (i > 0) {
            MyGrp.step(idx1, idx2, curr_bnd, i);
        }
    }

    return false;
}

double numCmbGrpGen(const std::vector<int> &grp, int n) {

    double result = 1;
    int curr_size = n;
    std::unordered_map<int, int> table;

    for (auto g: grp) {
        result    *= nChooseK(curr_size, g);
        curr_size -= g;
        ++table[g];
    }

    if (result < std::numeric_limits<double>::max()) {
        double myDiv = 1;

        for (auto const &t: table) {
            myDiv *= std::tgamma(t.second + 1);
        }

        result /= myDiv;
        return std::round(result);
    } else {
        return std::numeric_limits<double>::infinity();
    }
}

mpz_class numCmbGrpGenGmp(const std::vector<int> &grp, int n) {

    mpz_class result(1);
    mpz_class temp(1);

    int curr_size = n;
    std::unordered_map<int, int> table;

    for (auto g: grp) {
        nChooseKGmp(temp, curr_size, g);
        result    *= temp;
        curr_size -= g;
        ++table[g];
    }

    mpz_class myDiv(1);

    for (auto const &t: table) {
        mpz_fac_ui(temp.get_mpz_t(), t.second);
        myDiv *= temp;
    }

    mpz_divexact(result.get_mpz_t(), result.get_mpz_t(), myDiv.get_mpz_t());
    return result;
}

ComboGroupGeneral::ComboGroupGeneral(
    int n_, int numGroups, int i1, int i2, int bnd, Group MyGrp_
) : ComboGroup(n_, numGroups, i1, i2, bnd), MyGrp(MyGrp_) {};

bool ComboGroupGeneral::nextComboGroup(std::vector<int> &z) {
    return nextCmbGrpGen(z, idx1, idx2, curr_bnd, MyGrp);
}

double ComboGroupGeneral::numGroupCombs() {
    return numCmbGrpGen(MyGrp.grp, n);
}

mpz_class ComboGroupGeneral::numGroupCombsGmp() {
    return numCmbGrpGenGmp(MyGrp.grp, n);
}

#include <iostream>

void removeFirstSet(std::vector<int> &v, int &p) {

    int nSame = 1;

    for (int i = 1, v_s = v.size(); i < v_s && v.front() == v[i]; ++i) {
        ++nSame;
    }

    if (v.size() >= nSame) {
        p -= (v.front() * nSame);
        v.erase(v.begin(), v.begin() + nSame);
    }
}

double intermediate(int nGrps, int grpS, int n) {

    double result = 1;

    for (int i = 0; i < nGrps; ++i, n -= grpS) {
        result *= nChooseK(n, grpS);
    }

    if (nGrps > 1) result /= std::tgamma(nGrps + 1);
    return result;
}

void ResolveSet(std::vector<int> &v, std::vector<int> &res,
                std::vector<int> &idx_used, const mpz_class &mpzIdx,
                int n, int q, int g, int k, int idx, int setSize) {

    int grp_idx = 0;
    int q_g = q - g;
    int p   = q - 1;

    const int g1   = g - 1;
    std::int64_t m = nChooseK(p, g1);
    int curr_nGrps = setSize - 1;

    for (int i = 0; i < (setSize - 1); ++i, grp_idx = 0) {

        std::int64_t secLen = static_cast<std::int64_t>(
            intermediate(curr_nGrps, g, q_g)
        );

        while ((secLen * m) < idx) {
            idx     -= (secLen * m);
            grp_idx += m;

            --q_g;
            --p;

            m = nChooseK(p, g1);
            secLen = static_cast<std::int64_t>(
                intermediate(curr_nGrps, g, q_g)
            );
        }

        grp_idx += (idx / secLen);
        SettleRes(v, res, idx_used, mpzIdx, n, q, g, k, grp_idx);

        for (int j = 0; j < res[k]; ++j) {
            idx_used[j] = 1;
        }

        CleanV(v, idx_used, n);

        q   = v.size();
        k  += g;
        q_g = q - g;
        p   = q - 1;
        m   = nChooseK(p, g1);

        --curr_nGrps;
        idx -= ((idx / secLen) * secLen);
    }

    while (p > 0 && p < idx) {
        idx     -= p;
        grp_idx += p;
        --p;
    }

    grp_idx += idx;
    SettleRes(v, res, idx_used, mpzIdx, n, q, g, k, grp_idx);
    k += g;

    std::fill(idx_used.begin(), idx_used.end(), 0);

    for (int i = 0; i < k; ++i) {
        idx_used[res[i]] = 1;
    }

    CleanV(v, idx_used, n);
}

std::vector<int> ComboGroupGeneral::nthComboGroup(double myIndex) {

    int p = n;
    int q = n;

    std::vector<int> grpCopy(MyGrp.grp.begin(), MyGrp.grp.end());

    // This vector represents the "group of groups." That is, if:
    //
    //    MyGrp.grp = {2, 2, 2, 3, 4, 4, 5};
    //
    // Note, the vector will always be sorted (see GroupPrep in
    // ComboGroupClass.cpp). Thus grpSets will be:
    //
    //    grpSets = {3, 1, 2, 1};
    std::vector<int> grpSets;

    for (int i = 0, j = 0; i < r; ++i) {
        if (i > 0 && MyGrp.grp[i] == MyGrp.grp[i - 1]) {
            ++grpSets[j - 1];
        } else {
            grpSets.push_back(1);
            ++j;
        }
    }

    const int nSets = grpSets.size();
    std::int64_t intIdx = myIndex;
    mpz_class mpzDefault;

    std::vector<int>  res(n, 0);
    std::vector<int> idx_used(n, 0);
    std::vector<int>  v(n);
    std::iota(v.begin(), v.end(), 0);

    for (int i = 0, j = 0, k = 0; i < nSets; j += grpSets[i], ++i) {

        removeFirstSet(grpCopy, p);
        const std::int64_t secLen = static_cast<std::int64_t>(
            numCmbGrpGen(grpCopy, p)
        );

        const std::int64_t idx = intIdx / secLen;
        const int g = MyGrp.grp[j];

        if (grpSets[i] == 1) {
            SettleRes(v, res, idx_used, mpzDefault, n, q, g, k, idx);
        } else {
            ResolveSet(v, res, idx_used, mpzDefault,
                       n, q, g, k, idx, grpSets[i]);
        }

        intIdx -= (idx * secLen);
        k += (g * grpSets[i]);
        q = p;
    }

    return res;
}

std::vector<int> ComboGroupGeneral::nthComboGroupGmp(
    const mpz_class &lowerMpz
) {

    mpz_class ind1(lowerMpz);
    mpz_class ind2(lowerMpz);

    int s = n - 1;
    const int g = MyGrp.grp.front() - 1;

    mpz_class temp(1);
    mpz_class secLen(1);

    nChooseKGmp(temp, s, g);
    mpz_divexact(secLen.get_mpz_t(), computedRowsMpz.get_mpz_t(),
                 temp.get_mpz_t());

    std::vector<int>  res(n, 0);
    std::vector<char> idx_used(n, '0');
    std::vector<int>  v(s);
    std::iota(v.begin(), v.end(), 1);

    int myMin = 0;
    constexpr double dblDefault = 0;

    for (int j = 0; j < (r - 1); ++j) {
        ind2 /= secLen;
        res[j * MyGrp.grp[j]] = myMin;
        idx_used[myMin] = '1';

        const std::vector<int> comb = (g == 1) ?
        std::vector<int>(1, ind2.get_si()) :
            nthCombGmp(s, g, dblDefault, ind2, v);

        for (int k = j * MyGrp.grp[j] + 1, i = 0; i < g; ++k, ++i) {
            res[k] = v[comb[i]];
            idx_used[res[k]] = '1';
        }

        v.clear();

        for (int i = 1; i < n; ++i) {
            if (!idx_used[i]) {
                v.push_back(i);
            }
        }

        myMin = v.front();
        v.erase(v.begin());
        temp = ind2 * secLen;
        ind1 -= temp;
        ind2 = ind1;

        s -= MyGrp.grp[j];
        nChooseKGmp(temp, s, g);
        mpz_divexact(secLen.get_mpz_t(), secLen.get_mpz_t(), temp.get_mpz_t());
    }

    res[(r - 1) * MyGrp.grp.back()] = myMin;

    for (int k = (r - 1) * MyGrp.grp.back() + 1, i = 0; i < g; ++k, ++i) {
        res[k] = v[i];
    }

    return res;
}

void ComboGroupGeneral::FinalTouch(
    SEXP res, bool IsArray, int nRows, bool IsNamed,
    const std::vector<double> &mySample,
    const std::vector<mpz_class> &myBigSamp, bool IsSample
) {
    FinalTouchMisc(res, IsArray, nRows, IsNamed, MyGrp.grp,
                   mySample, myBigSamp, IsSample, IsGmp, r, n);
}
