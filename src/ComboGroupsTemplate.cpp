#include "ComboGroups/ComboGroupsGeneral.h"
#include "ComboGroups/ComboGroupsUnique.h"
#include "ComboGroups/ComboGroupsSame.h"
#include <algorithm>
#include <numeric>
#include <memory>

ComboGroupsTemplate::ComboGroupsTemplate(
    int n_, int numGroups, int i1, int i2, int bnd
) : n(n_), r(numGroups), idx1(i1), idx2(i2), curr_bnd(bnd) {

    OneGrp = false;
}

void ComboGroupsTemplate::SetCount() {
    computedRows = numGroupCombs();
    IsGmp = computedRows > Significand53;

    if (IsGmp) {
        computedRowsMpz = numGroupCombsGmp();
    } else {
        computedRowsMpz = computedRows;
    }
}

std::unique_ptr<ComboGroupsTemplate> MakeComboGroup(
    const std::vector<int> &vGrpSize, const GroupHelper &MyGrp, int i1,
    int i2, int bnd, int grpSize, bool IsGen, bool IsUni, bool OneGrp
) {

    const int n = std::accumulate(vGrpSize.begin(), vGrpSize.end(), 0);
    const int r = vGrpSize.size();

    if (IsGen) {
        return std::make_unique<ComboGroupsGeneral>(
            n, r, i1, i2, bnd, MyGrp, OneGrp
        );
    } else if (IsUni) {
        return std::make_unique<ComboGroupsUnique>(n, r, i1, i2, bnd, vGrpSize);
    } else {
        return std::make_unique<ComboGroupsSame>(
            n, r, i1, i2, bnd + 1, grpSize
        );
    }
}

std::unique_ptr<ComboGroupsTemplate> GroupPrep(SEXP Rv, SEXP RNumGroups,
                                               SEXP RGrpSize, int n) {

    bool IsGen  = false;
    bool IsUni  = false;
    bool IsSame = true;
    bool OneGrp = false;

    int grpSize = 0;
    int numGroups = 0;
    std::vector<int> vGrpSize;

    if (Rf_isNull(RNumGroups) && Rf_isNull(RGrpSize)) {
        cpp11::stop("numGroups and grpSize cannot both be NULL");
    }

    if (!Rf_isNull(RNumGroups)) {
        CppConvert::convertPrimitive(RNumGroups, numGroups,
                                     VecType::Integer, "numGroups");
        grpSize = n / numGroups;
    }

    if (!Rf_isNull(RGrpSize)) {
        CppConvert::convertVector(RGrpSize, vGrpSize,
                                  VecType::Integer, "grpSizes");

        if (!Rf_isNull(RNumGroups) &&
            static_cast<int>(vGrpSize.size()) != numGroups) {
            cpp11::stop("numGroups and grpSizes are incompatible");
        } else {
            numGroups = vGrpSize.size();
        }

        std::vector<int> grpUni(vGrpSize.cbegin(), vGrpSize.cend());
        std::sort(grpUni.begin(), grpUni.end());
        grpUni.erase(
            std::unique(grpUni.begin(), grpUni.end()), grpUni.end()
        );

        IsUni   = grpUni.size() == vGrpSize.size() && grpUni.size() > 1;
        IsGen   = grpUni.size() != vGrpSize.size() && grpUni.size() > 1;
        IsSame  = !IsUni && !IsGen;
        grpSize = IsSame ? vGrpSize.front() : 0;
    } else {
        vGrpSize.assign(numGroups, grpSize);
    }

    std::sort(vGrpSize.begin(), vGrpSize.end());

    if (IsSame && (n % numGroups != 0)) {
        cpp11::stop("The length of v (if v is a vector) or v (if v"
                    " is a scalar) must be divisible by numGroups");
    }

    if (std::accumulate(vGrpSize.cbegin(), vGrpSize.cend(), 0) != n) {
        cpp11::stop("The sum of all group sizes must equal the length "
                    "of v (if v is a vector) or v (if v is a scalar)");
    }

    int count_one = std::count(vGrpSize.cbegin(), vGrpSize.cend(), 1);

    if (IsGen && count_one > 1) {
        OneGrp = true;
        vGrpSize.erase(vGrpSize.begin(), vGrpSize.begin() + (count_one - 1));
        vGrpSize.front() = count_one;
        numGroups -= (count_one - 1);
    }

    std::vector<int> ubound(numGroups);
    std::partial_sum(vGrpSize.begin(), vGrpSize.end(), ubound.begin());

    std::vector<int> lbound(1, 0);

    if (ubound.size() > 1) {
        lbound.insert(lbound.end(), ubound.begin(), ubound.end() - 1);
    }

    std::for_each(ubound.begin(), ubound.end(), [](int& u) {--u;});
    std::vector<bool> same(numGroups, false);

    for (int i = numGroups - 2; i >= 0; --i) {
        same[i] = vGrpSize[i] == vGrpSize[i + 1L];
    }

    GroupHelper MyGrp(vGrpSize, ubound, lbound, same);

    const int idx1 = (vGrpSize.size() > 1) ?
        std::accumulate(vGrpSize.begin(), vGrpSize.end() - 1, 0) - 1 : 0;

    const int idx2 = n - 1;
    const int lbound_cnst = (vGrpSize.size() > 2) ?
        std::accumulate(vGrpSize.begin(), vGrpSize.end() - 2, 0) : 0;

    return MakeComboGroup(vGrpSize, MyGrp, idx1, idx2, lbound_cnst,
                          grpSize, IsGen, IsUni, OneGrp);
}

void CleanV(std::vector<int> &v, const std::vector<int> &idx_used, int n) {

    v.clear();

    for (int i = 0; i < n; ++i) {
        if (!idx_used[i]) {
            v.push_back(i);
        }
    }
}

void FinishUp(const std::vector<int> &comb, std::vector<int> &v,
              std::vector<int> &res, std::vector<int> &idx_used,
              int n, int g, int j) {

    for (int k = j, i = 0; i < g; ++k, ++i) {
        res[k] = v[comb[i]];
        idx_used[res[k]] = 1;
    }

    CleanV(v, idx_used, n);
}

void SettleRes(std::vector<int> &v, std::vector<int> &res,
               std::vector<int> &idx_used, const mpz_class &mpzIdx,
               int n, int q, int g, int j, int idx) {

    const std::vector<int> comb = (g == 1) ?
                  std::vector<int>(1, idx) :
                  nthComb(q, g, idx, mpzIdx, v);

    FinishUp(comb, v, res, idx_used, n, g, j);
}

void SettleResGmp(std::vector<int> &v, std::vector<int> &res,
                  std::vector<int> &idx_used, const mpz_class &mpzIdx,
                  int n, int q, int g, int j) {

    constexpr double dblDefault = 0;

    const std::vector<int> comb = (g == 1) ?
            std::vector<int>(1, mpzIdx.get_si()) :
            nthCombGmp(q, g, dblDefault, mpzIdx, v);

    FinishUp(comb, v, res, idx_used, n, g, j);
}

void FinalTouchMisc(SEXP res, bool IsArray, int nRows,
                    bool IsNamed, const std::vector<int> &vGrpSizes,
                    const std::vector<double> &mySample,
                    const std::vector<mpz_class> &myBigSamp,
                    bool IsSample,  bool IsGmp, int r, int n) {

    std::vector<std::string> myColNames(r, "Grp");

    for (int j = 0; j < r; ++j) {
        myColNames[j] += std::to_string(j + 1);
    }

    cpp11::writable::strings myNames(n);

    for (int i = 0, k = 0; i < r; ++i) {
        for (int j = 0; j < vGrpSizes[i]; ++j, ++k) {
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
