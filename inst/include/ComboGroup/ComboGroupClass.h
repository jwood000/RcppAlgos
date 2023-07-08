#pragma once

#include "cpp11/strings.hpp"
#include "cpp11/list.hpp"

#include "Combinations/NthCombination.h"
#include "Combinations/BigComboCount.h"
#include "Combinations/ComboCount.h"
#include "ComboGroup/GroupClass.h"
#include "SetUpUtils.h"
#include <numeric>
#include <memory>
#include <cmath>

class ComboGroup {
protected:

    const int n; // Size of vector which is also the size of z (i.e. z.size())
    const int r; // Number of groups

    const int idx1;
    const int idx2;
    const int curr_bnd;

    bool IsGmp;
    double computedRows;
    mpz_class computedRowsMpz;

public:

    virtual ~ComboGroup() = default;
    ComboGroup(int n_, int numGroups, int i1, int i2, int bnd);

    virtual bool nextComboGroup(std::vector<int> &z) = 0;
    virtual double numGroupCombs() = 0;
    virtual mpz_class numGroupCombsGmp() = 0;
    virtual std::vector<int> nthComboGroup(double myIndex) = 0;
    virtual std::vector<int> nthComboGroupGmp(const mpz_class &lowerMpz) = 0;

    virtual void FinalTouch(
        SEXP res, bool IsArray, int nRows, bool IsNamed,
        const std::vector<double> &mySample,
        const std::vector<mpz_class> &myBigSamp, bool IsSample
    ) = 0;

    void SetCount();
    double GetDblCount() const {return computedRows;}
    mpz_class GetMpzCount() const {return computedRowsMpz;}
    bool GetIsGmp() const {return IsGmp;}
    int GetNumGrps() const {return r;}

    SEXP GetCount() const {
        return CppConvert::GetCount(IsGmp, computedRowsMpz, computedRows);
    };
};

std::unique_ptr<ComboGroup> MakeComboGroup(
    const std::vector<int> &vGrpSize, Group *const MyGrp,
    int i1, int i2, int bnd, int grpSize, bool IsGen, bool IsUni
);

std::unique_ptr<ComboGroup> GroupPrep(
    std::vector<int> &vInt, std::vector<double> &vNum, int &n,
    VecType &myType, SEXP Rv, SEXP RNumGroups, SEXP RGrpSize
);

void CleanV(std::vector<int> &v, const std::vector<int> &idx_used, int n);

void SettleRes(std::vector<int> &v, std::vector<int> &res,
               std::vector<int> &idx_used, const mpz_class &mpzIdx,
               int n, int q, int g, int j, int idx);

void SettleResGmp(std::vector<int> &v, std::vector<int> &res,
                  std::vector<int> &idx_used, const mpz_class &mpzIdx,
                  int n, int q, int g, int j);

void FinalTouchMisc(SEXP res, bool IsArray, int nRows,
                    bool IsNamed, const std::vector<int> &vGrpSizes,
                    const std::vector<double> &mySample,
                    const std::vector<mpz_class> &myBigSamp,
                    bool IsSample,  bool IsGmp, int r, int n);
