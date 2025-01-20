#pragma once

#include "cpp11/strings.hpp"

#include "ComboGroups/GroupHelperClass.h"
#include "Combinations/NthCombination.h"
#include "Combinations/BigComboCount.h"
#include "Combinations/ComboCount.h"
#include "SetUpUtils.h"
#include <numeric>
#include <memory>
#include <cmath>

typedef std::function<std::vector<int>(const mpz_class &)> nthFuncGmp;
typedef std::function<std::vector<int>(double)>            nthFuncDbl;
typedef std::function<bool(std::vector<int>&)>             nextGrpFunc;

typedef std::function<void(
    SEXP, bool, int, bool, const std::vector<double>&,
    const std::vector<mpz_class>&, bool
)> finalTouchFunc;

struct CmbGrpClsFuncs {
    const nthFuncDbl nthDbl;
    const nthFuncGmp nthGmp;
    const nextGrpFunc next;
    const finalTouchFunc finishing;
};

class ComboGroupsTemplate {
protected:

    std::string GroupType;

    bool OneGrp; // Used only for General case, but we need to be able to
                 // to access it, so that we can properly name the output
                 // when using iterables (e.g. nextIter method)
    const int n; // Size of vector which is also the size of z (i.e. z.size())
    const int r; // Number of groups

    const int idx1;
    const int idx2;
    const int curr_bnd;

    bool IsGmp;
    double computedRows;
    mpz_class computedRowsMpz;

public:

    virtual ~ComboGroupsTemplate() = default;
    ComboGroupsTemplate(int n_, int numGroups, int i1, int i2, int bnd);

    virtual bool nextComboGroup(std::vector<int> &z) = 0;
    virtual double numGroupCombs() = 0;
    virtual mpz_class numGroupCombsGmp() = 0;
    virtual std::vector<int> nthComboGroup(double myIndex) = 0;
    virtual std::vector<int> nthComboGroupGmp(const mpz_class &lowerMpz) = 0;
    virtual std::vector<int> GetGroupSizes() = 0;

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
    std::string GetType() const {return GroupType;}
    bool GetOneGrp() const {return OneGrp;}

    SEXP GetCount() const {
        return CppConvert::GetCount(IsGmp, computedRowsMpz, computedRows);
    };
};

std::unique_ptr<ComboGroupsTemplate> MakeComboGroup(
    const std::vector<int> &vGrpSize, GroupHelper *const MyGrp,
    int i1, int i2, int bnd, int grpSize, bool IsGen, bool IsUni
);

std::unique_ptr<ComboGroupsTemplate> GroupPrep(SEXP Rv, SEXP RNumGroups,
                                               SEXP RGrpSize, int n);

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
