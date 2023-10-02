#pragma once

#include "ComboGroups/ComboGroupsTemplate.h"

class ComboGroupsUnique : public ComboGroupsTemplate {
private:

    const std::vector<int> grp;

public:

    ComboGroupsUnique(int n_, int numGroups, int i1, int i2,
                      int bnd, const std::vector<int> &grp_);

    bool nextComboGroup(std::vector<int> &z);
    double numGroupCombs();
    mpz_class numGroupCombsGmp();
    std::vector<int> nthComboGroup(double myIndex);
    std::vector<int> nthComboGroupGmp(const mpz_class &lowerMpz);
    std::vector<int> GetGroupSizes() {return grp;}

    void FinalTouch(
        SEXP res, bool IsArray, int nRows, bool IsNamed,
        const std::vector<double> &mySample,
        const std::vector<mpz_class> &myBigSamp, bool IsSample
    );
};
