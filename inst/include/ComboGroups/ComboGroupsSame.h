#pragma once

#include "ComboGroups/ComboGroupsTemplate.h"

class ComboGroupsSame : public ComboGroupsTemplate {
private:

    const int grpSize;

public:

    ComboGroupsSame(int n_, int numGroups, int i1, int i2, int bnd, int size);

    bool nextComboGroup(std::vector<int> &z);
    double numGroupCombs();
    mpz_class numGroupCombsGmp();
    std::vector<int> nthComboGroup(double myIndex);
    std::vector<int> nthComboGroupGmp(const mpz_class &lowerMpz);

    std::vector<int> GetGroupSizes() {
        std::vector<int> grp(1, grpSize);
        return grp;
    }

    void FinalTouch(
        SEXP res, bool IsArray, int nRows, bool IsNamed,
        const std::vector<double> &mySample,
        const std::vector<mpz_class> &myBigSamp, bool IsSample
    );
};
