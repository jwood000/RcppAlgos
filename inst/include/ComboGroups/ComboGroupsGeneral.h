#pragma once

#include "ComboGroups/ComboGroupsTemplate.h"

class ComboGroupsGeneral : public ComboGroupsTemplate {
private:

    int genGrps; // Could be different than r defined in the template
                 // class when there are groups of size 1.
    const GroupHelper MyGrp;
    std::vector<int> realGrps;

public:

    ComboGroupsGeneral(int n_, int numGroups, int i1, int i2,
                       int bnd, GroupHelper MyGrp_, bool OneGrp_);

    bool nextComboGroup(std::vector<int> &z);
    double numGroupCombs();
    mpz_class numGroupCombsGmp();
    std::vector<int> nthComboGroup(double myIndex);
    std::vector<int> nthComboGroupGmp(const mpz_class &lowerMpz);
    std::vector<int> GetGroupSizes() {return MyGrp.grp;}

    void FinalTouch(
        SEXP res, bool IsArray, int nRows, bool IsNamed,
        const std::vector<double> &mySample,
        const std::vector<mpz_class> &myBigSamp, bool IsSample
    );
};
