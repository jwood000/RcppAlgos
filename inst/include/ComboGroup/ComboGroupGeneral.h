#pragma once

#include "ComboGroup/ComboGroupClass.h"

class ComboGroupGeneral : public ComboGroup {
private:

    const Group MyGrp;
    const bool OneGrp;

public:

    ComboGroupGeneral(int n_, int numGroups, int i1, int i2,
                      int bnd, Group MyGrp_, bool OneGrp_);

    bool nextComboGroup(std::vector<int> &z);
    double numGroupCombs();
    mpz_class numGroupCombsGmp();
    std::vector<int> nthComboGroup(double myIndex);
    std::vector<int> nthComboGroupGmp(const mpz_class &lowerMpz);

    void FinalTouch(
        SEXP res, bool IsArray, int nRows, bool IsNamed,
        const std::vector<double> &mySample,
        const std::vector<mpz_class> &myBigSamp, bool IsSample
    );
};
