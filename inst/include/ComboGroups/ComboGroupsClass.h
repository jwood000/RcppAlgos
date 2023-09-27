#pragma once

#include "ComboGroups/ComboGroupsTemplate.h"
#include "ComboGroups/GetComboGroups.h"
#include "ClassUtils/ComboClass.h"

class ComboGroupsClass : public Combo {
private:

    cpp11::integers dim;
    cpp11::writable::list dimNames;
    cpp11::writable::strings myNames;

    bool IsArray;
    int r; // Number of groups
    const std::unique_ptr<ComboGroupsTemplate> CmbGrp;

    SEXP SingleReturn();
    SEXP GeneralReturn(int numResults);

public:

    ComboGroupsClass(
        SEXP Rv, int Rm, SEXP RcompRows, const std::vector<int> &bVec,
        const std::vector<int> &Rreps, const std::vector<int> &Rfreqs,
        const std::vector<int> &RvInt, const std::vector<double> &RvNum,
        VecType typePass, int RmaxThreads, SEXP RnumThreads, bool Rparallel,
        SEXP RNumGroups, SEXP RGrpSize, SEXP RRetType
    );

    void startOver();
    SEXP nextComb();
    SEXP nextNumCombs(SEXP RNum);
    SEXP nextGather();
    SEXP currComb();
    SEXP randomAccess(SEXP RindexVec);
    SEXP front();
    SEXP back();
    SEXP summary();
};
