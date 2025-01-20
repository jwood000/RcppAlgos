#pragma once

#include "ComboGroups/ComboGroupsTemplate.h"
#include "ComboGroups/GetComboGroups.h"
#include "ClassUtils/ComboClass.h"

class ComboGroupsClass : public Combo {
private:

    cpp11::integers dim;
    cpp11::writable::list dimNames;
    cpp11::writable::strings myNames;

    const std::unique_ptr<ComboGroupsTemplate> CmbGrp;

    nextGrpFunc nextCmbGrp;
    nthFuncDbl nthCmbGrp;
    nthFuncGmp nthCmbGrpGmp;
    finalTouchFunc FinalTouch;

    std::string grpSizeDesc;
    bool IsArray;
    int rDisp; // This will differ in the General case when OneGrp = true
    int r;     // Number of groups

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
    SEXP nextIter();
    SEXP nextNumIters(SEXP RNum);
    SEXP nextGather();
    SEXP currIter();
    SEXP randomAccess(SEXP RindexVec);
    SEXP front();
    SEXP back();
    SEXP summary();
};
