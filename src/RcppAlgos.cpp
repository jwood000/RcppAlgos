#include "NumbersUtils/PrimeSieveCount.h"
#include "NumbersUtils/MotleyPrimes.h"
#include "NumbersUtils/DivNumSieve.h"
#include "NumbersUtils/PollardRho.h"
#include "Sample/SamplePartitions.h"
#include "ClassUtils/ExposeClass.h"
#include "Sample/SampleCombPerm.h"
#include "CombinatoricsCnstrt.h"
#include "CombinatoricsCount.h"
#include "CombinatoricsApply.h"
#include "CartesianContainer.h"
#include "CombinatoricsMain.h"
#include "GetClassValues.h"
#include "ComboGroups.h"
#include "CheckReturn.h"
#include "NumThreads.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef CallEntries[] = {
    CALLDEF(CheckReturn, 6),
    CALLDEF(CheckConstrndCpp, 3),
    CALLDEF(CombinatoricsCount, 5),
    CALLDEF(PartitionsCount, 10),
    CALLDEF(CombinatoricsApply, 10),
    CALLDEF(CombinatoricsStndrd, 10),
    CALLDEF(CombinatoricsCnstrt, 15),
    CALLDEF(SamplePartitions, 16),
    CALLDEF(SampleCombPerm, 16),
    CALLDEF(ComboGroupsCountCpp, 2),
    CALLDEF(ComboGroupsCpp, 15),
    CALLDEF(ComboGridCpp, 2),
    CALLDEF(PrimeSieveCpp, 5),
    CALLDEF(PrimeCountCpp, 3),
    CALLDEF(MotleyContainer, 6),
    CALLDEF(DivNumSieveCpp, 6),
    CALLDEF(PollardRhoContainer, 6),
    CALLDEF(cpp11GetNumThreads, 0),
    CALLDEF(GetClassVals, 9),
    CALLDEF(CombClassNew, 14),
    CALLDEF(StartOverGlue, 1),
    CALLDEF(NextCombGlue, 1),
    CALLDEF(NextNumCombGlue, 2),
    CALLDEF(NextGatherGlue, 1),
    CALLDEF(PrevCombGlue, 1),
    CALLDEF(PrevNumCombGlue, 2),
    CALLDEF(PrevGatherGlue, 1),
    CALLDEF(CurrCombGlue, 1),
    CALLDEF(SourceVectorGlue, 1),
    CALLDEF(RandomAccessGlue, 2),
    CALLDEF(FrontGlue, 1),
    CALLDEF(BackGlue, 1),
    CALLDEF(SummaryGlue, 1),
    {NULL, NULL, 0}
};

extern "C"
void R_init_RcppAlgos(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, (Rboolean) FALSE);
}
