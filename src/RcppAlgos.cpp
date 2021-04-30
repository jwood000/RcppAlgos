#include "CombinatoricsCnstrt.h"
#include "CombinatoricsCount.h"
#include "CombinatoricsApply.h"
#include "CombinatoricsMain.h"
#include "CheckReturn.h"
#include "NumThreads.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef CallEntries[] = {
    CALLDEF(CheckReturn, 6),
    CALLDEF(CombinatoricsCount, 5),
    CALLDEF(PartitionsCount, 9),
    CALLDEF(CombinatoricsApply, 10),
    CALLDEF(CombinatoricsStndrd, 10),
    CALLDEF(CombinatoricsCnstrt, 15),
    CALLDEF(cpp11GetNumThreads, 0),
    {NULL, NULL, 0}
};

extern "C"
void R_init_RcppAlgos(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, (Rboolean) FALSE);
}
