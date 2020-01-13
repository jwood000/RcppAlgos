#include "GmpDependUtils.h"

// [[Rcpp::export]]
SEXP CombinatoricsCount(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs, bool IsComb) {
    
    int n, m = 0, lenFreqs = 0;
    VecType myType = VecType::Integer;
    bool IsMult = false;
    
    std::vector<double> vNum;
    std::vector<int> vInt, myReps, freqs;
    
    bool IsRep = CleanConvert::convertLogical(RisRep, "repetition");
    SetType(myType, Rv);
    SetValues(myType, vInt, vNum, n, Rv);
    SetFreqsAndM(RFreqs, IsMult, myReps, IsRep, lenFreqs, freqs, Rm, n, m);
    
    const double computedRows = GetComputedRows(IsMult, IsComb, IsRep, n,
                                                m, Rm, lenFreqs, freqs, myReps);
    
    const bool IsGmp = (computedRows > Significand53);
    mpz_t computedRowsMpz;
    mpz_init(computedRowsMpz);
    
    if (IsGmp) {
        GetComputedRowMpz(computedRowsMpz, IsMult, 
                          IsComb, IsRep, n, m, Rm, freqs, myReps);
    }
    
    return GetCount(IsGmp, computedRowsMpz, computedRows);
}
