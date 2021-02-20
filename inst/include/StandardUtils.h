#ifndef STANDARD_UTILS_H
#define STANDARD_UTILS_H

#include "CleanConvert.h"

std::vector<int> zUpdateIndex(SEXP x, SEXP y, int m);

SEXP CopyRv(const SEXP &Rv, const std::vector<int> &vInt,
            const std::vector<double> &vNum, VecType myType, bool IsFactor = false);

void SetRandomSample(SEXP RindexVec, SEXP RNumSamp, std::size_t &sampSize,
                     bool IsGmp, double computedRows, std::vector<double> &mySample,
                     Rcpp::Function baseSample);
    
#endif
