#include "Cartesian/GetProduct.h"

[[cpp11::register]]
SEXP ExpandGridCpp(
    cpp11::list RList, SEXP Rlow, SEXP Rhigh, SEXP RNumThreads,
    SEXP RmaxThreads, bool IsSample, SEXP RindexVec, SEXP RmySeed,
    SEXP RNumSamp, SEXP baseSample, SEXP RNamed, SEXP myEnv, bool Force_DF
) {

    bool IsNamed  = (IsSample) ?
        CppConvert::convertFlag(RNamed, "namedSample") : false;

    const int nCols = RList.size();
    std::vector<std::vector<int>> myVec(nCols);
    std::vector<int> typeCheck(N_TYPES, 0);

    std::vector<int> IsFactor(nCols);
    std::vector<int> lenGrps(nCols);
    CartesianInitialPrep(RList, IsFactor, lenGrps, nCols);

    const int sumLength = std::accumulate(
        lenGrps.begin(), lenGrps.end(), 0
    );

    cpp11::writable::strings charVec(sumLength);
    std::vector<Rcomplex> cmplxVec(sumLength);
    std::vector<Rbyte> rawVec(sumLength);
    std::vector<double> dblVec(sumLength);
    std::vector<int> intVec(sumLength);
    std::vector<int> boolVec(sumLength);

    VecType myType = VecType::Integer;
    bool IsDF = true;

    ProductPrepare(
        RList, IsFactor, lenGrps, myVec, charVec, cmplxVec, rawVec,
        dblVec, intVec, boolVec, typeCheck, myType, nCols, IsDF
    );

    IsDF = IsDF || Force_DF;

    int nRows = 0;
    int nThreads = 1;
    int maxThreads = 1;

    CppConvert::convertPrimitive(RmaxThreads, maxThreads,
                                 VecType::Integer, "maxThreads");

    const double computedRows = CartesianCount(lenGrps);
    const bool IsGmp = IsSample ? computedRows > SampleLimit :
        computedRows > Significand53;

    mpz_class computedRowsMpz;

    if (IsGmp) {
        CartesianCountGmp(computedRowsMpz, lenGrps);
    }

    double lower = 0;
    double upper = 0;
    bool Parallel = false;

    bool bLower = false;
    bool bUpper = false;

    mpz_class lowerMpz;
    mpz_class upperMpz;

    if (!IsSample) {
        SetBounds(Rlow, Rhigh, IsGmp, bLower, bUpper, lower, upper,
                  lowerMpz, upperMpz, computedRowsMpz, computedRows);

        double userNumRows = 0;
        SetNumResults(IsGmp, bLower, bUpper, true, upperMpz,
                      lowerMpz, lower, upper, computedRows,
                      computedRowsMpz, nRows, userNumRows);
    }

    int sampSize;
    std::vector<double> mySample;

    if (IsSample) {
        SetRandomSample(RindexVec, RNumSamp, sampSize, IsGmp,
                        computedRows, mySample, baseSample, myEnv);
    }

    const int bigSampSize = (IsSample && IsGmp) ? sampSize : 1;
    std::vector<mpz_class> myBigSamp(bigSampSize);

    if (IsSample) {
        SetRandomSampleMpz(RindexVec, RmySeed, sampSize,
                           IsGmp, computedRowsMpz, myBigSamp);
    }

    const int numResults = (IsSample) ? sampSize : nRows;
    const int limit = (IsSample) ? 2 : 20000;

    SetThreads(Parallel, maxThreads, numResults,
               myType, nThreads, RNumThreads, limit);

    std::vector<int> startZ;
    GetStartProd(lenGrps, startZ, lowerMpz, lower, 0, IsGmp);

    const int maxLen = *std::max_element(lenGrps.begin(), lenGrps.end());
    std::vector<int> mat_idx(maxLen * nCols);

    // transposing myVec so that each row represents
    // each element of the given list
    for (int i = 0; i < nCols; ++i) {
        for (int j = 0; j < lenGrps[i]; ++j) {
            mat_idx[i + j * nCols] = myVec[i][j];
        }
    }

    // Transform lenGrps to be used in nextProduct
    for (auto &v_i: lenGrps) {
        v_i = nCols * (v_i - 1);
    }

    cpp11::sexp res = GetProduct(
        mat_idx, typeCheck, IsFactor, RList, intVec, dblVec, boolVec, cmplxVec,
        rawVec, charVec, lenGrps, startZ, mySample, myBigSamp, lower, lowerMpz,
        numResults, nCols, IsDF, nThreads, Parallel, IsGmp, IsSample
    );

    if (IsSample && IsNamed) {
        SetSampleNames(res, IsGmp, numResults, mySample,
                       myBigSamp, IsNamed, RList.names(), 1);
    } else {
        SetMatrixColnames(res, RList.names());
    }

    return res;
}
