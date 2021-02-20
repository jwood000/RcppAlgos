// // #include "CleanConvert.h"
// 
// void InitialSetup(std::vector<int> &vInt, std::vector<double> &vNum,
//                   int &n, int &m, bool &IsRep, bool &IsMult, bool &IsGmp,
//                   bool &generalRet, std::vector<int> &freqs,
//                   std::vector<int> &startZ, std::vector<int> &myReps, 
//                   double &lower, mpz_t *const lowerMpz, int &nRows,
//                   int &nThreads, bool &Parallel, int &phaseOne,
//                   VecType &myType, SEXP Rm, SEXP Rv, bool IsComb) {
//     
//     Parallel = CleanConvert::convertLogical(Rparallel, "Parallel");
//     IsRep = CleanConvert::convertLogical(RisRep, "repetition");
//     
//     SetType(myType, Rv);
//     SetValues(myType, vInt, vNum, n, Rv);
//     SetFreqsAndM(RFreqs, IsMult, myReps, IsRep, freqs, Rm, n, m);
//     
//     const double computedRows = GetCombCount(n, m, IsMult, IsRep, freqs, myReps);
//     const bool IsGmp = (computedRows > Significand53);
//     
//     mpz_t computedRowsMpz;
//     mpz_init(computedRowsMpz);
//     
//     if (IsGmp) {
//         GetCombCountMpz(computedRowsMpz, n, m, IsMult, IsRep, myReps);
//     }
//     
//     double lower = 0;
//     double upper = 0;
//     
//     bool bLower = false;
//     bool bUpper = false;
//     
//     auto upperMpz = FromCpp14::make_unique<mpz_t[]>(1);
//     mpz_init(upperMpz[0]);
//     
//     SetBounds(Rlow, Rhigh, IsGmp, bLower, bUpper, lower, upper,
//               lowerMpz.get(), upperMpz.get(), computedRowsMpz, computedRows);
//     
//     std::vector<int> startZ(m);
//     SetCombStartZ(myReps, freqs, startZ, n, m,
//                   lower, lowerMpz[0], IsRep, IsMult, IsGmp);
//     
//     double userNumRows = 0;
//     SetNumResults(IsGmp, bLower, bUpper, false, upperMpz.get(), lowerMpz.get(),
//                   lower, upper, computedRows, computedRowsMpz, nRows, userNumRows);
//     
//     const int limit = 20000;
//     SetThreads(Parallel, maxThreads, nRows, myType, nThreads, RNumThreads, limit);
// }