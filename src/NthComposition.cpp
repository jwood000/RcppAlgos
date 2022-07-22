// #include "Partitions/PartitionsCountMultiset.h"
// #include "Partitions/PartitionsCountDistinct.h"
// #include "Partitions/PartitionsCountRep.h"
// #include "Partitions/BigPartsCountRep.h"
// #include "Partitions/PartitionsTypes.h"
// #include "Partitions/PartitionsCount.h"
// #include <numeric>  // std::accumulate
//
// #include "cpp11/R.hpp"
// #include "cpp11/protect.hpp"
//
// // The variable k is strtLen
// using nthCompsPtr = std::vector<int> (*const)(int n, int m, bool includeZero,
//                                       double dblIdx, mpz_t mpzIdx);
//
// std::vector<int> nthCompsRepLen(int n, int m, bool includeZero,
//                                 double dblIdx, mpz_t mpzIdx) {
//
//     const int width = m;
//     const int max_n = n;
//
//     std::vector<int> res(width);
//     --n;
//     --m;
//
//     for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m) {
//         for (double temp = CountPartsPermRep(n, m);
//              temp <= dblIdx; ++j) {
//             n -= (m + 1);
//             dblIdx -= temp;
//             temp = CountPartsPermRep(n, m);
//         }
//
//         res[i] = j;
//     }
//
//     res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
//     return res;
// }
//
// std::vector<int> nthCompsRepShort(int n, int m, bool includeZero,
//                                   double dblIdx, mpz_t mpzIdx) {
//
//     return nthCompsRepLen(n, m, includeZero, dblIdx, mpzIdx);
// }
//
// std::vector<int> nthCompsRep(int n, int m, bool includeZero,
//                              double dblIdx, mpz_t mpzIdx) {
//
//     return nthCompsRepLen(n, m, includeZero, dblIdx, mpzIdx);
// }
//
// std::vector<int> nthCompsRepLenGmp(int n, int m, bool includeZero,
//                                    double dblIdx, mpz_t mpzIdx) {
//
//     const int width = m;
//     const int max_n = n;
//
//     std::vector<int> res(width);
//     --n;
//     --m;
//
//     mpz_t temp;
//     mpz_t index;
//
//     mpz_init(temp);
//     mpz_init(index);
//
//     mpz_set(index, mpzIdx);
//     const PartitionType ptype = PartitionType::RepShort;
//
//     std::unique_ptr<CountClass> myClass = MakeCount(ptype);
//     myClass->SetArrSize(ptype, n, m, cap);
//     myClass->InitializeMpz();
//
//     for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m) {
//         myClass->GetCount(temp, n, m, cap, k);
//
//         for (; mpz_cmp(temp, index) <= 0; ++j) {
//             n -= (m + 1);
//             mpz_sub(index, index, temp);
//             myClass->GetCount(temp, n, m, cap, k);
//         }
//
//         res[i] = j;
//     }
//
//     res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
//
//     mpz_clear(temp);
//     mpz_clear(index);
//     myClass->ClearMpz();
//     return res;
// }
//
// std::vector<int> nthCompsRepShortGmp(int n, int m, bool includeZero,
//                                      double dblIdx, mpz_t mpzIdx) {
//
//     return nthCompsRepLenGmp(n, m, includeZero, dblIdx, mpzIdx);
// }
//
// std::vector<int> nthCompsRepGmp(int n, int m, bool includeZero,
//                                 double dblIdx, mpz_t mpzIdx) {
//
//     return nthCompsRepLenGmp(n, m,includeZero, dblIdx, mpzIdx);
// }
//
// nthCompsPtr GetNthCompsFunc(PartitionType ptype, bool IsGmp) {
//
//     if (IsGmp) {
//         switch (ptype) {
//             case PartitionType::RepNoZero: {
//                 return(nthCompsPtr(nthCompsRepLenGmp));
//             } case PartitionType::RepShort : {
//                 return(nthCompsPtr(nthCompsRepShortGmp));
//             } case PartitionType::RepStdAll : {
//                 return(nthCompsPtr(nthCompsRepGmp));
//             } default : {
//                 cpp11::stop("No algorithm available");
//             }
//         }
//     } else {
//         switch (ptype) {
//             case PartitionType::RepNoZero: {
//                 return(nthCompsPtr(nthCompsRepLen));
//             } case PartitionType::RepShort: {
//                 return(nthCompsPtr(nthCompsRepShort));
//             } case PartitionType::RepStdAll: {
//                 return(nthCompsPtr(nthCompsRep));
//             }  default : {
//                 cpp11::stop("No algorithm available");
//             }
//         }
//     }
// }
