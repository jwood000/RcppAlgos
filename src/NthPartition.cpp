// #inCclude "Partitions/BigPartsCount.h"
// #include "Partitions/PartsCount.h"
//
// using nthPartsPtr = std::vector<int> (*const)(int n, int r,
//                                       double dblIdx, mpz_t mpzIdx);
//
// std::vector<int> nthPartsDistinct(int n, int r, double dblIdx, mpz_t mpzIdx) {
//
//     double index1 = dblIdx, index2 = dblIdx;
//     std::vector<int> res(r);
//     double temp = nChooseK(n - 1, r - 1);
//
//     for (int k = 0, j = 0, n1 = n - 1,
//          r1 = r - 1; k < r; ++k, --n1, --r1, ++j) {
//         double test = temp;
//
//         for (int s = n1 - r1; test <= index1;
//                     --n1, ++j, --s, test += temp) {
//             index2 -= temp;
//             temp *= s;
//             temp /= n1;
//         }
//
//         temp *= r1;
//         temp /= n1;
//         res[k] = j;
//         index1 = index2;
//     }
//
//     return res;
// }
//
// std::vector<int> nthPartsRep(int n, int r, double dblIdx, mpz_t mpzIdx) {
//
//     double index1 = dblIdx, index2 = dblIdx;
//     std::vector<int> res(r);
//     double temp = NumCombsWithRep(n, r - 1);
//
//     for (int k = 0, j = 0, n1 = n, r1 = r - 1; k < r; ++k, --r1) {
//         double test = temp;
//
//         for (; test <= index1; --n1, ++j, test += temp) {
//             index2 -= temp;
//             temp *= (n1 - 1);
//             temp /= (n1 + r1 - 1);
//         }
//
//         temp *= r1;
//         temp /= (n1 + r1 - 1);
//         res[k] = j;
//         index1 = index2;
//     }
//
//     return res;
// }
//
// std::vector<int> nthPartsDistinctGmp(int n, int r,
//                                      double dblIdx, mpz_t mpzIdx) {
//
//     mpz_t test, temp, index1, index2;
//     mpz_init(test); mpz_init(temp);
//     mpz_init(index1); mpz_init(index2);
//     mpz_set(index1, mpzIdx);
//     mpz_set(index2, mpzIdx);
//
//     std::vector<int> res(r);
//     nChooseKGmp(temp, n - 1, r - 1);
//
//     for (int k = 0, j = 0, n1 = n - 1, r1 = r - 1;
//          k < r; ++k, --n1, --r1, ++j) {
//         mpz_set(test, temp);
//
//         for (int s = n1 - r1; mpz_cmp(test, index1) <= 0; --s, ++j, --n1) {
//             mpz_sub(index2, index2, temp);
//             mpz_mul_ui(temp, temp, s);
//             mpz_divexact_ui(temp, temp, n1);
//             mpz_add(test, test, temp);
//         }
//
//         mpz_mul_ui(temp, temp, r1);
//         if (n1 > 0) mpz_divexact_ui(temp, temp, n1);
//         res[k] = j;
//         mpz_set(index1, index2);
//     }
//
//     mpz_clear(index1); mpz_clear(index2);
//     mpz_clear(temp); mpz_clear(test);
//     return res;
// }
//
// std::vector<int> nthPartsRepGmp(int n, int r,
//                                 double dblIdx, mpz_t mpzIdx) {
//
//     mpz_t test, temp, index1, index2;
//     mpz_init(test); mpz_init(temp);
//     mpz_init(index1); mpz_init(index2);
//     mpz_set(index1, mpzIdx);
//     mpz_set(index2, mpzIdx);
//
//     std::vector<int> res(r);
//     NumCombsWithRepGmp(temp, n, r - 1);
//
//     for (int k = 0, j = 0, n1 = n, r1 = r - 1; k < r; ++k, --r1) {
//         mpz_set(test, temp);
//
//         for (; mpz_cmp(test, index1) <= 0; ++j, --n1) {
//             mpz_sub(index2, index2, temp);
//             mpz_mul_ui(temp, temp, n1 - 1);
//             mpz_divexact_ui(temp, temp, n1 + r1 - 1);
//             mpz_add(test, test, temp);
//         }
//
//         mpz_mul_ui(temp, temp, r1);
//         if ((n1 + r1) > 1) mpz_divexact_ui(temp, temp, n1 + r1 - 1);
//         res[k] = j;
//         mpz_set(index1, index2);
//     }
//
//     mpz_clear(index1); mpz_clear(index2);
//     mpz_clear(temp); mpz_clear(test);
//     return res;
// }
//
// nthPartsPtr GetNthPartsFunc(bool IsRep, bool IsGmp) {
//
//     if (IsGmp) {
//         if (IsRep) {
//             return(nthPartsPtr(nthPartsRepGmp));
//         } else {
//             return(nthPartsPtr(nthPartsDistinctGmp));
//         }
//     } else {
//         if (IsRep) {
//             return(nthPartsPtr(nthPartsRep));
//         } else {
//             return(nthPartsPtr(nthPartsDistinct));
//         }
//     }
// }
