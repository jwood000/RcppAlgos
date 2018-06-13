#ifndef GET_FACS_UTILS_H
#define GET_FACS_UTILS_H

const double Significand53 = 9007199254740991.0;

/**
 * Get Prime Factorization
 * t: number to factorize
 */

template <typename typeReturn>
void getPrimefactors (int64_t &t, std::vector<typeReturn> &factors);

#endif
