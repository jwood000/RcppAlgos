#pragma once

#include "cpp11/sexp.hpp"

#include "Constants.h"
#include <gmpxx.h>
#include <vector>
#include <string>

namespace CppConvert {

    void convertMPZVector(SEXP input, std::vector<mpz_class> &myVec,
                          std::size_t vecSize, const std::string &nameOfObject,
                          bool negPoss = false);

    void convertMpzClass(SEXP input, mpz_class &result,
                         const std::string &nameOfObject,
                         bool negPoss = false);
}
