#include "CppConvert/GmpConvert.h"
#include "cpp11/function.hpp"
#include "cpp11/strings.hpp"
#include "check.h"

// negPoss preserves zero and signed bigz conversions
[[cpp11::register]]
void check_signed_gmp_conversions() {
    const cpp11::function asBigz = cpp11::package("gmp")["as.bigz"];
    const cpp11::writable::strings values = {
        "0", "-1", "9007199254740993", "-9007199254740993"
    };
    const cpp11::sexp input = asBigz(values);
    std::vector<mpz_class> actual(values.size());
    CppConvert::convertMPZVector(
        input, actual, actual.size(), "values", true
    );

    for (R_xlen_t i = 0; i < values.size(); ++i) {
        const std::string value = cpp11::r_string(values[i]);
        const mpz_class expected(value);
        const bool vectorMatches = actual[i] == expected;
        NATIVE_CHECK(vectorMatches);

        const cpp11::sexp scalarInput = asBigz(value);
        mpz_class scalar;
        CppConvert::convertMpzClass(scalarInput, scalar, "value", true);
        const bool scalarMatches = scalar == expected;
        NATIVE_CHECK(scalarMatches);
    }
}
