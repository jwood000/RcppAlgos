#include "NumbersUtils/PollardRhoUtils.h"
#include <limits>
#include <stdexcept>
#include <testthat.h>

context("Native primality bounds") {
    test_that("native primality accepts its inclusive upper bound") {
        // INT_MAX is prime; 46337 is prime and exceeds the trial-division
        // table, so both inputs satisfy IsPrime's small-factor precondition.
        const std::int64_t largest = std::numeric_limits<int>::max();
        expect_true(IsPrime(largest) == 1);
        expect_true(IsPrime(std::int64_t(46337) * 46337) == 0);
    }

    test_that("native primality rejects inputs beyond its supported bound") {
        const std::int64_t firstUnsupported =
            std::int64_t(std::numeric_limits<int>::max()) + 1;
        const std::int64_t inputs[] = {
            firstUnsupported, std::numeric_limits<std::int64_t>::max()
        };

        for (const std::int64_t input : inputs) {
            bool rejected = false;
            try {
                IsPrime(input);
            } catch (const std::out_of_range&) {
                rejected = true;
            }
            expect_true(rejected);
        }
    }
}
