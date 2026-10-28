#include "Partitions/MultisetCountClass.h"
#include "Partitions/NthPartition.h"
#include "Partitions/PartitionsCountMultiset.h"
#include <numeric>
#include <stdexcept>
#include <string>
#include <testthat.h>

context("GMP multiset internals") {
    test_that("changing reps preserve stored reps and both counting paths") {
        mpz_class large = 1;
        large <<= 100;
        std::size_t counts = 0;

        for (int cap = 1; cap <= 5; ++cap) {
            const std::vector<int> stored(cap, 3);
            PartsMultiset counter(stored);
            counter.SetArrSize(PartitionType::Multiset, 30, 6);
            counter.InitializeMpz();
            std::vector<int> allowed(cap);
            std::iota(allowed.begin(), allowed.end(), 1);
            int combinations = 1;
            for (int i = 0; i < cap; ++i) combinations *= 3;

            for (int width = 1; width <= 6; ++width) {
                for (int combination = 0; combination < combinations;
                     ++combination) {
                    int digits = combination;
                    std::vector<int> reps(cap);
                    for (int& f : reps) {
                        f = 1 + digits % 3;
                        digits /= 3;
                    }
                    for (int target = width; target <= cap * width; ++target) {
                        const std::string inputs = " cap=" +
                            std::to_string(cap) + " width=" +
                            std::to_string(width) + " target=" +
                            std::to_string(target) + " reps-code=" +
                            std::to_string(combination);
                        const double expected = CountPartsMultiset(
                            target, width, allowed, reps
                        );
                        for (bool forceGmp : {false, true}) {
                            mpz_class count = forceGmp ? large : mpz_class(0);
                            counter.GetCount(
                                count, target, width, allowed, reps
                            );
                            if (count != expected) {
                                throw std::runtime_error(
                                    "Changing reps mismatch; forceGmp=" +
                                    std::to_string(forceGmp) + inputs
                                );
                            }
                            ++counts;
                        }

                        // The inherited interface still uses stored Reps.
                        mpz_class count = large;
                        CountClass& base = counter;
                        base.GetCount(count, target, width, allowed);
                        if (count != CountPartsMultiset(
                                target, width, allowed, stored)) {
                            throw std::runtime_error(
                                "Stored reps changed" + inputs
                            );
                        }
                    }
                }
            }
        }
        expect_true(counts == 58032);
    }

    test_that("GMP dispatch agrees with double unranking on small counts") {
        const auto nth = GetNthPartsFunc(PartitionType::Multiset, true);
        expect_true(nth == &nthPartsMultisetGmp);
        if (!nth) return;
        std::size_t checked = 0;

        // Public R tests independently enumerate these small cases. Direct
        // dispatch here forces GMP unranking despite their small counts.
        for (int cap = 1; cap <= 5; ++cap) {
            std::vector<int> allowed(cap);
            std::iota(allowed.begin(), allowed.end(), 1);
            int combinations = 1;
            for (int i = 0; i < cap; ++i) combinations *= 3;

            for (int width = 1; width <= 6; ++width) {
                for (int combination = 0; combination < combinations;
                     ++combination) {
                    int digits = combination;
                    std::vector<int> reps(cap);
                    for (int& f : reps) {
                        f = 1 + digits % 3;
                        digits /= 3;
                    }
                    for (int target = width; target <= cap * width; ++target) {
                        const int total = CountPartsMultiset(
                            target, width, allowed, reps
                        );
                        for (int rank = 0; rank < total; ++rank) {
                            const mpz_class index(rank);
                            const auto actual = nth(
                                target, width, cap, 0, 0, index, reps
                            );
                            const auto expected = nthPartsMultiset(
                                target, width, cap, 0, rank, index, reps
                            );
                            if (actual != expected || index != rank) {
                                throw std::runtime_error(
                                    "GMP unranking mismatch: cap=" +
                                    std::to_string(cap) + " width=" +
                                    std::to_string(width) + " target=" +
                                    std::to_string(target) + " reps-code=" +
                                    std::to_string(combination) + " rank=" +
                                    std::to_string(rank)
                                );
                            }
                            ++checked;
                        }
                    }
                }
            }
        }
        expect_true(checked == 45659);
    }
}
