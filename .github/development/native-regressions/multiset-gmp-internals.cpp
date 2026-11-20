#include "Partitions/MultisetCountClass.h"
#include "Partitions/NthPartition.h"
#include "Partitions/RankPartition.h"
#include "Partitions/PartitionsCountMultiset.h"
#include <numeric>
#include <stdexcept>
#include <string>
#include "check.h"

// changing reps preserve stored reps and both counting paths
[[cpp11::register]]
void check_multiset_changing_reps() {
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
    NATIVE_CHECK(counts == 58032);
}

// GMP ranking and unranking agree with double on small counts
[[cpp11::register]]
void check_multiset_gmp_rank_unrank() {
    const auto nth = GetNthPartsFunc(PartitionType::Multiset, true);
    NATIVE_CHECK(nth == &nthPartsMultisetGmp);
    const auto rankGmp = GetRankPartsFunc(PartitionType::Multiset, true);
    const auto rankDbl = GetRankPartsFunc(PartitionType::Multiset, false);
    NATIVE_CHECK(rankGmp == &rankPartsMultisetGmp);
    NATIVE_CHECK(rankDbl == &rankPartsMultiset);
    if (!nth || !rankGmp || !rankDbl) return;
    std::size_t checked = 0;

    // Public R tests independently enumerate these small cases. Direct
    // dispatch here forces GMP entry points despite their small counts.
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
                        auto actual = nth(
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
                        double dblRank = -1;
                        mpz_class mpzRank = -1;
                        rankGmp(actual.begin(), target, width, cap, 0,
                                dblRank, mpzRank, reps);
                        rankDbl(actual.begin(), target, width, cap, 0,
                                dblRank, mpzRank, reps);
                        if (mpzRank != rank || dblRank != rank) {
                            throw std::runtime_error(
                                "Multiset ranking mismatch: cap=" +
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
    NATIVE_CHECK(checked == 45659);
}
