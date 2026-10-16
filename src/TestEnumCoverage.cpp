#include "cpp11/R.hpp"
#include "cpp11/protect.hpp"

#include "Partitions/MultisetCountClass.h"
#include "Partitions/PartitionsUtils.h"
#include "Partitions/PartitionsCount.h"
#include "Partitions/RankPartition.h"
#include "Partitions/NthPartition.h"
#include <memory>

[[cpp11::register]]
bool testRankCoverage() {
    for (int i = 0; i < static_cast<int>(PartitionType::NumTypes); ++i) {
        const auto ptype = static_cast<PartitionType>(i);

        if (!HasRankAlgo(ptype)) continue;

        for (const bool isGmp : {false, true}) {
            const auto nth = GetNthPartsFunc(ptype, isGmp);

            if (!nth) {
                cpp11::stop(
                    "Missing unranking algorithm for %s (%s)",
                    GetPTypeName(ptype).c_str(),
                    isGmp ? "GMP" : "double"
                );
            }

            // This factory throws if no implementation exists.
            const auto rank = GetRankPartsFunc(ptype, isGmp);

            if (!rank) {
                cpp11::stop(
                    "Missing ranking algorithm for %s (%s)",
                    GetPTypeName(ptype).c_str(),
                    isGmp ? "GMP" : "double"
                );
            }
        }
    }

    return true;
}


[[cpp11::register]]
bool testCountCoverage() {
    for (int i = 0; i < static_cast<int>(PartitionType::NumTypes); ++i) {
        const auto ptype = static_cast<PartitionType>(i);

        if (!HasCountAlgo(ptype)) continue;

        // Handled directly by PartitionsCount().
        if (ptype == PartitionType::LengthOne) continue;

        std::unique_ptr<CountClass> Counter = MakeCount(ptype);

        if (!Counter) {
            const std::vector<int> Reps = {1, 2, 3};
            Counter = MakeMultisetCount(ptype, Reps);
        }

        if (!Counter) {
            cpp11::stop(
                "Missing counting algorithm for %s",
                GetPTypeName(ptype).c_str()
            );
        }
    }

    return true;
}
