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

        for (const bool isGmp : {false, true}) {
            if (ptype == PartitionType::NoSolution) continue;

            const auto nth = GetNthPartsFunc(ptype, isGmp);
            const auto rank = GetRankPartsFunc(ptype, isGmp);
            const bool advertised = HasRankAlgo(ptype);

            if (static_cast<bool>(nth) != advertised) {
                cpp11::stop(
                    "Unranking coverage mismatch for %s (%s)",
                    GetPTypeName(ptype).c_str(),
                    isGmp ? "GMP" : "double"
                );
            }

            if (static_cast<bool>(rank) != advertised) {
                cpp11::stop(
                    "Ranking coverage mismatch for %s (%s)",
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

        // Handled directly by PartitionsCount().
        if (ptype == PartitionType::LengthOne) continue;

        const bool advertised = HasCountAlgo(ptype);
        std::unique_ptr<CountClass> Counter = MakeCount(ptype);

        if (!Counter) {
            const std::vector<int> Reps = {1, 2, 3};
            Counter = MakeMultisetCount(ptype, Reps);
        }

        if (static_cast<bool>(Counter) != advertised) {
            cpp11::stop(
                "Counting coverage mismatch for %s",
                GetPTypeName(ptype).c_str()
            );
        }
    }

    return true;
}
