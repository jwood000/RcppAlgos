#include "Partitions/NthPartition.h"
#include "Partitions/NthComposition.h"
#include "Partitions/PartitionsUtils.h"
#include "cpp11/R.hpp"
#include "cpp11/protect.hpp"

namespace {
    std::vector<int> nthLengthOne(
        int n, int m, int cap, int k, double dblIdx,
        const mpz_class &mpzIdx, const std::vector<int> & /* Reps */
    ) {

        std::vector<int> res(m, n - 1);
        return res;
    }

    std::vector<int> EmptyReturn(
        int n, int m, int cap, int k, double dblIdx,
        const mpz_class &mpzIdx, const std::vector<int> & /* Reps */
    ) {
        std::vector<int> res(m);
        return(res);
    }
} // namespace

nthPartsPtr GetNthPartsFunc(PartitionType ptype, bool IsGmp) {

    if (IsGmp) {
        switch (ptype) {
            case PartitionType::LengthOne: {
                return nthLengthOne;
            } case PartitionType::DistinctCapped: {
                return nthPartsDistinctCapGmp;
            } case PartitionType::DistinctCappedMZ: {
                return nthPartsDistinctCapMZGmp;
            } case PartitionType::DistinctNoZero: {
                return nthPartsDistinctLenGmp;
            } case PartitionType::DistinctOneZero: {
                return nthPartsDistinctLenGmp;
            } case PartitionType::DistinctMZ: {
                return nthPartsDistinctMZGmp;
            } case PartitionType::DistinctStdAll: {
                return nthPartsDistinctMZGmp;
            } case PartitionType::RepCapped: {
                return nthPartsRepCapGmp;
            } case PartitionType::RepNoZero: {
                return nthPartsRepLenGmp;
            } case PartitionType::RepShort: {
                return nthPartsRepShortGmp;
            } case PartitionType::RepStdAll: {
                return nthPartsRepGmp;
            } case PartitionType::Multiset: {
                return nthPartsMultisetGmp;
            } case PartitionType::CompRepNoZero: {
                return nthCompsRepGmp;
            } case PartitionType::CompRepCapped: {
                return nthCompsRepCappedGmp;
            } case PartitionType::CompRepWeak: {
                return nthCompsRepGmp;
            } case PartitionType::CompRepWeakCap: {
                return nthCompsRepCappedGmp;
            } case PartitionType::CompRepZero: {
                return nthCompsRepZeroGmp;
            } case PartitionType::CompRepCapZero: {
                return nthCompsRepZeroCapGmp;
            } case PartitionType::CompDistinctNoZero: {
                return nthCompsDistinctGmp;
            } case PartitionType::CompDistinctCapped: {
                return nthCompsDistinctGmp;
            } case PartitionType::CompDistinctWeak: {
                return nthCompsDistinctGmp;
            } case PartitionType::CompDistinctCapWeak: {
                return nthCompsDistinctGmp;
            } case PartitionType::CompDistinctZero: {
                return nthCompsDistinctMZGmp;
            } case PartitionType::CompDistinctCapMZ: {
                return nthCompsDistinctMZGmp;
            } case PartitionType::CompDistinctMZWeak: {
                return nthCompsDistinctMZWeakGmp;
            } case PartitionType::CompDistinctCapMZWeak: {
                return nthCompsDistinctMZWeakGmp;
            } case PartitionType::NoSolution: {
                return EmptyReturn;
            } default : {
                return nullptr;
            }
        }
    } else {
        switch (ptype) {
            case PartitionType::LengthOne: {
                return nthLengthOne;
            } case PartitionType::DistinctCapped: {
                return nthPartsDistinctCap;
            } case PartitionType::DistinctCappedMZ: {
                return nthPartsDistinctCapMZ;
            } case PartitionType::DistinctNoZero: {
                return nthPartsDistinctLen;
            } case PartitionType::DistinctOneZero: {
                return nthPartsDistinctLen;
            } case PartitionType::DistinctMZ: {
                return nthPartsDistinctMZ;
            } case PartitionType::DistinctStdAll: {
                return nthPartsDistinctMZ;
            } case PartitionType::RepCapped: {
                return nthPartsRepCap;
            } case PartitionType::RepNoZero: {
                return nthPartsRepLen;
            } case PartitionType::RepShort: {
                return nthPartsRepShort;
            } case PartitionType::RepStdAll: {
                return nthPartsRep;
            } case PartitionType::Multiset: {
                return nthPartsMultiset;
            } case PartitionType::CompRepNoZero: {
                return nthCompsRep;
            } case PartitionType::CompRepCapped: {
                return nthCompsRepCapped;
            } case PartitionType::CompRepWeak: {
                return nthCompsRep;
            } case PartitionType::CompRepWeakCap: {
                return nthCompsRepCapped;
            } case PartitionType::CompRepZero: {
                return nthCompsRepZero;
            } case PartitionType::CompRepCapZero: {
                return nthCompsRepZeroCap;
            } case PartitionType::CompDistinctNoZero: {
                return nthCompsDistinct;
            } case PartitionType::CompDistinctCapped: {
                return nthCompsDistinct;
            } case PartitionType::CompDistinctWeak: {
                return nthCompsDistinct;
            } case PartitionType::CompDistinctCapWeak: {
                return nthCompsDistinct;
            } case PartitionType::CompDistinctZero: {
                return nthCompsDistinctMZ;
            } case PartitionType::CompDistinctCapMZ: {
                return nthCompsDistinctMZ;
            } case PartitionType::CompDistinctMZWeak: {
                return nthCompsDistinctMZWeak;
            } case PartitionType::CompDistinctCapMZWeak: {
                return nthCompsDistinctMZWeak;
            } case PartitionType::NoSolution: {
                return EmptyReturn;
            } default : {
                return nullptr;
            }
        }
    }
}

nthPartsPtr GetNthPartsFuncOrStop(PartitionType ptype, bool IsGmp) {
    if (auto res = GetNthPartsFunc(ptype, IsGmp)) {
        return res;
    }

    cpp11::stop(
        "No algorithm available for PartitionType = " + GetPTypeName(ptype)
    );

    return nullptr;
}
