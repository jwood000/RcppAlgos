#include "Partitions/RankPartition.h"
#include "Partitions/RankComposition.h"

namespace {
    void rankLengthOne(
        std::vector<int>::iterator iter, int n, int m, int cap, int k,
        double &dblIdx, mpz_class &mpzIdx, const std::vector<int> & /* Reps */
    ) {
        dblIdx = 0;
        mpzIdx = 0;
    }
} // namespace

rankPartsPtr GetRankPartsFunc(PartitionType ptype, bool IsGmp) {

    if (IsGmp) {
        switch (ptype) {
            case PartitionType::LengthOne: {
                return(rankPartsPtr(rankLengthOne));
            } case PartitionType::DistinctCapped: {
                return(rankPartsPtr(rankPartsDistinctCapGmp));
            } case PartitionType::DistinctCappedMZ: {
                return(rankPartsPtr(rankPartsDistinctCapMZGmp));
            } case PartitionType::DistinctNoZero: {
                return(rankPartsPtr(rankPartsDistinctLenGmp));
            } case PartitionType::DistinctOneZero: {
                return(rankPartsPtr(rankPartsDistinctOneZeroGmp));
            } case PartitionType::DistinctMZ: {
                return(rankPartsPtr(rankPartsDistinctMZGmp));
            } case PartitionType::DistinctStdAll: {
                return(rankPartsPtr(rankPartsDistinctMZGmp));
            } case PartitionType::RepCapped: {
                return(rankPartsPtr(rankPartsRepCapGmp));
            } case PartitionType::RepNoZero: {
                return(rankPartsPtr(rankPartsRepLenGmp));
            } case PartitionType::RepShort: {
                return(rankPartsPtr(rankPartsRepShortGmp));
            } case PartitionType::RepStdAll: {
                return(rankPartsPtr(rankPartsRepGmp));
            } case PartitionType::Multiset: {
                return(rankPartsPtr(rankPartsMultisetGmp));
            } case PartitionType::CompRepNoZero: {
                return(rankPartsPtr(rankCompsRepGmp));
            } case PartitionType::CompRepCapped: {
                return(rankPartsPtr(rankCompsRepCappedGmp));
            } case PartitionType::CompRepWeak: {
                return(rankPartsPtr(rankCompsRepGmp));
            } case PartitionType::CompRepWeakCap: {
                return(rankPartsPtr(rankCompsRepCappedGmp));
            } case PartitionType::CompRepZero: {
                return(rankPartsPtr(rankCompsRepZeroGmp));
            } case PartitionType::CompRepCapZero: {
                return(rankPartsPtr(rankCompsRepZeroCapGmp));
            } case PartitionType::CompDistinctNoZero: {
                return(rankPartsPtr(rankCompsDistinctGmp));
            } case PartitionType::CompDistinctCapped: {
                return(rankPartsPtr(rankCompsDistinctGmp));
            } case PartitionType::CompDistinctWeak: {
                return(rankPartsPtr(rankCompsDistinctGmp));
            } case PartitionType::CompDistinctCapWeak: {
                return(rankPartsPtr(rankCompsDistinctGmp));
            } case PartitionType::CompDistinctZero: {
                return(rankPartsPtr(rankCompsDistinctMZGmp));
            } case PartitionType::CompDistinctCapMZ: {
                return(rankPartsPtr(rankCompsDistinctMZGmp));
            } case PartitionType::CompDistinctMZWeak: {
                return(rankPartsPtr(rankCompsDistinctMZWeakGmp));
            } case PartitionType::CompDistinctCapMZWeak: {
                return(rankPartsPtr(rankCompsDistinctMZWeakGmp));
            } default : {
                return nullptr;
            }
        }
    } else {
        switch (ptype) {
            case PartitionType::LengthOne: {
                return(rankPartsPtr(rankLengthOne));
            } case PartitionType::DistinctCapped: {
                return(rankPartsPtr(rankPartsDistinctCap));
            } case PartitionType::DistinctCappedMZ: {
                return(rankPartsPtr(rankPartsDistinctCapMZ));
            } case PartitionType::DistinctNoZero: {
                return(rankPartsPtr(rankPartsDistinctLen));
            } case PartitionType::DistinctOneZero: {
                return(rankPartsPtr(rankPartsDistinctOneZero));
            } case PartitionType::DistinctMZ: {
                return(rankPartsPtr(rankPartsDistinctMZ));
            } case PartitionType::DistinctStdAll: {
                return(rankPartsPtr(rankPartsDistinctMZ));
            } case PartitionType::RepCapped: {
                return(rankPartsPtr(rankPartsRepCap));
            } case PartitionType::RepNoZero: {
                return(rankPartsPtr(rankPartsRepLen));
            } case PartitionType::RepShort: {
                return(rankPartsPtr(rankPartsRepShort));
            } case PartitionType::RepStdAll: {
                return(rankPartsPtr(rankPartsRep));
            } case PartitionType::Multiset: {
                return(rankPartsPtr(rankPartsMultiset));
            } case PartitionType::CompRepNoZero: {
                return(rankPartsPtr(rankCompsRep));
            } case PartitionType::CompRepCapped: {
                return(rankPartsPtr(rankCompsRepCapped));
            } case PartitionType::CompRepWeak: {
                return(rankPartsPtr(rankCompsRep));
            } case PartitionType::CompRepWeakCap: {
                return(rankPartsPtr(rankCompsRepCapped));
            } case PartitionType::CompRepZero: {
                return(rankPartsPtr(rankCompsRepZero));
            } case PartitionType::CompRepCapZero: {
                return(rankPartsPtr(rankCompsRepZeroCap));
            } case PartitionType::CompDistinctNoZero: {
                return(rankPartsPtr(rankCompsDistinct));
            } case PartitionType::CompDistinctCapped: {
                return(rankPartsPtr(rankCompsDistinct));
            } case PartitionType::CompDistinctWeak: {
                return(rankPartsPtr(rankCompsDistinct));
            } case PartitionType::CompDistinctCapWeak: {
                return(rankPartsPtr(rankCompsDistinct));
            } case PartitionType::CompDistinctZero: {
                return(rankPartsPtr(rankCompsDistinctMZ));
            } case PartitionType::CompDistinctCapMZ: {
                return(rankPartsPtr(rankCompsDistinctMZ));
            } case PartitionType::CompDistinctMZWeak: {
                return(rankPartsPtr(rankCompsDistinctMZWeak));
            } case PartitionType::CompDistinctCapMZWeak: {
                return(rankPartsPtr(rankCompsDistinctMZWeak));
            } default : {
                return nullptr;
            }
        }
    }
}
