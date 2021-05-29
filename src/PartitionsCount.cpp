#include "Partitions/PartitionsCountMultiset.h"
#include "Partitions/PartitionsCountDistinct.h"
#include "Partitions/PartitionsCountRep.h"
#include "Partitions/PartitionsTypes.h"
#include <cmath>

double GetSpecialCount(const std::vector<int> &z, int target, int m) {

    double count = 0;
    const int startLen = std::count_if(z.cbegin(), z.cend(),
                                       [](int i){return i > 0;});

    for (int i = startLen; i <= m; ++i)
        count += CountPartDistinctLen(target, i);

    return count;
}

double PartitionsCount(const std::vector<int> &Reps, const PartDesign &part,
                       int lenV, bool bCalcMultiset) {

    switch (part.ptype) {
        case PartitionType::RepStdAll: {
            return CountPartRep(part.mapTar);
        } case PartitionType::RepNoZero: {
            return CountPartRepLen(part.mapTar, part.width);
        } case PartitionType::RepShort: {
            return CountPartRepLen(part.mapTar, part.width);
        } case PartitionType::RepCapped: {
            return CountPartRepLenCap(part.mapTar, part.width, lenV);
        } case PartitionType::DstctStdAll: {
            return CountPartDistinct(part.mapTar);
        } case PartitionType::DstctShort: {
            return GetSpecialCount(part.startZ, part.mapTar, part.width);
        } case PartitionType::DstctSpecial: {
            return GetSpecialCount(part.startZ, part.mapTar, part.width);
        } case PartitionType::DstctOneZero: {
            return CountPartDistinctLen(part.mapTar, part.width);
        } case PartitionType::DstctNoZero: {
            return CountPartDistinctLen(part.mapTar, part.width);
        } case PartitionType::DistCapped: {
            return CountPartDistinctLenCap(part.mapTar, part.width, lenV);
        } case PartitionType::Multiset: {
            if (bCalcMultiset && part.solnExist) {
                return CountPartMultiset(Reps, part.startZ);
            }
        } default: {
            return 0.0;
        }
    }
}
