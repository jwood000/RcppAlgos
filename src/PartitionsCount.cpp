#include "Partitions/PartitionsCountMultiset.h"
#include "Partitions/PartitionsCountDistinct.h"
#include "Partitions/BigPartsCountDistinct.h"
#include "Partitions/PartitionsCountRep.h"
#include "Partitions/BigPartsCountRep.h"
#include "Partitions/PartitionsTypes.h"
#include "CleanConvert.h"  // Significand53
#include <algorithm>       // std::count_if, std::find

void PartitionsCount(const std::vector<int> &Reps, PartDesign &part,
                     int lenV, bool bCalcMultiset, bool IsComb) {

    mpz_init(part.bigCount);

    if (IsComb) {
        switch (part.ptype) {
            case PartitionType::RepStdAll: {
                part.count = CountPartsRep(part.mapTar);

                if (part.count > Significand53) {
                    part.isGmp = true;
                    CountPartsRep(part.bigCount, part.mapTar);
                }

                break;
            } case PartitionType::RepNoZero: {
                part.count = CountPartsRepLen(part.mapTar, part.width);

                if (part.count > Significand53) {
                    part.isGmp = true;
                    CountPartsRepLen(part.bigCount, part.mapTar, part.width);
                }

                break;
            } case PartitionType::RepShort: {
                part.count = CountPartsRepLen(part.mapTar, part.width);

                if (part.count > Significand53) {
                    part.isGmp = true;
                    CountPartsRepLen(part.bigCount, part.mapTar, part.width);
                }

                break;
            } case PartitionType::RepCapped: {
                part.count = CountPartsRepLenCap(part.mapTar,
                                                part.width, lenV);

                if (part.count > Significand53) {
                    part.isGmp = true;
                    CountPartsRepLenCap(part.bigCount, part.mapTar,
                                        part.width, lenV);
                }

                break;
            } case PartitionType::DstctStdAll: {
                part.count = CountPartsDistinct(part.mapTar);

                if (part.count > Significand53) {
                    part.isGmp = true;
                    CountPartsDistinct(part.bigCount, part.mapTar);
                }

                break;
            } case PartitionType::DstctSpecial: {
                const int strtLen = std::count_if(part.startZ.cbegin(),
                                                  part.startZ.cend(),
                                                  [](int i){return i > 0;});

                part.count = CountPartsDistinctMultiZero(
                    part.mapTar, part.width, strtLen
                );

                if (part.count > Significand53) {
                    part.isGmp = true;
                    CountPartsDistinctMultiZero(
                        part.bigCount, part.mapTar, part.width, strtLen
                    );
                }

                break;
            } case PartitionType::DstctOneZero: {
                part.count = CountPartsDistinctLen(part.mapTar, part.width);

                if (part.count > Significand53) {
                    part.isGmp = true;
                    CountPartsDistinctLen(part.bigCount,
                                         part.mapTar, part.width);
                }

                break;
            } case PartitionType::DstctNoZero: {
                part.count = CountPartsDistinctLen(part.mapTar, part.width);

                if (part.count > Significand53) {
                    part.isGmp = true;
                    CountPartsDistinctLen(part.bigCount,
                                          part.mapTar, part.width);
                }

                break;
            } case PartitionType::DstctCapped: {
                part.count = CountPartsDistinctLenCap(part.mapTar,
                                                      part.width, lenV);

                if (part.count > Significand53) {
                    part.isGmp = true;
                    CountPartsDistinctLenCap(part.bigCount, part.mapTar,
                                             part.width, lenV);
                }

                break;
            } case PartitionType::Multiset: {
                if (bCalcMultiset && part.solnExist) {
                    part.count = CountPartsMultiset(Reps, part.startZ);
                }

                break;
            } default: {
                part.count = 0.0;
                break;
            }
        }
    } else {
        const auto it = std::find(DistPTypeArr.cbegin(),
                                  DistPTypeArr.cend(), part.ptype);
        if (part.isRep) {
            if (part.ptype != PartitionType::RepCapped) {
                part.count = CountPartsPermRep(part.mapTar, part.width,
                                               part.includeZero);
            } else {
                part.count = 0.0;
            }
        } else if (it != DistPTypeArr.cend()) {
            part.count = CountPartsPermDistinct(part.startZ, part.mapTar,
                                                part.width,
                                                part.includeZero);
        } else {
            part.count = 0.0;
        }
    }
}
