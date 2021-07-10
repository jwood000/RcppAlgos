#include "Partitions/PartitionsCountMultiset.h"
#include "Partitions/PartitionsCountDistinct.h"
#include "Partitions/BigPartsCountDistinct.h"
#include "Partitions/PartitionsCountRep.h"
#include "Partitions/BigPartsCountRep.h"
#include "Partitions/PartitionsTypes.h"
#include "Combinations/ComboCount.h"
#include "CleanConvert.h"  // Significand53
#include <algorithm>       // std::count_if, std::find

void PartitionsCount(const std::vector<int> &Reps, PartDesign &part,
                     int lenV, bool bCalcDifficult, bool IsComb) {

    mpz_init(part.bigCount);
    constexpr double cutOff = 3.0;
    const double capNumIters = static_cast<double>(part.mapTar + 1) *
                               static_cast<double>(part.width - 1) *
                               static_cast<double>(lenV + 1);

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
                const double theBar = NumCombsWithRep(lenV, part.width);

                if (bCalcDifficult || (theBar / capNumIters) > cutOff) {
                    part.count = CountPartsRepLenCap(part.mapTar,
                                                     part.width, lenV);

                    if (part.count > Significand53) {
                        part.isGmp = true;
                        CountPartsRepLenCap(part.bigCount, part.mapTar,
                                            part.width, lenV);
                    }
                } else {
                    part.numUnknown = true;
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
                const double theBar = nChooseK(lenV, part.width);

                if (bCalcDifficult || (theBar / capNumIters) > cutOff) {
                    part.count = CountPartsDistinctLenCap(part.mapTar,
                                                          part.width, lenV);

                    if (part.count > Significand53) {
                        part.isGmp = true;
                        CountPartsDistinctLenCap(part.bigCount, part.mapTar,
                                                 part.width, lenV);
                    }
                } else {
                    part.numUnknown = true;
                }

                break;
            } case PartitionType::Multiset: {
                if (bCalcDifficult && part.solnExist) {
                    part.count = CountPartsMultiset(Reps, part.startZ);
                } else {
                    part.numUnknown = true;
                }

                break;
            } default: {
                part.numUnknown = true;
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
                                               part.mapIncZero);
            } else {
                part.numUnknown = true;
                part.count = 0.0;
            }
        } else if (part.ptype == PartitionType::DstctCapped) {
            const double theBar = nChooseK(lenV, part.width);

            if (bCalcDifficult || (theBar / capNumIters) > cutOff) {
                part.count = CountPartsPermDistinctCap(part.startZ, lenV,
                                                       part.mapTar, part.width,
                                                       part.mapIncZero);
            } else {
                part.numUnknown = true;
                part.count = 0.0;
            }
        } else if (it != DistPTypeArr.cend()) {
            part.count = CountPartsPermDistinct(part.startZ, part.mapTar,
                                                part.width,
                                                part.mapIncZero);
        } else {
            part.numUnknown = true;
            part.count = 0.0;
        }
    }
}
