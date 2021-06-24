#include "Partitions/PartitionsCountMultiset.h"
#include "Partitions/PartitionsCountDistinct.h"
#include "Partitions/PartitionsCountRep.h"
#include "Partitions/PartitionsTypes.h"
#include "CleanConvert.h"
#include <cmath>

double GetSpecialCount(const std::vector<int> &z, int target, int m) {

    double count = 0;
    const int startLen = std::count_if(z.cbegin(), z.cend(),
                                       [](int i){return i > 0;});

    for (int i = startLen; i <= m; ++i) {
        count += CountPartDistinctLen(target, i);
    }

    return count;
}

void GetSpecialCount(mpz_t res, const std::vector<int> &z,
                     int target, int m) {

    const int startLen = std::count_if(z.cbegin(), z.cend(),
                                       [](int i){return i > 0;});
    mpz_t temp;
    mpz_init(temp);

    for (int i = startLen; i <= m; ++i) {
        CountPartDistinctLen(temp, target, i);
        mpz_add(res, res, temp);
    }
}

void PartitionsCount(const std::vector<int> &Reps, PartDesign &part,
                     int lenV, bool bCalcMultiset, bool IsComb) {

    mpz_init(part.bigCount);

    if (IsComb) {
        switch (part.ptype) {
            case PartitionType::RepStdAll: {
                part.count = CountPartRep(part.mapTar);

                if (part.count > Significand53) {
                    part.isGmp = true;
                    CountPartRep(part.bigCount, part.mapTar);
                }

                break;
            } case PartitionType::RepNoZero: {
                part.count = CountPartRepLen(part.mapTar, part.width);

                if (part.count > Significand53) {
                    part.isGmp = true;
                    CountPartRepLen(part.bigCount, part.mapTar, part.width);
                }

                break;
            } case PartitionType::RepShort: {
                part.count = CountPartRepLen(part.mapTar, part.width);

                if (part.count > Significand53) {
                    part.isGmp = true;
                    CountPartRepLen(part.bigCount, part.mapTar, part.width);
                }

                break;
            } case PartitionType::RepCapped: {
                part.count = CountPartRepLenCap(part.mapTar,
                                                part.width, lenV);

                if (part.count > Significand53) {
                    part.isGmp = true;
                }

                break;
            } case PartitionType::DstctStdAll: {
                part.count = CountPartDistinct(part.mapTar);

                if (part.count > Significand53) {
                    part.isGmp = true;
                    CountPartDistinct(part.bigCount, part.mapTar);
                }

                break;
            } case PartitionType::DstctShort: {
                part.count = GetSpecialCount(part.startZ,
                                             part.mapTar, part.width);

                if (part.count > Significand53) {
                    part.isGmp = true;
                    GetSpecialCount(part.bigCount, part.startZ,
                                    part.mapTar, part.width);
                }

                break;
            } case PartitionType::DstctSpecial: {
                part.count = GetSpecialCount(part.startZ,
                                             part.mapTar, part.width);

                if (part.count > Significand53) {
                    part.isGmp = true;
                    GetSpecialCount(part.bigCount, part.startZ,
                                    part.mapTar, part.width);
                }

                break;
            } case PartitionType::DstctOneZero: {
                part.count = CountPartDistinctLen(part.mapTar, part.width);

                if (part.count > Significand53) {
                    part.isGmp = true;
                    CountPartDistinctLen(part.bigCount,
                                         part.mapTar, part.width);
                }

                break;
            } case PartitionType::DstctNoZero: {
                part.count = CountPartDistinctLen(part.mapTar, part.width);

                if (part.count > Significand53) {
                    part.isGmp = true;
                    CountPartDistinctLen(part.bigCount,
                                         part.mapTar, part.width);
                }

                break;
            } case PartitionType::DstctCapped: {
                part.count = CountPartDistinctLenCap(part.mapTar,
                                               part.width, lenV);

                if (part.count > Significand53) {
                    part.isGmp = true;
                }

                break;
            } case PartitionType::Multiset: {
                if (bCalcMultiset && part.solnExist) {
                    part.count = CountPartMultiset(Reps, part.startZ);
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
                part.count = CountPartPermRep(part.mapTar, part.width,
                                        part.includeZero);
            } else {
                part.count = 0.0;
            }
        } else if (it != DistPTypeArr.cend()) {
            part.count = CountPartPermDistinct(part.startZ, part.mapTar,
                                         part.width, part.includeZero);
        } else {
            part.count = 0.0;
        }
    }
}
