#ifndef PARTITIONS_COUNT_H
#define PARTITIONS_COUNT_H

#include "Partitions/PartitionsTypes.h"

double PartitionsCount(const std::vector<int> &Reps, const PartDesign &part,
                       int lenV, bool bCalcMultiset);

#endif
