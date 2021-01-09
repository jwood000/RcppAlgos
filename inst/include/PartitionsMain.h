#ifndef PARTITIONS_MAIN_H
#define PARTITIONS_MAIN_H

#include "PartitionEnums.h"

namespace Partitions {
    template <typename typeRcpp, typename typeVector>
    typeRcpp PartitionsMain(const std::vector<typeVector> &v, std::vector<int> &z, const std::vector<int> &Reps,
                            PartitionType PartType, typeVector target, int lenV, int m, bool IsRep, bool IsMult,
                            double numParts, bool IsComb, bool xtraCol, bool bUserRows, bool mIsNull, bool getAll);
}

#endif
