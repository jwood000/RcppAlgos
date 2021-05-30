#ifndef PARTITION_UTILS_H
#define PARTITION_UTILS_H

#include "Partitions/PartitionsTypes.h"

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

bool CheckPartition(const std::vector<std::string> &compFunVec,
                    const std::vector<double> &v, const std::string &mainFun,
                    const std::vector<double> &target, PartDesign &part,
                    ConstraintType &ctype, SEXP Rlow, int lenV, int m,
                    double tolerance, bool IsBetween);

void GetPartitionDesign(const std::vector<int> &Reps,
                        const std::vector<double> &v, PartDesign &part,
                        ConstraintType &ctype, int lenV, int &m,
                        bool bCalcMultiset);

#endif
