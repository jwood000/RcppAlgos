#ifndef PARTITION_UTILS_H
#define PARTITION_UTILS_H

#include "Constraints/ConstraintsTypes.h"
#include "Partitions/PartitionsCount.h"
#include "Cpp14MakeUnique.h"
#include "ImportExportMPZ.h"
#include <numeric>
#include <cmath>

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

void CheckPartition(const std::vector<std::string> &compFunVec,
                    const std::vector<double> &v, const std::string &mainFun,
                    const std::vector<double> &target, PartDesign &part,
                    SEXP Rlow, int lenV, int m, double tolerance,
                    bool IsBetween);

void SetPartitionDesign(const std::vector<int> &Reps,
                        const std::vector<double> &v, PartDesign &part,
                        ConstraintType &ctype, int lenV, int &m,
                        bool bCalcMultiset, bool IsComb);

#endif
