#ifndef PARTITION_DESIGN_H
#define PARTITION_DESIGN_H

#include "cpp11/R.hpp"

#include "Constraints/ConstraintsTypes.h"
#include "Partitions/PartitionsTypes.h"

SEXP GetDesign(const PartDesign &part, ConstraintType ctype,
               int lenV, bool verbose);

#endif
