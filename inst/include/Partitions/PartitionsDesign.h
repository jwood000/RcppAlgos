#ifndef PARTITION_DESIGN_H
#define PARTITION_DESIGN_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

#include "Constraints/ConstraintsTypes.h"
#include "Partitions/PartitionsTypes.h"

SEXP GetDesign(const PartDesign &part, ConstraintType ctype,
               int lenV, bool verbose);

#endif
