#pragma once

#include "cpp11/R.hpp"

#include "Constraints/ConstraintsTypes.h"
#include "Partitions/PartitionsTypes.h"

SEXP GetDesign(const PartDesign &part, ConstraintType ctype,
               int lenV, bool verbose);
