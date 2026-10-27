#pragma once

#include "Partitions/NthPartition.h"

// Composition implementations use the shared unranking signature.
nthPartsFunc nthCompsRep;
nthPartsFunc nthCompsRepGmp;
nthPartsFunc nthCompsRepCapped;
nthPartsFunc nthCompsRepCappedGmp;
nthPartsFunc nthCompsRepZero;
nthPartsFunc nthCompsRepZeroGmp;
nthPartsFunc nthCompsRepZeroCap;
nthPartsFunc nthCompsRepZeroCapGmp;
nthPartsFunc nthCompsDistinct;
nthPartsFunc nthCompsDistinctGmp;
nthPartsFunc nthCompsDistinctMZ;
nthPartsFunc nthCompsDistinctMZGmp;
nthPartsFunc nthCompsDistinctMZWeak;
nthPartsFunc nthCompsDistinctMZWeakGmp;
