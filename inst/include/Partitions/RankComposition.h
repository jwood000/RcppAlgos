#pragma once

#include "Partitions/RankPartition.h"

// Composition implementations use the shared ranking signature.
rankPartsFunc rankCompsRep;
rankPartsFunc rankCompsRepGmp;
rankPartsFunc rankCompsRepCapped;
rankPartsFunc rankCompsRepCappedGmp;
rankPartsFunc rankCompsRepZero;
rankPartsFunc rankCompsRepZeroGmp;
rankPartsFunc rankCompsRepZeroCap;
rankPartsFunc rankCompsRepZeroCapGmp;
rankPartsFunc rankCompsDistinct;
rankPartsFunc rankCompsDistinctGmp;
rankPartsFunc rankCompsDistinctMZ;
rankPartsFunc rankCompsDistinctMZGmp;
rankPartsFunc rankCompsDistinctMZWeak;
rankPartsFunc rankCompsDistinctMZWeakGmp;
