#include "NumbersUtils/PrimeCount.h"

namespace PrimeCounting {

    std::int64_t PiAuto(std::int64_t n, int nThreads,
                        const PhiTables& tables) {
        switch (SelectPrimeCountMethod(n, nThreads)) {
            case PrimeCountMethod::Meissel:
                return PiMeissel(n, nThreads, tables);
            case PrimeCountMethod::Lehmer:
                return PiLehmer(n, nThreads, tables);
            default:
                return PiLegendre(n, nThreads, tables);
        }
    }
}
