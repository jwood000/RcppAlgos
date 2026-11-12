#pragma once

// Cutoff definitions shared only by the temporary selector v3 comparison experiment.
#include "NumbersUtils/PartialSieveUtils.h"
#include "NumbersUtils/Phi.h"
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace P2SelectorV3Cases {
    inline std::int64_t ExtraTailTerms(const std::string& name) {
        if (name == "tail2") return 2;
        if (name == "tail4") return 4;
        if (name == "tail8") return 8;
        if (name == "tail16") return 16;
        if (name == "tail31") return 31;
        if (name == "tail64") return 64;
        if (name == "tail256") return 256;
        if (name == "tail1024") return 1024;
        if (name == "tail2048") return 2048;
        if (name == "tail4094") return 4094;
        if (name == "tail4097") return 4097;
        if (name == "tail8192") return 8192;
        return -1;
    }

    inline bool Valid(const std::string& name) {
        return name == "meissel" || name == "lehmer" || name == "zero" ||
            name == "one" || name == "half-meissel" || name == "double-meissel" ||
            name == "half-root" || name == "near-root" || name == "last" ||
            name == "empty" || name == "tail4095" || name == "tail4096" ||
            name == "quarter-root" || name == "three-quarter-root" || ExtraTailTerms(name) >= 0;
    }

    inline const char* Definition(const std::string& name) {
        if (name == "meissel") return "pi(floor(cuberoot(x)))";
        if (name == "lehmer") return "pi(floor(fourthroot(x)))";
        if (name == "zero") return "0";
        if (name == "one") return "1";
        if (name == "half-meissel") return "floor(meissel_a/2)";
        if (name == "double-meissel") return "min(b;2*meissel_a)";
        if (name == "half-root") return "floor(b/2)";
        if (name == "near-root") return "max(0;b-32)";
        if (name == "last") return "max(0;b-1)";
        if (name == "empty") return "b";
        if (name == "tail4095") return "max(0;b-4095)";
        if (name == "tail4096") return "max(0;b-4096)";
        if (name == "quarter-root") return "floor(b/4)";
        if (name == "three-quarter-root") return "floor(3*b/4)";
        if (name == "tail2") return "max(0;b-2)";
        if (name == "tail4") return "max(0;b-4)";
        if (name == "tail8") return "max(0;b-8)";
        if (name == "tail16") return "max(0;b-16)";
        if (name == "tail31") return "max(0;b-31)";
        if (name == "tail64") return "max(0;b-64)";
        if (name == "tail256") return "max(0;b-256)";
        if (name == "tail1024") return "max(0;b-1024)";
        if (name == "tail2048") return "max(0;b-2048)";
        if (name == "tail4094") return "max(0;b-4094)";
        if (name == "tail4097") return "max(0;b-4097)";
        if (name == "tail8192") return "max(0;b-8192)";
        throw std::runtime_error("unknown cutoff: " + name);
    }

    struct Values {
        std::int64_t root, b, meissel, lehmer;

        Values(std::int64_t sqrtX, std::int64_t piRoot,
               std::int64_t piCube, std::int64_t piFourth)
            : root(sqrtX), b(piRoot), meissel(piCube), lehmer(piFourth) {}

        std::int64_t get(const std::string& name) const {
            if (name == "meissel") return meissel;
            if (name == "lehmer") return lehmer;
            if (name == "zero") return 0;
            if (name == "one") return 1;
            if (name == "half-meissel") return meissel / 2;
            if (name == "double-meissel") return std::min(b, 2 * meissel);
            if (name == "half-root") return b / 2;
            if (name == "near-root") return std::max<std::int64_t>(0, b - 32);
            if (name == "last") return std::max<std::int64_t>(0, b - 1);
            if (name == "empty") return b;
            if (name == "tail4095") return std::max<std::int64_t>(0, b - 4095);
            if (name == "tail4096") return std::max<std::int64_t>(0, b - 4096);
            if (name == "quarter-root") return b / 4;
            if (name == "three-quarter-root") return 3 * b / 4;
            const auto terms = ExtraTailTerms(name);
            if (terms >= 0) return std::max<std::int64_t>(0, b - terms);
            throw std::runtime_error("unknown cutoff: " + name);
        }
    };

    inline Values Describe(std::int64_t x, const PrimeCounting::PhiTables& tables) {
        using namespace PrimeCounting::PartialSieveDetail;
        const auto root = IntegerSqrt(x);
        return Values(root, tables.pi()[root], tables.pi()[IntegerCbrt(x)],
                      tables.pi()[IntegerSqrt(root)]);
    }

    // Prepare only the requested cutoff before timing. Computing the other
    // definitions and diagnostic decisions waits until after the P2 call.
    inline std::int64_t Cutoff(const std::string& name, std::int64_t x,
                               const PrimeCounting::PhiTables& tables) {
        using namespace PrimeCounting::PartialSieveDetail;
        if (name == "zero") return 0;
        if (name == "one") return 1;
        if (name == "meissel" || name == "half-meissel") {
            const std::int64_t a = tables.pi()[IntegerCbrt(x)];
            return name == "meissel" ? a : a / 2;
        }
        const auto root = IntegerSqrt(x);
        if (name == "lehmer") return tables.pi()[IntegerSqrt(root)];
        const std::int64_t b = tables.pi()[root];
        if (name == "double-meissel") {
            return std::min(b, 2 * std::int64_t(tables.pi()[IntegerCbrt(x)]));
        }
        if (name == "half-root") return b / 2;
        if (name == "near-root") return std::max<std::int64_t>(0, b - 32);
        if (name == "last") return std::max<std::int64_t>(0, b - 1);
        if (name == "empty") return b;
        if (name == "tail4095") return std::max<std::int64_t>(0, b - 4095);
        if (name == "tail4096") return std::max<std::int64_t>(0, b - 4096);
        if (name == "quarter-root") return b / 4;
        if (name == "three-quarter-root") return 3 * b / 4;
        const auto terms = ExtraTailTerms(name);
        if (terms >= 0) return std::max<std::int64_t>(0, b - terms);
        throw std::runtime_error("unknown cutoff: " + name);
    }
}
