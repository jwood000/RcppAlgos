#include "NumbersUtils/CacheInspector.h"

namespace CacheInspection {

    const std::array<const char*, PhaseCount> PhaseNames = {{
        "setup", "offsets", "buckets", "dense_mark", "bucket_mark",
        "extract", "reset"
    }};

    void PhaseTimer::Record(Clock::time_point end) {
        const std::size_t idx = static_cast<std::size_t>(phase);
        record->seconds[idx] += std::chrono::duration<double>(end - start).count();
        ++record->calls[idx];
    }

    KernelRecord* BeginKernel(std::vector<KernelRecord>* records,
                              const char* kernel, std::int_fast64_t lower,
                              std::int_fast64_t upper) {
        if (!records) return nullptr;

        records->emplace_back();
        KernelRecord& record = records->back();
        record.kernel = kernel;
        record.lower = lower;
        record.upper = upper;
        return &record;
    }
}
