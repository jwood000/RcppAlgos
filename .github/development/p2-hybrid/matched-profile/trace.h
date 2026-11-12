#pragma once
// Experiment-only diagnostics. Workers own disjoint slots, read after join.
#include <array>
#include <chrono>
#include <cstdint>
#include <vector>
#ifdef MATCHED_COUNTERS
#include "phi_metrics.h"
#endif

namespace Matched {
using Clock = std::chrono::steady_clock;
struct Chunk { std::int64_t first, last; double begin, end; };
struct Worker {
    double begin = 0, end = 0, seed = 0, setup = 0, fills = 0;
    std::uint64_t fillCount = 0;
    std::vector<Chunk> chunks;
#ifdef MATCHED_COUNTERS
    std::array<MatchedPhi::Stats, 3> stats; // prefix, Phi tail, stream seed
#endif
};
struct Trace {
    Clock::time_point origin;
    std::int64_t split = 0;
    int group = 0;
    std::array<std::array<Worker, 8>, 2> workers;
    std::array<double, 2> phaseBegin{{0,0}}, phaseEnd{{0,0}};
    explicit Trace(std::int64_t s) : split(s) {
#ifdef MATCHED_TRACE
        for (auto& group : workers) for (auto& worker : group)
            worker.chunks.reserve(1024);
#endif
    }
};
extern Trace* current;
extern thread_local Worker* worker;
inline double now() {
    return std::chrono::duration<double, std::micro>(Clock::now()-current->origin).count();
}
struct WorkerScope {
    Worker* previous;
    explicit WorkerScope(std::size_t slot) : previous(worker) {
        worker = &current->workers[current->group][slot];
        worker->begin = now();
    }
    ~WorkerScope() { worker->end = now(); worker = previous; }
};
struct ChunkScope {
    Chunk value;
    ChunkScope(std::int64_t first, std::int64_t last) : value{first,last,now(),0} {}
    ~ChunkScope() { value.end=now(); worker->chunks.push_back(value); }
};
struct PhaseScope {
    int group;
    explicit PhaseScope(int g) : group(g) { current->group=g; current->phaseBegin[g]=now(); }
    ~PhaseScope() { current->phaseEnd[group]=now(); }
};
struct FillScope {
    double begin;
    FillScope() : begin(now()) {}
    ~FillScope() { worker->fills += now()-begin; ++worker->fillCount; }
};
}

#ifdef MATCHED_TRACE
#define MP_WORKER(slot) Matched::WorkerScope mpWorker(slot)
#define MP_CHUNK(first,last) Matched::ChunkScope mpChunk(first,last)
#define MP_PHASE(group) Matched::PhaseScope mpPhase(group)
#define MP_FILL() Matched::FillScope mpFill
#else
#define MP_WORKER(slot)
#define MP_CHUNK(first,last)
#define MP_PHASE(group)
#define MP_FILL()
#endif
#ifdef MATCHED_COUNTERS
#define MP_QUERY(i) MatchedPhi::Scope mpQuery(Matched::worker->stats[(i) < Matched::current->split ? 0 : 1])
#else
#define MP_QUERY(i)
#endif
