#pragma once
// Temporary diagnostics for the unmodified P2 scheduling algorithm.
// Each query/chunk is claimed once; workers publish only disjoint record slots.
#include <array>
#include <chrono>
#include <cstdint>
#include <cstddef>

namespace QueryProfile {
using Clock = std::chrono::steady_clock;
struct Query {
    std::int64_t index=0, prime=0, quotient=0, cutoff=0, pi=0, term=0;
    int worker=-1;
    double begin=0, end=0;
};
struct Chunk {
    std::int64_t first=0, last=0;
    int worker=-1;
    double begin=0, end=0;
};
struct Worker {
    double begin=0, end=0;
    std::int64_t chunks=0, queries=0, finalFirst=0, finalLast=0;
};
struct Trace {
    Clock::time_point origin;
    std::int64_t first=0, last=0;
    std::size_t workerCount=0;
    std::array<Query,32> queries;
    std::array<Chunk,64> chunks;
    std::array<Worker,32> workers;
};
extern Trace* current;
extern thread_local int workerSlot;
inline double now() {
    return std::chrono::duration<double,std::micro>(Clock::now()-current->origin).count();
}
inline void range(std::int64_t first,std::int64_t last) {
    current->first=first; current->last=last;
}
inline bool measured(std::int64_t index) {
    return index-current->first<32;
}
inline double queryBegin(std::int64_t index) {
    return measured(index)?now():0;
}
inline void queryEnd(std::int64_t index,std::int64_t prime,std::int64_t quotient,
                     std::int64_t cutoff,std::int64_t pi,double begin) {
    if(!measured(index)) return;
    const double end=now();
    auto& q=current->queries[static_cast<std::size_t>(index-current->first)];
    q.index=index; q.prime=prime; q.quotient=quotient; q.cutoff=cutoff;
    q.pi=pi; q.term=pi-(index-1); q.worker=workerSlot; q.begin=begin; q.end=end;
}
struct WorkerScope {
    int previous;
    explicit WorkerScope(std::size_t slot):previous(workerSlot) {
        workerSlot=static_cast<int>(slot);
        current->workers[slot].begin=now();
    }
    ~WorkerScope() { current->workers[workerSlot].end=now(); workerSlot=previous; }
};
struct ChunkScope {
    Chunk* record;
    ChunkScope(std::int64_t first,std::int64_t last):record(nullptr) {
        auto& worker=current->workers[workerSlot];
        ++worker.chunks; worker.queries+=last-first;
        worker.finalFirst=first; worker.finalLast=last;
        const auto index=static_cast<std::size_t>((first-current->first)/32);
        if(index>=current->chunks.size()) return;
        record=&current->chunks[index];
        record->first=first; record->last=last; record->worker=workerSlot;
        record->begin=now();
    }
    ~ChunkScope() { if(record) record->end=now(); }
};
}
