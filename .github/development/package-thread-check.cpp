// Current-source launcher control-flow tests. Run the preparer first.
// Mathematical kernels below are stand-ins: package tests check their results.
// Exact launcher definitions and ThreadBatch are copied by the preparer.
// Gated real workers dereference borrowed matrices/vectors only after the caller
// has entered join(), exposing premature scratch destruction to ASan and TSan.
#include "RMatrix.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <exception>
#include <functional>
#include <iostream>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <string>
#include <system_error>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

namespace PackageThreadFaults {
enum class Fault { None, Worker, Allocation, NonStandard, Launch, Preparation };
struct WorkerError : std::runtime_error {
    WorkerError() : std::runtime_error("injected package worker failure") {}
};
struct PreparationError : std::runtime_error {
    PreparationError() : std::runtime_error("injected next-state preparation failure") {}
};
struct LaunchError : std::system_error {
    LaunchError() : std::system_error(
        std::make_error_code(std::errc::resource_unavailable_try_again),
        "injected package thread launch failure") {}
};
struct State {
    int attempts = 0, started = 0, joined = 0, joinAttempts = 0;
    int exited = 0, bodies = 0, preparations = 0, finalTouches = 0;
    bool returned = false, released = false, timedOut = false;
};
std::mutex mutex;
std::condition_variable changed;
State state;
Fault fault = Fault::None;
int failureSlot = 0;
bool gated = false;
thread_local int workerSlot = -1;
const std::chrono::seconds TIMEOUT(15);

void Check(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}
void Configure(Fault selected, int slot, bool gate) {
    std::lock_guard<std::mutex> lock(mutex);
    state = State(); fault = selected; failureSlot = slot; gated = gate;
}
State Inspect() {
    std::lock_guard<std::mutex> lock(mutex); return state;
}
void Preparation() {
    // All preparation calls occur on the calling thread.
    std::lock_guard<std::mutex> lock(mutex);
    const int slot = state.preparations++;
    if (fault == Fault::Preparation && slot == failureSlot) throw PreparationError();
}
void Body() {
    std::lock_guard<std::mutex> lock(mutex);
    ++state.bodies;
    if (workerSlot == failureSlot) {
        if (fault == Fault::Worker) throw WorkerError();
        if (fault == Fault::Allocation) throw std::bad_alloc();
        if (fault == Fault::NonStandard) throw 7162;
    }
}
void FinalTouch() {
    std::lock_guard<std::mutex> lock(mutex);
    Check(state.started == state.joined, "result processing before all joins");
    ++state.finalTouches;
}
class Thread {
    std::thread actual;
public:
    template <typename Function, typename... Args>
    explicit Thread(Function&& work, Args&&... args) {
        int slot;
        {
            std::lock_guard<std::mutex> lock(mutex);
            slot = state.attempts++;
            if (fault == Fault::Launch && slot == failureSlot) throw LaunchError();
        }
        actual = std::thread([slot](typename std::decay<Function>::type&& call,
                                   typename std::decay<Args>::type&&... values) {
            workerSlot = slot;
            {
                std::unique_lock<std::mutex> lock(mutex);
                ++state.started; changed.notify_all();
                if (gated && !changed.wait_for(lock, TIMEOUT, [] { return state.released; })) {
                    state.timedOut = true;
                }
            }
            call(std::move(values)...); // The real ThreadBatch catches worker exceptions.
            {
                std::lock_guard<std::mutex> lock(mutex);
                ++state.exited; changed.notify_all();
            }
        }, std::forward<Function>(work), std::forward<Args>(args)...);
    }
    Thread(Thread&&) = default;
    Thread& operator=(Thread&&) = default;
    Thread(const Thread&) = delete;
    Thread& operator=(const Thread&) = delete;
    bool joinable() const { return actual.joinable(); }
    void join() {
        {
            std::lock_guard<std::mutex> lock(mutex);
            ++state.joinAttempts; changed.notify_all();
        }
        actual.join();
        {
            std::lock_guard<std::mutex> lock(mutex);
            ++state.joined; changed.notify_all();
        }
    }
};
void Returned() {
    std::lock_guard<std::mutex> lock(mutex);
    state.returned = true; changed.notify_all();
}
bool ObserveAndRelease(int expected) {
    std::unique_lock<std::mutex> lock(mutex);
    const bool ready = changed.wait_for(lock, TIMEOUT, [&] {
        return state.returned || (state.started == expected && state.joinAttempts > 0);
    });
    const bool correct = ready && (expected == 0 ? state.returned :
        !state.returned && state.exited == 0 && state.joined == 0);
    state.released = true; changed.notify_all();
    return correct;
}
void CheckJoined() {
    const auto s = Inspect();
    Check(!s.timedOut && s.started == s.bodies && s.started == s.exited &&
          s.started == s.joined, "workers not joined before propagation");
}
}

// Lightweight type/preparation stand-ins keep the launcher tests independent of
// the R runtime and GMP. IsGmp branches are exercised as control flow; exact GMP
// mathematics is covered by the package regression tests, not this fixture.
using mpz_class = long long;
#ifndef PACKAGE_THREAD_R_BRIDGE
using SEXP = void*;
#endif
using nthResultPtr = int;
using rankResultPtr = int;
using rankPartsPtr = int;
using nextIterPtr = int;
using nextGrpFunc = int;
using nthFuncDbl = int;
using nthFuncGmp = int;
using finalTouchFunc = void(*)(SEXP, bool, int, bool,
    const std::vector<double>&, const std::vector<mpz_class>&, bool);
template <typename T> using funcPtr = int;
template <typename T> using compPtr = int;
enum class PartitionType { Fixture };
struct PartDesign {
    bool isComp = false, isGmp = false;
    int width = 2, mapTar = 20, maxZeros = 0;
    PartitionType ptype = PartitionType::Fixture;
};
using nthPartsPtr = std::vector<int>(*)(int, int, int, int, double,
                                      const mpz_class&, const std::vector<int>&);

namespace LauncherStandins {
using PackageThreadFaults::Check;
template <typename T> void Touch(const std::vector<T>& values) {
    // Force an actual read after the join gate opens.
    volatile std::size_t count = values.size();
    if (count) { volatile T value = values.front(); (void)value; }
}
template <typename T> void Matrix(RcppParallel::RMatrix<T>& matrix,
                                  int begin, int end, const std::vector<int>* z = nullptr) {
    Check(begin >= 0 && end <= static_cast<int>(matrix.nrow()) && begin <= end,
          "invalid chunk bounds");
    if (z) { Touch(*z); Check(!z->empty() && z->front() == begin, "incorrect starting-state capture"); }
    for (int row = begin; row < end; ++row) matrix(row, 0) = row + 1;
    PackageThreadFaults::Body();
}
template <typename T> void Rank(std::vector<T>& result, std::vector<int>& indices,
                                const std::vector<int>& reps, int rows) {
    Touch(reps); Touch(indices);
    Check(rows == static_cast<int>(result.size()), "incorrect rank section size");
    for (int row = 0; row < rows; ++row) result[row] = indices[row * 2] + 1;
    PackageThreadFaults::Body();
}
}

int GetNthResultFunc(bool, bool, bool, bool) { return 0; }
template <typename T> int GetFuncPtr(const std::string&) { return 0; }
template <typename T> int GetCompPtr(const std::string&) { return 0; }
int GetNextIterPtr(bool, bool, bool, bool) { return 0; }
void SetNextIter(const std::vector<int>&, std::vector<int>& z, int,
                 double& lower, mpz_class& lowerMpz, int step,
                 int, int, bool isGmp, bool, bool, bool) {
    PackageThreadFaults::Preparation();
    if (isGmp) lowerMpz += step; else lower += step;
    z[0] = isGmp ? lowerMpz : lower;
}
std::vector<int> NthParts(int, int width, int, int, double lower,
                        const mpz_class& lowerMpz, const std::vector<int>&) {
    PackageThreadFaults::Preparation();
    std::vector<int> z(width); z[0] = std::max<mpz_class>(lower, lowerMpz); return z;
}
nthPartsPtr GetNthPartsFuncOrStop(PartitionType, bool) { return NthParts; }
void GetStartProd(const std::vector<int>&, std::vector<int>& z,
                  mpz_class& lowerMpz, double& lower, int step, bool isGmp) {
    PackageThreadFaults::Preparation();
    if (isGmp) lowerMpz += step; else lower += step;
    z[0] = isGmp ? lowerMpz : lower;
}
void GetStartGrp(int, int, std::vector<int>& z, mpz_class& lowerMpz,
                 double& lower, int step, bool isGmp) {
    GetStartProd({}, z, lowerMpz, lower, step, isGmp);
}
void Create2D(const std::vector<int>& idx, std::vector<std::vector<int>>& sections,
              int step, int width, int threads) {
    for (int slot = 0; slot < threads; ++slot) {
        const int begin = slot * step * width;
        const int end = slot + 1 == threads ? idx.size() : begin + step * width;
        sections[slot].assign(idx.begin() + begin, idx.begin() + end);
    }
}
void SetSampleNames(SEXP, bool, int, const std::vector<double>&,
                    const std::vector<mpz_class>&, bool) { PackageThreadFaults::FinalTouch(); }
void GroupFinalTouch(SEXP, bool, int, bool, const std::vector<double>&,
                      const std::vector<mpz_class>&, bool) { PackageThreadFaults::FinalTouch(); }

// Serial branches must compile but are deliberately outside this launch test.
template <typename... Args> void ComboManager(Args&&...) {}
template <typename... Args> void PermuteManager(Args&&...) {}
template <typename... Args> void ComboResStd(Args&&...) {}
template <typename... Args> void PermuteResStd(Args&&...) {}
template <typename... Args> void PartsStdManager(Args&&...) {}
template <typename... Args> void PartsGenManager(Args&&...) {}
template <typename... Args> void SampleResults(Args&&...) {}
template <typename... Args> void SerialGlue(Args&&...) {}
template <typename... Args> void CnstrntSpcWorker(Args&&...) {}

template <typename T>
void ComboParallel(RcppParallel::RMatrix<T>& mat, const std::vector<T>& v,
    std::vector<int>& z, int, int, int begin, int end,
    const std::vector<int>& freqs, bool, bool) {
    LauncherStandins::Touch(v); LauncherStandins::Touch(freqs);
    LauncherStandins::Matrix(mat, begin, end, &z);
}
template <typename T, typename... Args> void PermuteParallel(RcppParallel::RMatrix<T>& mat, Args&&... args) {
    ComboParallel<T>(mat, std::forward<Args>(args)...);
}
template <typename T>
void ComboResPar(RcppParallel::RMatrix<T>& mat, const std::vector<T>& v,
    std::vector<int>& z, int n, int m, int begin, int end,
    const std::vector<int>& freqs, int, bool mult, bool rep) {
    ComboParallel(mat, v, z, n, m, begin, end, freqs, mult, rep);
}
template <typename T, typename... Args> void PermuteResPar(RcppParallel::RMatrix<T>& mat, Args&&... args) {
    ComboResPar<T>(mat, std::forward<Args>(args)...);
}
void PartsStdParallel(RcppParallel::RMatrix<int>& mat, std::vector<int>& z,
    int begin, int, int, int, int end, PartitionType, int) {
    LauncherStandins::Matrix(mat, begin, end, &z);
}
template <typename T>
void PartsGenParallel(RcppParallel::RMatrix<T>& mat, const std::vector<T>& v,
    const std::vector<int>& reps, std::vector<int>& z,
    int begin, int, int, int, int end, PartitionType, int) {
    LauncherStandins::Touch(v); LauncherStandins::Touch(reps);
    LauncherStandins::Matrix(mat, begin, end, &z);
}
template <typename T>
void CnstrntLowerWorker(const std::vector<T>& v, const std::vector<T>& targets,
    const std::vector<int>& freqs, const std::vector<std::string>& comps,
    std::vector<T>& constraint, std::vector<T>& result, std::vector<int>& z,
    int, int, int, int, int, int, int rows, bool) {
    LauncherStandins::Touch(v); LauncherStandins::Touch(targets);
    LauncherStandins::Touch(freqs);
    PackageThreadFaults::Check(!comps.empty(), "lost comparison storage");
    const int begin = z.at(0);
    for (int i = 0; i < rows; ++i) { result.push_back(begin + i + 1); constraint.push_back(i); }
    PackageThreadFaults::Body();
}
template <typename T>
void ParallelGlue(RcppParallel::RMatrix<T>& mat, const std::vector<T>& v,
    const std::vector<double>& sample, const std::vector<mpz_class>& big,
    const std::vector<int>& reps, int, int, int begin, int end, int, bool) {
    LauncherStandins::Touch(v); LauncherStandins::Touch(sample);
    LauncherStandins::Touch(big); LauncherStandins::Touch(reps);
    LauncherStandins::Matrix(mat, begin, end);
}
template <typename T>
void ParallelGlue(RcppParallel::RMatrix<T>& mat, const std::vector<T>& v,
    const std::vector<double>& sample, const std::vector<mpz_class>& big,
    const std::vector<int>& reps, nthPartsPtr, int, int begin, int end, int, int, int, bool) {
    LauncherStandins::Touch(v); LauncherStandins::Touch(sample);
    LauncherStandins::Touch(big); LauncherStandins::Touch(reps);
    LauncherStandins::Matrix(mat, begin, end);
}
template <typename T>
void ParallelGlue(RcppParallel::RMatrix<T>& mat, const std::vector<int>& idx,
    const std::vector<int>& lengths, const std::vector<T>& v,
    const std::vector<int>& next, const std::vector<double>& sample,
    const std::vector<mpz_class>& big, std::vector<int> z,
    int, int begin, int end, bool, bool) {
    LauncherStandins::Touch(idx); LauncherStandins::Touch(lengths); LauncherStandins::Touch(v);
    LauncherStandins::Touch(next); LauncherStandins::Touch(sample); LauncherStandins::Touch(big);
    LauncherStandins::Matrix(mat, begin, end, &z);
}
template <typename T>
void ParallelGlue(RcppParallel::RMatrix<T>& mat, const std::vector<T>& v,
    int, int, int, const std::vector<double>& sample,
    const std::vector<mpz_class>& big, std::vector<int> z,
    int, int begin, int end, bool, bool) {
    LauncherStandins::Touch(v); LauncherStandins::Touch(sample); LauncherStandins::Touch(big);
    LauncherStandins::Matrix(mat, begin, end, &z);
}
template <typename T>
void RankCPResultsGeneric(std::vector<T>& result, std::vector<int>& idx,
    const std::vector<int>& reps, int, int, int, int rows) {
    LauncherStandins::Rank(result, idx, reps, rows);
}
template <typename T>
void RankPartsResultsGeneric(std::vector<T>& result, std::vector<int>& idx,
    const std::vector<int>& reps, int, int, int, int, int, int rows) {
    LauncherStandins::Rank(result, idx, reps, rows);
}
namespace PrimeSieve {
template <typename T> void sqrtBigPrimes(int, bool, bool, bool, std::vector<T>& primes) {
    primes = {2, 3, 5};
}
}
namespace MotleyPrimes {
template <typename T, typename U>
void EulerPhiSieve(T lower, U upper, T offset, std::vector<T>& primes,
                   std::vector<T>& nums, U* output) {
    LauncherStandins::Touch(primes); LauncherStandins::Touch(nums);
    for (T i = 0; i <= upper - lower; ++i) output[offset + i] = lower + i;
    PackageThreadFaults::Body();
}
template <typename T>
void PrimeFactorizationSieve(T lower, T upper, T offset,
    const std::vector<T>& primes, std::vector<std::vector<T>>& output) {
    LauncherStandins::Touch(primes);
    for (T i = 0; i <= upper - lower; ++i) output[offset + i] = {lower + i};
    PackageThreadFaults::Body();
}
}
template <typename T, typename U>
void DivisorsSieve(T lower, U upper, T offset, std::vector<std::vector<U>>& output) {
    for (T i = 0; i <= upper - lower; ++i) output[offset + i] = {lower + i};
    PackageThreadFaults::Body();
}
template <typename T, typename U>
void NumDivisorsSieve(T lower, T upper, T offset, U* output) {
    for (T i = 0; i <= upper - lower; ++i) output[offset + i] = lower + i;
    PackageThreadFaults::Body();
}
template <typename T>
void PrimeFacList(std::size_t begin, std::size_t end,
    const std::vector<double>& nums, std::vector<std::vector<T>>& output) {
    for (std::size_t i = begin; i < end; ++i) output[i] = {static_cast<T>(nums[i])};
    PackageThreadFaults::Body();
}
template <typename T>
void FactorList(std::size_t begin, std::size_t end,
    const std::vector<double>& nums, std::vector<std::vector<T>>& output) {
    PrimeFacList<T>(begin, end, nums, output);
}
void IsPrimeVec(std::size_t begin, std::size_t end, const std::vector<double>& nums, int* output) {
    for (std::size_t i = begin; i < end; ++i) output[i] = nums[i];
    PackageThreadFaults::Body();
}

#include "package-thread-core.inc"

namespace PackageThreadChecks {
using PackageThreadFaults::Check;
const int ROWS = 11;
int THREADS = 3;
bool HasPreparation(int fixture) { return (fixture >= 3 && fixture <= 9) || fixture >= 14; }
int Variants(int fixture) { return fixture == 2 ? 3 : 2; }
void Invoke(int fixture, int variant) {
    std::vector<int> matrix(ROWS * 2, 0), v{1, 2}, reps{1, 1}, z(2, 0);
    std::vector<double> sample{0, 1}; std::vector<mpz_class> big{0, 1};
    mpz_class lowerMpz = 0;
    const bool isGmp = variant != 0;
    if (fixture <= 2) {
        const int count = fixture == 2 ? ROWS : 30005;
        std::vector<int> out(count, 0), nums(count, 1);
        std::vector<std::vector<int>> lists(count);
        if (fixture == 0) Launcher0::MotleyMain(1, count, variant == 0,
            out.data(), nums, lists, THREADS, THREADS);
        if (fixture == 1) Launcher1::DivisorMain(1, count, variant == 0,
            out.data(), lists, count, THREADS, THREADS);
        if (fixture == 2) {
            std::vector<double> inputs(count); std::iota(inputs.begin(), inputs.end(), 1);
            Launcher2::PollardRhoMain(inputs, count, variant == 0, variant == 1,
                lists, out.data(), count, THREADS, THREADS);
        }
        const bool scalar = fixture == 0 ? variant == 0 : fixture == 1 ? variant != 0 : variant == 2;
        for (int i = 0; i < count; ++i) {
            Check(scalar ? out[i] == i + 1 : lists[i] == std::vector<int>{i + 1},
                  "number-theory chunk coverage/capture changed");
        }
        return;
    }
    switch (fixture) {
    case 3: Launcher3::ThreadSafeCombinations(matrix.data(), v, 2, 2, true,
        false, false, isGmp, reps, z, reps, 0, lowerMpz, ROWS, THREADS); break;
    case 4: Launcher4::ThreadSafePermutations(matrix.data(), v, 2, 2, 0, false,
        true, false, false, isGmp, reps, z, reps, 0, lowerMpz, ROWS, THREADS); break;
    case 5: Launcher5::ComboResMain(matrix.data(), v, 0, 2, 2, true, false, false,
        isGmp, reps, z, reps, 0, lowerMpz, ROWS, THREADS); break;
    case 6: Launcher6::PermuteResMain(matrix.data(), v, 0, 2, 2, true, false, false,
        isGmp, reps, z, reps, 0, lowerMpz, ROWS, THREADS); break;
    case 7: {
        std::vector<int> constraint, output;
        Launcher7::ConstraintsSpecial(v, v, {"=="}, reps, reps, constraint, output,
            "sum", z, 0, lowerMpz, 2, 2, ROWS, THREADS, false, false, true,
            false, isGmp, true);
        Check(output.size() == ROWS, "special result size changed");
        for (int i = 0; i < ROWS; ++i) Check(output[i] == i + 1, "special starting-state capture changed");
        return;
    }
    case 8: Launcher8::StandardPartitions(matrix.data(), z, PartitionType::Fixture,
        0, lowerMpz, 2, 2, ROWS, THREADS, 1, 1, 20, 2, 2, isGmp,
        true, true, false, 0); break;
    case 9: {
        PartDesign part; part.isGmp = isGmp;
        Launcher9::GeneralPartitions(matrix.data(), v, reps, z, part, 0,
            lowerMpz, 2, ROWS, THREADS, 1, 1, 2, 2, true); break;
    }
    case 10: Launcher10::ThreadSafeSample(matrix.data(), nullptr, v, sample, big,
        reps, 0, 2, ROWS, THREADS, true, true, isGmp, 2); break;
    case 11: Launcher11::ThreadSafeSample(matrix.data(), nullptr, v, sample, big,
        reps, NthParts, 2, ROWS, THREADS, true, true, 20, 2, 2, isGmp); break;
    case 12: case 13: {
        const int step = ROWS / THREADS, prefix = (THREADS - 1) * step;
        const int last = ROWS - prefix;
        std::vector<std::vector<int>> sections(THREADS - 1, std::vector<int>(step));
        std::vector<int> tail(last), idx(ROWS * 2);
        for (int i = 0; i < ROWS; ++i) idx[2 * i] = i;
        if (fixture == 12) Launcher12::RankCPParallel(sections, tail, idx, reps,
            0, 2, 2, ROWS, THREADS, step, last);
        else Launcher13::RankPartsParallel(sections, tail, idx, reps, 0, 20,
            2, 2, 2, ROWS, THREADS, step, last);
        for (int i = 0; i < ROWS; ++i) Check((i < prefix ? sections[i / step][i % step] :
            tail[i - prefix]) == i + 1, "rank chunk capture changed");
        return;
    }
    case 14: Launcher14::PureOutputMain(matrix.data(), reps, reps, v, reps, sample,
        big, z, 2, ROWS, THREADS, true, lowerMpz, 0, false, isGmp); break;
    case 15: Launcher15::GroupsMain(matrix.data(), nullptr, 0, 0, 0, GroupFinalTouch,
        v, z, sample, big, lowerMpz, 0, 2, ROWS, THREADS, false, true, true,
        false, isGmp); break;
    default: throw std::logic_error("unknown launcher fixture");
    }
    for (int i = 0; i < ROWS; ++i) Check(matrix[i] == i + 1, "matrix chunk coverage changed");
}

void Run(int fixture, int variant, PackageThreadFaults::Fault fault, int slot) {
    using namespace PackageThreadFaults;
    Configure(fault, slot, true);
    std::exception_ptr caught;
    std::thread caller([&] {
        try { Invoke(fixture, variant); } catch (...) { caught = std::current_exception(); }
        Returned();
    });
    const int expected = fault == Fault::Launch ? slot :
        fault == Fault::Preparation ? slot + 1 : THREADS;
    const bool blocked = ObserveAndRelease(expected);
    caller.join();
    Check(blocked, "caller escaped before gated workers joined");
    CheckJoined();
    const auto s = Inspect();
    Check(s.started == expected && s.attempts == expected + (fault == Fault::Launch),
          "incorrect number of launches around failure");
    if (fault == Fault::None) { if (caught) std::rethrow_exception(caught); }
    else {
        Check(static_cast<bool>(caught), "injected failure did not reach caller");
        bool correct = false;
        try { std::rethrow_exception(caught); }
        catch (const WorkerError&) { correct = fault == Fault::Worker; }
        catch (const LaunchError&) { correct = fault == Fault::Launch; }
        catch (const PreparationError&) { correct = fault == Fault::Preparation; }
        catch (const std::bad_alloc&) { correct = fault == Fault::Allocation; }
        catch (int value) { correct = fault == Fault::NonStandard && value == 7162; }
        Check(correct, "wrong exception propagated");
        Check(s.finalTouches == 0, "partial results processed after failure");
    }
}
}

#ifndef PACKAGE_THREAD_NO_MAIN
int main(int argc, char** argv) {
    try {
        using namespace PackageThreadFaults;
        using namespace PackageThreadChecks;
        if (argc == 3 && std::string(argv[1]) == "--threads" &&
            (std::string(argv[2]) == "2" || std::string(argv[2]) == "3")) {
            THREADS = std::stoi(argv[2]);
        } else if (argc != 1) {
            throw std::runtime_error("usage: checker [--threads 2|3]");
        }
        int cases = 0, recoveries = 0, preparations = 0;
        for (int fixture = 0; fixture < 16; ++fixture) {
            for (int variant = 0; variant < Variants(fixture); ++variant) {
                Run(fixture, variant, Fault::None, 0); ++recoveries;
                for (auto fault : {Fault::Worker, Fault::Allocation, Fault::NonStandard, Fault::Launch}) {
                    for (int slot : {0, THREADS - 1}) {
                        Run(fixture, variant, fault, slot); ++cases;
                        Run(fixture, variant, Fault::None, 0); ++recoveries;
                    }
                }
                if (HasPreparation(fixture)) for (int slot = 0; slot < THREADS - 1; ++slot) {
                    Run(fixture, variant, Fault::Preparation, slot); ++cases; ++preparations;
                    Run(fixture, variant, Fault::None, 0); ++recoveries;
                }
            }
        }
        std::cout << "PASS launchers=16 variants=33 threads=" << THREADS << " failure_scenarios=" << cases
                  << " preparation_failures=" << preparations << " recoveries=" << recoveries << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n'; return 1;
    }
}
#endif
