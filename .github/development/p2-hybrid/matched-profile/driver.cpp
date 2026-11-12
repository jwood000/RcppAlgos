// Native profiling harness for one frozen matched pair. Not a package API.
#include "P2Instrumented.cpp"
#include <iomanip>
#include <iostream>
#include <string>
#include <stdexcept>

namespace Matched {
Trace* current = nullptr;
thread_local Worker* worker = nullptr;
}
namespace PrimeCounting {
std::int64_t P2Production(std::int64_t, std::int64_t, int, const PhiTables&);
}
namespace {
constexpr std::int64_t X = 80000000000LL;
constexpr std::int64_t LIMIT = 339410;
using PrimeCounting::PhiTables;
struct Geometry {
    std::int64_t a, first, split, last, b, direct;
    Geometry(const std::string& cutoff, const PhiTables& tables) {
        if (cutoff != "zero" && cutoff != "meissel") throw std::runtime_error("bad cutoff");
        const auto& primes = tables.primes();
        const auto& pi = tables.pi();
        a = cutoff == "zero" ? 0 : pi[PrimeCounting::PartialSieveDetail::IntegerCbrt(X)];
        first = a+1;
        b = pi[PrimeCounting::PartialSieveDetail::IntegerSqrt(X)];
        last = std::upper_bound(primes.begin()+first,primes.begin()+b+1,X/(LIMIT+1))-primes.begin();
        split = std::upper_bound(primes.begin()+first,primes.begin()+last,X/(16*LIMIT+1))-primes.begin();
        direct = 0;
        for (auto i=last;i<=b;++i) direct += static_cast<std::int64_t>(pi[X/primes[i]])-(i-1);
        if (split!=1725 || last!=20886 || b!=24668 || (a!=0 && a!=590))
            throw std::runtime_error("matched geometry changed");
    }
};
std::int64_t run(const std::string& mode, const Geometry& g, int threads, const PhiTables& t) {
    Matched::current->group=0;
    const auto budget=static_cast<std::size_t>(threads);
    if (mode=="baseline") return PrimeCounting::P2Production(X,g.a,threads,t);
    if (mode=="hybrid") return PrimeCounting::P2(X,g.a,threads,t);
    if (mode=="prefix") return PrimeCounting::countPhiPrefix(X,g.first,g.split,budget,t);
    if (mode=="tail-phi") return PrimeCounting::countPhiPrefix(X,g.split,g.last,budget,t);
    if (mode=="tail-stream") return PrimeCounting::countStreamingTail(X,g.split,g.last,threads,t);
    if (mode=="split-phi") {
        // Match the direct suffix work inside baseline/hybrid's timed call.
        std::int64_t direct=0;
        for(auto i=g.last;i<=g.b;++i)
            direct += static_cast<std::int64_t>(t.pi()[X/t.primes()[i]])-(i-1);
        auto result=PrimeCounting::countPhiPrefix(X,g.first,g.split,budget,t);
        Matched::current->group=1;
        return direct+result+PrimeCounting::countPhiPrefix(X,g.split,g.last,budget,t);
    }
    throw std::runtime_error("bad mode");
}
void writeTrace(const Matched::Trace& trace) {
    std::cout << ",\"phases\":[";
    for (int g=0;g<2;++g) {
        if(g) std::cout << ',';
        std::cout << '[' << trace.phaseBegin[g] << ',' << trace.phaseEnd[g] << ']';
    }
    std::cout << "],\"workers\":[";
    bool comma=false;
    for (int g=0;g<2;++g) for (int slot=0;slot<8;++slot) {
        const auto& w=trace.workers[g][slot];
        if(!w.end) continue;
        if(comma) std::cout << ',';
        comma=true;
        std::cout << "{\"group\":" << g << ",\"slot\":" << slot
            << ",\"begin_us\":" << w.begin << ",\"end_us\":" << w.end
            << ",\"seed_us\":" << w.seed << ",\"setup_us\":" << w.setup
            << ",\"fill_us\":" << w.fills << ",\"fills\":" << w.fillCount
            << ",\"chunks\":[";
        for (std::size_t k=0;k<w.chunks.size();++k) {
            if(k) std::cout << ',';
            const auto& c=w.chunks[k];
            std::cout << '[' << c.first << ',' << c.last << ',' << c.begin << ',' << c.end << ']';
        }
        std::cout << ']';
#ifdef MATCHED_COUNTERS
        std::cout << ",\"stats\":[";
        for (int i=0;i<3;++i) {
            if(!w.stats[i].consistent()) throw std::runtime_error("Phi counter accounting failed");
            if(i) std::cout << ',';
            w.stats[i].writeJson(std::cout);
        }
        std::cout << ']';
#endif
        std::cout << '}';
    }
    std::cout << ']';
}
}
int main(int argc,char** argv) {
    try {
        if(argc!=7) throw std::runtime_error("usage: profile MODE CUTOFF THREADS HISTORY REPEATS EXPECTED");
        const std::string mode(argv[1]),cutoff(argv[2]),history(argv[4]);
        const int threads=std::stoi(argv[3]),repeats=std::stoi(argv[5]);
        const std::int64_t expected=std::stoll(argv[6]);
        if((threads!=1 && threads!=2) || repeats<1 || repeats>50) throw std::runtime_error("invalid limits");
        PhiTables tables(LIMIT);
        Geometry g(cutoff,tables);
        if(history!="none") {
            const bool zero=history=="zero-prefix" || history=="zero-whole";
            if(!zero && history!="meissel-prefix" && history!="meissel-whole" && history!="tail")
                throw std::runtime_error("bad history");
            Geometry prior(zero?"zero":"meissel",tables);
            Matched::Trace trace(prior.split);
            Matched::current=&trace; trace.origin=Matched::Clock::now();
            const auto priorMode=history=="tail"?"tail-phi":(history.find("whole")!=std::string::npos?"baseline":"prefix");
            const auto primed=run(priorMode,prior,threads,tables);
            const auto priorExpected=history=="tail"?1336311383LL:
                (history.find("whole")!=std::string::npos?
                    (zero?11238880589LL:2069687562LL):(zero?9887491223LL:718298196LL));
            if(primed!=priorExpected) throw std::runtime_error("history result mismatch");
        }
        std::cout << std::fixed << std::setprecision(6);
        for(int repeat=0;repeat<repeats;++repeat) {
            Matched::Trace trace(g.split);
            Matched::current=&trace; trace.origin=Matched::Clock::now();
            const auto result=run(mode,g,threads,tables);
            const double elapsed=Matched::now();
            if(result!=expected) throw std::runtime_error("result mismatch");
            std::cout << "{\"mode\":\"" << mode << "\",\"cutoff\":\"" << cutoff
                << "\",\"threads\":" << threads << ",\"history\":\"" << history
                << "\",\"repeat\":" << repeat << ",\"result\":" << result
                << ",\"elapsed_us\":" << elapsed;
            writeTrace(trace);
            std::cout << "}\n";
        }
        Matched::current=nullptr;
    } catch(const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
