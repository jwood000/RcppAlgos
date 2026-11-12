// Native diagnostic driver. Compile with a generated copy of original P2.cpp.
#include "NumbersUtils/PartialSieve.h"
#include "NumbersUtils/PartialSieveUtils.h"
#include "query_trace.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace QueryProfile {
Trace* current=nullptr;
thread_local int workerSlot=-1;
}
namespace {
std::int64_t integer(const char* value) {
    std::size_t used=0; const std::string text(value);
    const auto result=std::stoll(text,&used);
    if(used!=text.size()) throw std::runtime_error("invalid integer argument");
    return result;
}
void writeTrace(const QueryProfile::Trace& trace) {
    std::cout << ",\"first\":" << trace.first << ",\"last\":" << trace.last
              << ",\"workers\":" << trace.workerCount << ",\"worker_records\":[";
    for(std::size_t slot=0;slot<trace.workerCount;++slot) {
        if(slot) std::cout << ',';
        const auto& w=trace.workers[slot];
        std::cout << "{\"slot\":" << slot << ",\"begin_us\":" << w.begin
                  << ",\"end_us\":" << w.end << ",\"chunks\":" << w.chunks
                  << ",\"queries\":" << w.queries << ",\"final_first\":" << w.finalFirst
                  << ",\"final_last\":" << w.finalLast << '}';
    }
    std::cout << "],\"chunks\":[";
    bool comma=false;
    for(const auto& c:trace.chunks) {
        if(c.worker<0) continue;
        if(comma) std::cout << ','; comma=true;
        std::cout << "{\"first\":" << c.first << ",\"last\":" << c.last
                  << ",\"worker\":" << c.worker << ",\"begin_us\":" << c.begin
                  << ",\"end_us\":" << c.end << '}';
    }
    std::cout << "],\"queries\":[";
    comma=false;
    for(const auto& q:trace.queries) {
        if(q.worker<0) continue;
        if(comma) std::cout << ','; comma=true;
        std::cout << "{\"index\":" << q.index << ",\"prime\":" << q.prime
                  << ",\"quotient\":" << q.quotient << ",\"cutoff\":" << q.cutoff
                  << ",\"pi\":" << q.pi << ",\"term\":" << q.term
                  << ",\"worker\":" << q.worker << ",\"begin_us\":" << q.begin
                  << ",\"end_us\":" << q.end << '}';
    }
    std::cout << ']';
}
}
int main(int argc,char** argv) {
    try {
        if(argc!=8) throw std::runtime_error("usage: query-profile CUTOFF X THREADS REPEATS EXPECTED TABLE_NUM TABLE_DEN");
        const std::string cutoff(argv[1]);
        const auto x=integer(argv[2]),threads=integer(argv[3]),repeats=integer(argv[4]);
        const auto expected=integer(argv[5]),numerator=integer(argv[6]),denominator=integer(argv[7]);
        if((cutoff!="zero" && cutoff!="meissel" && cutoff!="lehmer") ||
            x<2 || x>PrimeCounting::PartialSieveDetail::MAX_X || threads<1 || threads>32 ||
            repeats<1 || repeats>32 || expected<0 || expected>x ||
            denominator<1 || numerator<denominator || numerator>std::numeric_limits<int>::max())
            throw std::runtime_error("invalid diagnostic arguments");
        const auto root=PrimeCounting::PartialSieveDetail::IntegerSqrt(x);
        if(root>std::numeric_limits<std::int64_t>::max()/numerator)
            throw std::runtime_error("table ratio overflow");
        const auto limit=root*numerator/denominator;
        PrimeCounting::PhiTables tables(limit);
        const auto a=cutoff=="zero"?0:tables.pi()[cutoff=="meissel"?
            PrimeCounting::PartialSieveDetail::IntegerCbrt(x):
            PrimeCounting::PartialSieveDetail::IntegerSqrt(root)];
        std::cout << std::fixed << std::setprecision(6);
        for(int repeat=0;repeat<repeats;++repeat) {
            QueryProfile::Trace trace;
            QueryProfile::current=&trace;trace.origin=QueryProfile::Clock::now();
            const auto result=PrimeCounting::P2(x,a,static_cast<int>(threads),tables);
            const double elapsed=QueryProfile::now();
            if(result!=expected) throw std::runtime_error("P2 result differs from expected value");
            std::cout << "{\"cutoff\":\"" << cutoff << "\",\"x\":" << x
                      << ",\"a\":" << a << ",\"threads\":" << threads
                      << ",\"table_limit\":" << limit << ",\"table_numerator\":" << numerator
                      << ",\"table_denominator\":" << denominator << ",\"repeat\":" << repeat
                      << ",\"expected\":" << expected << ",\"result\":" << result
                      << ",\"elapsed_us\":" << elapsed;
            writeTrace(trace);std::cout << "}\n";
        }
        QueryProfile::current=nullptr;
    } catch(const std::exception& error) {
        std::cerr << error.what() << '\n';return 1;
    }
}
