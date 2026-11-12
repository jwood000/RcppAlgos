#include "InitialChunks.h"
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

namespace {
std::uint64_t checks = 0;
void check(bool ok) { ++checks; if (!ok) throw std::runtime_error("queue check failed"); }
using Range = std::pair<std::int64_t, std::int64_t>;

template <std::int64_t K>
std::vector<Range> expected(std::int64_t first, std::int64_t last) {
    std::vector<Range> out;
    const auto headEnd = first + std::min<std::int64_t>(32, last-first);
    auto begin=first;
    while(begin<headEnd) {
        const auto end=begin+std::min(K,headEnd-begin);
        out.emplace_back(begin,end); begin=end;
    }
    while(begin<last) {
        const auto end=begin+std::min<std::int64_t>(32,last-begin);
        out.emplace_back(begin,end); begin=end;
    }
    return out;
}

template <std::int64_t K>
void serial(std::int64_t first,std::int64_t last) {
    InitialChunks::Queue<K> q(first,last);
    const auto ref=expected<K>(first,last);
    check(q.size()==ref.size());
    for(std::uint64_t i=0;i<ref.size();++i) {
        std::int64_t begin=-1,end=-1;
        check(q.range(i,begin,end));check(Range(begin,end)==ref[i]);
        begin=end=-1;check(q.claim(begin,end));check(Range(begin,end)==ref[i]);
    }
    for(const auto ticket : {q.size(),q.size()+1,std::numeric_limits<std::uint64_t>::max()}) {
        std::int64_t begin=-1,end=-1;
        check(!q.range(ticket,begin,end));check(begin==-1 && end==-1);
    }
    std::int64_t begin=-1,end=-1;check(!q.claim(begin,end));
}

template <std::int64_t K>
void concurrent(std::int64_t first,std::int64_t last,int workers) {
    InitialChunks::Queue<K> q(first,last);
    std::vector<std::vector<Range>> slots(workers);
    std::vector<std::thread> threads;
    for(int slot=0;slot<workers;++slot) threads.emplace_back([&,slot] {
        std::int64_t begin,end;
        while(q.claim(begin,end)) slots[slot].emplace_back(begin,end);
    });
    for(auto& thread:threads) thread.join();
    std::vector<Range> actual;
    for(const auto& slot:slots) actual.insert(actual.end(),slot.begin(),slot.end());
    std::sort(actual.begin(),actual.end());
    check(actual==expected<K>(first,last));
}

template <std::int64_t K>
void noBarrier() {
    // A worker holds the initial claim unfinished. The other worker must be
    // able to claim tail work without a completion handshake from the first.
    InitialChunks::Queue<K> q(1,98);
    std::atomic<bool> held(false),tailClaimed(false);
    bool firstOK=false,secondOK=false;
    std::thread first([&] {
        std::int64_t begin,end;
        firstOK=q.claim(begin,end) && begin==1;
        held.store(true,std::memory_order_release);
        while(!tailClaimed.load(std::memory_order_acquire)) std::this_thread::yield();
    });
    std::thread second([&] {
        while(!held.load(std::memory_order_acquire)) std::this_thread::yield();
        std::int64_t begin,end;
        while(q.claim(begin,end)) if(begin>=33) { secondOK=true;break; }
        tailClaimed.store(true,std::memory_order_release);
    });
    first.join();second.join();check(firstOK && secondOK);
}

template <std::int64_t K>
void run() {
    const auto maximum=std::numeric_limits<std::int64_t>::max();
    for(std::int64_t n=0;n<=130;++n)
        for(const auto first:{std::int64_t(0),std::int64_t(1),std::int64_t(590),maximum-130})
            serial<K>(first,first+n);
    for(const auto n:{2047,2048,4095,4096,8191,8192,8193}) serial<K>(591,591+n);
    for(const auto n:{0,1,2,31,32,33,63,64,65,2047,4095,4096,8193})
        for(const auto workers:{1,2,3,4,8}) concurrent<K>(maximum-8193,maximum-8193+n,workers);
    noBarrier<K>();

    // Allocation-free extreme ranges, including a partial last task.
    for(const auto first:{std::int64_t(0),std::int64_t(1),std::int64_t(31)}) {
        InitialChunks::Queue<K> q(first,maximum);
        for(const auto ticket:{std::uint64_t(0),std::uint64_t(31/K),std::uint64_t(32/K),q.size()-1}) {
            std::int64_t begin=-1,end=-1;
            check(q.range(ticket,begin,end));check(first<=begin && begin<end && end<=maximum);
            if(ticket==q.size()-1) check(end==maximum);
        }
        std::int64_t begin=-1,end=-1;check(!q.range(q.size(),begin,end));
    }
    for(const auto bad:{Range(-1,0),Range(2,1),Range(-2,-1)}) {
        bool caught=false;
        try{InitialChunks::Queue<K> q(bad.first,bad.second);}
        catch(const std::out_of_range&){caught=true;}
        check(caught);
    }
}
}
int main() {
    try {
        run<1>();run<2>();run<4>();run<8>();run<32>();
        std::cout << "PASS queue_checks=" << checks << '\n';
    }catch(const std::exception& e){std::cerr << e.what() << '\n';return 1;}
}
