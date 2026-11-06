// Serial warm-P2 reuse analysis. Link the copies produced by
// instrument_reuse.py, production Eratosthenes.cpp, and include reuse_metrics.h.
// This replays logical work, not execution time or a production cache design.
#include "NumbersUtils/PartialSieve.h"
#include "NumbersUtils/PartialSieveUtils.h"
#include "NumbersUtils/Phi.h"
#include "reuse_metrics.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
    using Record = P2ReuseMetrics::Record;
    constexpr std::uint32_t MAX_EVENTS = 12000000;
    static_assert(sizeof(Record) == 32, "reuse memory budget assumes 32-byte records");

    void Require(bool condition, const char* message) {
        if (!condition) throw std::runtime_error(message);
    }

    std::int64_t Integer(const char* text) {
        const std::string value(text);
        std::size_t used = 0;
        const auto result = std::stoll(value, &used);
        Require(used == value.size(), "invalid integer argument");
        return result;
    }

    bool IsTop(const std::vector<Record>& records, std::uint32_t index) {
        return index == 0 || records[index - 1].query != records[index].query;
    }

    struct Tally {
        std::uint64_t visits = 0, firstKeys = 0, sameQuery = 0, earlierQuery = 0;
        std::uint64_t topVisits = 0, descendantVisits = 0;
        std::uint64_t topFirstKeys = 0, descendantFirstKeys = 0;
        std::uint64_t kernels = 0, suffix = 0;

        void Add(bool first, bool earlier, bool top, const Record& record) {
            ++visits;
            if (first) ++firstKeys;
            else if (earlier) ++earlierQuery;
            else ++sameQuery;
            if (top) {
                ++topVisits; topFirstKeys += first;
                kernels += record.kernel; suffix += record.suffix;
            }
            else { ++descendantVisits; descendantFirstKeys += first; }
        }
    };

    struct Population {
        std::uint64_t keys = 0, visits = 0;
        void Add(std::uint64_t frequency) { ++keys; visits += frequency; }
    };

    struct Distribution {
        Tally total;
        std::array<Tally, 8> bins;
        std::array<Tally, 2> regions;
        std::array<Population, 3> reasons;
        std::array<Population, 6> aBands;
        std::array<Population, 6> xBands;
        std::map<std::uint32_t, std::uint64_t> frequency;
        std::uint32_t uniqueKeys = 0;
        std::uint64_t indexBytes = 0;
        std::uint64_t withinQueryDuplicates = 0, crossQueryFirstRepeats = 0;
        std::uint32_t repeatedKeys = 0;
        std::uint64_t repeatedMinX = std::numeric_limits<std::uint64_t>::max(), repeatedMaxX = 0;
        std::uint32_t repeatedMinA = std::numeric_limits<std::uint32_t>::max(), repeatedMaxA = 0;
        std::uint32_t maxRepeatedExpanded = 0, maxRepeatedKernels = 0, maxRepeatedSuffix = 0;
    };

    unsigned ABand(std::uint32_t a) {
        return a < 32 ? 0 : a < 64 ? 1 : a < 100 ? 2 : a < 256 ? 3 : a < 1024 ? 4 : 5;
    }

    unsigned XBand(std::uint64_t x) {
        return x <= 65535 ? 0 : x <= 262143 ? 1 : x <= 1048575 ? 2 :
            x <= 16777215 ? 3 : x <= 4294967295ULL ? 4 : 5;
    }

    unsigned Bin(std::uint32_t query, std::uint32_t first, std::uint32_t last) {
        Require(first < last && query >= first && query < last, "query outside trace interval");
        return static_cast<unsigned>((std::uint64_t(query - first) * 8) / (last - first));
    }

    Distribution AssignKeys(std::vector<Record>& records, std::uint32_t first,
                            std::uint32_t last, std::uint32_t commonFirst) {
        Distribution result;
        std::vector<std::uint32_t> order(records.size());
        std::iota(order.begin(), order.end(), std::uint32_t(0));
        result.indexBytes = order.capacity() * sizeof(std::uint32_t);
        std::sort(order.begin(), order.end(), [&](std::uint32_t left, std::uint32_t right) {
            const auto& l = records[left];
            const auto& r = records[right];
            if (l.x != r.x) return l.x < r.x;
            if (l.a != r.a) return l.a < r.a;
            if (l.query != r.query) return l.query < r.query;
            return left < right;
        });
        for (std::size_t begin = 0; begin < order.size();) {
            const auto& key = records[order[begin]];
            std::size_t end = begin + 1;
            while (end < order.size() && records[order[end]].x == key.x &&
                   records[order[end]].a == key.a) ++end;
            const auto frequency = static_cast<std::uint32_t>(end - begin);
            ++result.frequency[frequency];
            const unsigned reason = key.x > 65535 ? (key.a >= 100 ? 2 : 0) : 1;
            result.reasons[reason].Add(frequency);
            result.aBands[ABand(key.a)].Add(frequency);
            result.xBands[XBand(key.x)].Add(frequency);
            if (frequency > 1) {
                ++result.repeatedKeys;
                result.repeatedMinX = std::min(result.repeatedMinX, key.x);
                result.repeatedMaxX = std::max(result.repeatedMaxX, key.x);
                result.repeatedMinA = std::min(result.repeatedMinA, key.a);
                result.repeatedMaxA = std::max(result.repeatedMaxA, key.a);
            }
            const auto firstQuery = key.query;
            std::uint32_t distinctQueries = 0;
            for (std::size_t position = begin; position < end; ++position) {
                const auto index = order[position];
                auto& record = records[index];
                record.keyId = result.uniqueKeys;
                if (frequency > 1) {
                    result.maxRepeatedExpanded = std::max(result.maxRepeatedExpanded, record.end - index);
                    result.maxRepeatedKernels = std::max(result.maxRepeatedKernels, record.kernel);
                    result.maxRepeatedSuffix = std::max(result.maxRepeatedSuffix, record.suffix);
                }
                const bool firstVisit = position == begin;
                const bool earlier = record.query > firstQuery;
                const bool top = IsTop(records, index);
                if (position == begin || record.query != records[order[position - 1]].query) {
                    ++distinctQueries;
                }
                result.total.Add(firstVisit, earlier, top, record);
                result.bins[Bin(record.query, first, last)].Add(firstVisit, earlier, top, record);
                result.regions[record.query < commonFirst ? 0 : 1].Add(firstVisit, earlier, top, record);
            }
            result.withinQueryDuplicates += frequency - distinctQueries;
            result.crossQueryFirstRepeats += distinctQueries - 1;
            ++result.uniqueKeys;
            begin = end;
        }
        // The sorting vector is destroyed before any replay's key-state array.
        return result;
    }

    struct Policy {
        const char* name;
        std::uint64_t maxX;
        std::uint32_t exclusiveA;
        bool Accept(const Record& record) const {
            return record.x <= maxX && record.a < exclusiveA;
        }
    };

    struct Savings {
        std::uint64_t hits = 0, topHits = 0, descendantHits = 0;
        std::uint64_t expanded = 0, kernels = 0, suffix = 0, retained = 0;
        void Hit(const Record& record, std::uint32_t index, bool top) {
            ++hits;
            if (top) ++topHits;
            else ++descendantHits;
            expanded += record.end - index;
            kernels += record.kernel;
            suffix += record.suffix;
        }
    };

    struct Replay {
        Policy policy;
        std::uint64_t candidates = 0;
        std::uint32_t candidatePages = 0;
        bool hasPageCount = false;
        Savings total;
        std::array<Savings, 2> regions;
        explicit Replay(const Policy& value) : policy(value) {}
    };

    Replay ReplayPolicy(const std::vector<Record>& records, std::uint32_t keys,
                         const Policy& policy, std::uint32_t commonFirst) {
        Replay result(policy);
        std::vector<unsigned char> available(keys, 0);
        result.hasPageCount = policy.maxX <= 1048575 && policy.exclusiveA <= 256;
        const std::size_t pagesPerRow = result.hasPageCount ? policy.maxX / 4096 + 1 : 0;
        std::vector<unsigned char> pages(pagesPerRow *
            (result.hasPageCount ? policy.exclusiveA : 0), 0);
        // First count full-trace candidate keys. Some are never populated in
        // replay because a parent hit removes their first baseline appearance.
        for (const auto& record : records) {
            if (policy.Accept(record) && !available[record.keyId]) {
                available[record.keyId] = 1;
                ++result.candidates;
                if (result.hasPageCount) {
                    const auto page = static_cast<std::size_t>(record.a) * pagesPerRow + record.x / 4096;
                    if (!pages[page]) { pages[page] = 1; ++result.candidatePages; }
                }
            }
        }
        std::fill(available.begin(), available.end(), 0);
        std::vector<std::uint32_t> pending;
        std::uint32_t index = 0;
        while (index <= records.size()) {
            // Publish only when the original subtree has completed. In
            // particular, an active ancestor is never an available cache hit.
            while (!pending.empty() && records[pending.back()].end <= index) {
                const auto& completed = records[pending.back()];
                if (!available[completed.keyId]) {
                    available[completed.keyId] = 1;
                    ++result.total.retained;
                    ++result.regions[completed.query < commonFirst ? 0 : 1].retained;
                }
                pending.pop_back();
            }
            if (index == records.size()) break;
            const auto& record = records[index];
            if (policy.Accept(record)) {
                if (available[record.keyId]) {
                    const bool top = IsTop(records, index);
                    result.total.Hit(record, index, top);
                    result.regions[record.query < commonFirst ? 0 : 1].Hit(record, index, top);
                    // This single jump accounts for the whole subtree once;
                    // none of its descendants can become additional hits.
                    index = record.end;
                    continue;
                }
                pending.push_back(index);
            }
            ++index;
        }
        Require(pending.empty(), "unfinished replay completion stack");
        return result;
    }

    Record Node(std::uint64_t x, std::uint32_t query, std::uint32_t end,
                std::uint32_t kernel, std::uint32_t suffix) {
        Record record{};
        record.x = x; record.a = 7; record.query = query; record.end = end;
        record.kernel = kernel; record.suffix = suffix; record.keyId = 0;
        return record;
    }

    void SelfTest() {
        const Policy all = {"all", std::numeric_limits<std::uint64_t>::max(),
                             std::numeric_limits<std::uint32_t>::max()};
        // B is cached within the first query. The second A hit removes a
        // subtree containing another B, which must not be double-counted.
        std::vector<Record> records = {
            Node(900000, 1, 4, 6, 10), Node(700000, 1, 2, 2, 3),
            Node(700000, 1, 3, 2, 3), Node(800000, 1, 4, 1, 4),
            Node(900000, 2, 7, 4, 7), Node(700000, 2, 6, 2, 3),
            Node(800000, 2, 7, 1, 4)
        };
        const auto distribution = AssignKeys(records, 1, 3, 2);
        const auto replay = ReplayPolicy(records, distribution.uniqueKeys, all, 2);
        Require(distribution.uniqueKeys == 3 && distribution.total.sameQuery == 1 &&
                distribution.total.earlierQuery == 3 && distribution.withinQueryDuplicates == 1 &&
                distribution.crossQueryFirstRepeats == 3, "synthetic frequency categorization failed");
        Require(replay.total.hits == 2 && replay.total.expanded == 4 &&
                replay.total.kernels == 6 && replay.total.suffix == 10 &&
                replay.total.retained == 3 && replay.total.topHits == 1,
                "synthetic overlapping subtree replay failed");

        // Deliberately allow a repeated active key in this artificial trace:
        // it cannot hit before a completed occurrence has published a value.
        records = {Node(900000, 1, 2, 2, 0), Node(900000, 1, 2, 1, 0),
                   Node(900000, 2, 3, 1, 0)};
        const auto activeKeys = AssignKeys(records, 1, 3, 2);
        const auto active = ReplayPolicy(records, activeKeys.uniqueKeys, all, 2);
        Require(active.total.hits == 1 && active.total.kernels == 1 &&
                active.total.retained == 1, "synthetic completion publication failed");
    }

    void Validate(const P2ReuseMetrics::State& trace, const PrimeCounting::PhiTables& tables,
                  std::int64_t x) {
        const auto& records = trace.records;
        Require(trace.first >= 0 && trace.last >= trace.first &&
                trace.last <= std::numeric_limits<std::uint32_t>::max(), "invalid query interval");
        std::vector<std::uint32_t> ancestors;
        std::uint64_t topKernels = 0, topSuffix = 0, topCount = 0;
        auto nextQuery = trace.first;
        for (std::uint32_t index = 0; index < records.size(); ++index) {
            const auto& record = records[index];
            Require(record.a > 6 && (record.x > 65535 || record.a >= 100),
                    "trace contains a currently cache-eligible or tiny key");
            Require(record.end > index && record.end <= records.size() &&
                    record.kernel >= record.end - index, "invalid subtree interval/count");
            while (!ancestors.empty() && records[ancestors.back()].end == index) ancestors.pop_back();
            if (!ancestors.empty()) {
                const auto& parent = records[ancestors.back()];
                Require(record.end <= parent.end && record.query == parent.query,
                        "subtree intervals are not properly nested");
            }
            if (IsTop(records, index)) {
                Require(ancestors.empty() && record.query == nextQuery,
                        "missing/nonmonotonic query root in trace");
                const auto quotient = x / tables.primes()[record.query];
                const auto cutoff = tables.pi()[PrimeCounting::PartialSieveDetail::IntegerSqrt(quotient)];
                Require(record.x == static_cast<std::uint64_t>(quotient) &&
                        record.a == static_cast<std::uint32_t>(cutoff),
                        "first recorded node is not the public Phi query root");
                ++nextQuery; ++topCount;
                topKernels += record.kernel; topSuffix += record.suffix;
            } else {
                Require(!ancestors.empty(), "disconnected subtree within a query");
            }
            ancestors.push_back(index);
        }
        Require(nextQuery == trace.last && topCount == trace.queries &&
                topKernels == trace.totalKernels && topSuffix == trace.totalSuffix,
                "trace roots do not cover every query and all kernel/suffix work");
    }

    void PrintTally(const Tally& value) {
        std::cout << "{\"node_visits\":" << value.visits << ",\"unique_first_appearances\":" << value.firstKeys
                  << ",\"repeats_same_first_query\":" << value.sameQuery
                  << ",\"repeats_from_earlier_queries\":" << value.earlierQuery
                  << ",\"top_visits\":" << value.topVisits << ",\"descendant_visits\":" << value.descendantVisits
                  << ",\"top_first_keys\":" << value.topFirstKeys
                  << ",\"descendant_first_keys\":" << value.descendantFirstKeys
                  << ",\"kernel_entries\":" << value.kernels
                  << ",\"suffix_terms\":" << value.suffix << '}';
    }

    template <std::size_t Size>
    void PrintBands(const std::array<Population, Size>& values, const char* const (&names)[Size]) {
        std::cout << '[';
        for (std::size_t i = 0; i < Size; ++i) {
            if (i) std::cout << ',';
            std::cout << "{\"band\":\"" << names[i] << "\",\"keys\":" << values[i].keys
                      << ",\"node_visits\":" << values[i].visits << '}';
        }
        std::cout << ']';
    }

    void PrintSavings(const Savings& value) {
        std::cout << "{\"hits\":" << value.hits << ",\"top_hits\":" << value.topHits
                  << ",\"descendant_hits\":" << value.descendantHits
                  << ",\"saved_expanded\":" << value.expanded
                  << ",\"saved_kernels_including_hit\":" << value.kernels
                  << ",\"saved_kernel_descendants\":" << value.kernels - value.hits
                  << ",\"saved_suffix\":" << value.suffix << ",\"retained_keys\":" << value.retained
                  << ",\"logical_bytes_packed20\":" << value.retained * 20
                  << ",\"logical_bytes_aligned24\":" << value.retained * 24 << '}';
    }
}

int main(int argc, char** argv) {
    try {
        SelfTest();
        if (argc == 2 && std::string(argv[1]) == "--self-test") {
            std::cout << "{\"self_test\":true,\"cases\":2}\n";
            return 0;
        }
        Require(argc == 5, "usage: reuse meissel|lehmer x expectedP2 maxEvents");
        const std::string method(argv[1]);
        const auto x = Integer(argv[2]);
        const auto expected = Integer(argv[3]);
        const auto maximum = Integer(argv[4]);
        Require((method == "meissel" || method == "lehmer") && x >= 2 &&
                x <= PrimeCounting::PartialSieveDetail::MAX_X && expected >= 0 && expected <= x &&
                maximum >= 1 && maximum <= MAX_EVENTS, "invalid reuse arguments");
        std::int64_t a = 0, cubeA = 0, result = 0;
        {
            const PrimeCounting::PhiTables tables(PrimeCounting::PartialSieveDetail::IntegerSqrt(x));
            cubeA = tables.pi()[PrimeCounting::PartialSieveDetail::IntegerCbrt(x)];
            a = method == "meissel" ? cubeA :
                tables.pi()[PrimeCounting::PartialSieveDetail::IntegerSqrt(tables.limit())];
            Require(PrimeCounting::P2(x, a, 1, tables) == expected, "warm-up P2 mismatch");
            P2ReuseMetrics::begin(static_cast<std::size_t>(maximum));
            result = PrimeCounting::P2(x, a, 1, tables);
            P2ReuseMetrics::end();
            Require(result == expected, "recorded P2 mismatch");
            Validate(P2ReuseMetrics::state(), tables, x);
        }
        // Release the dense tables before allocating the sorting index.
        auto& trace = P2ReuseMetrics::state();
        const auto first = static_cast<std::uint32_t>(trace.first);
        const auto last = static_cast<std::uint32_t>(trace.last);
        const auto commonFirst = static_cast<std::uint32_t>(std::min<std::int64_t>(cubeA + 1, last));
        const auto distribution = AssignKeys(trace.records, first, last, commonFirst);
        const std::uint64_t unlimitedX = PrimeCounting::PartialSieveDetail::MAX_X;
        const auto unlimitedA = std::numeric_limits<std::uint32_t>::max();
        const Policy policies[] = {
            {"all_ineligible", unlimitedX, unlimitedA},
            {"a_lt_100", unlimitedX, 100}, {"a_lt_256", unlimitedX, 256},
            {"x_le_131071_a_lt_100", 131071, 100},
            {"x_le_262143_a_lt_64", 262143, 64},
            {"x_le_262143_a_lt_100", 262143, 100},
            {"x_le_1048575_a_lt_100", 1048575, 100},
            {"x_le_1048575_a_lt_256", 1048575, 256},
            {"rows_only_x_le_65535_a_lt_256", 65535, 256}
        };
        std::vector<Replay> replays;
        for (const auto& policy : policies) {
            replays.push_back(ReplayPolicy(trace.records, distribution.uniqueKeys, policy, commonFirst));
            const auto& saved = replays.back().total;
            Require(saved.expanded <= trace.records.size() && saved.kernels <= trace.totalKernels &&
                    saved.suffix <= trace.totalSuffix, "replay double-counted baseline work");
        }
        std::cout << "{\"method\":\"" << method << "\",\"x\":" << x << ",\"a\":" << a
                  << ",\"p2\":" << result << ",\"queries\":" << trace.queries
                  << ",\"first_query\":" << first << ",\"last_query_exclusive\":" << last
                  << ",\"common_first_query\":" << commonFirst
                  << ",\"extra_prefix_queries\":" << commonFirst - first
                  << ",\"total_kernels\":" << trace.totalKernels << ",\"total_suffix\":" << trace.totalSuffix
                  << ",\"records\":" << trace.records.size() << ",\"unique_keys\":" << distribution.uniqueKeys
                  << ",\"within_query_duplicates\":" << distribution.withinQueryDuplicates
                  << ",\"cross_query_first_repeats\":" << distribution.crossQueryFirstRepeats
                  << ",\"repeated_key_count\":" << distribution.repeatedKeys
                  << ",\"max_repeated_ineligible_subtree\":" << distribution.maxRepeatedExpanded
                  << ",\"max_repeated_kernel_subtree\":" << distribution.maxRepeatedKernels
                  << ",\"max_repeated_suffix_subtree\":" << distribution.maxRepeatedSuffix
                  << ",\"repeated_key_min_x\":";
        if (distribution.repeatedKeys) std::cout << distribution.repeatedMinX;
        else std::cout << "null";
        std::cout << ",\"repeated_key_max_x\":";
        if (distribution.repeatedKeys) std::cout << distribution.repeatedMaxX;
        else std::cout << "null";
        std::cout << ",\"repeated_key_min_a\":";
        if (distribution.repeatedKeys) std::cout << distribution.repeatedMinA;
        else std::cout << "null";
        std::cout << ",\"repeated_key_max_a\":";
        if (distribution.repeatedKeys) std::cout << distribution.repeatedMaxA;
        else std::cout << "null";
        std::cout
                  << ",\"trace_bytes\":" << trace.records.capacity() * sizeof(Record)
                  << ",\"index_bytes\":" << distribution.indexBytes << ",\"totals\":";
        PrintTally(distribution.total);
        std::cout << ",\"bins\":[";
        for (std::size_t bin = 0; bin < 8; ++bin) {
            if (bin) std::cout << ',';
            std::cout << "{\"bin\":" << bin << ",\"counts\":";
            PrintTally(distribution.bins[bin]); std::cout << '}';
        }
        std::cout << "],\"regions\":{\"extra_prefix\":";
        PrintTally(distribution.regions[0]); std::cout << ",\"common_suffix\":";
        PrintTally(distribution.regions[1]); std::cout << "},\"eligibility_reasons\":";
        const char* const reasons[] = {"x_only", "a_only", "both"};
        PrintBands(distribution.reasons, reasons);
        std::cout << ",\"a_histogram\":";
        const char* const aNames[] = {"7..31", "32..63", "64..99", "100..255", "256..1023", "1024+"};
        PrintBands(distribution.aBands, aNames);
        std::cout << ",\"x_histogram\":";
        const char* const xNames[] = {"0..65535", "65536..262143", "262144..1048575", "1048576..16777215",
                                    "16777216..4294967295", "4294967296+"};
        PrintBands(distribution.xBands, xNames);
        std::cout << ",\"frequency_histogram\":[";
        bool comma = false;
        for (const auto& entry : distribution.frequency) {
            if (comma) std::cout << ',';
            comma = true;
            std::cout << "{\"frequency\":" << entry.first << ",\"keys\":" << entry.second << '}';
        }
        std::cout << "],\"policies\":[";
        for (std::size_t i = 0; i < replays.size(); ++i) {
            if (i) std::cout << ',';
            const auto& replay = replays[i];
            std::cout << "{\"name\":\"" << replay.policy.name << "\",\"max_x\":" << replay.policy.maxX
                      << ",\"a_exclusive\":" << replay.policy.exclusiveA
                      << ",\"candidate_keys\":" << replay.candidates << ",\"candidate_pages4096\":";
            if (replay.hasPageCount) std::cout << replay.candidatePages;
            else std::cout << "null";
            std::cout << ",\"savings\":";
            PrintSavings(replay.total);
            std::cout << ",\"extra_prefix\":"; PrintSavings(replay.regions[0]);
            std::cout << ",\"common_suffix\":"; PrintSavings(replay.regions[1]);
            std::cout << '}';
        }
        std::cout << "]}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
