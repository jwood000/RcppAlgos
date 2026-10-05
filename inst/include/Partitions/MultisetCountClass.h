#pragma once

#include "Partitions/PartitionsCount.h"

class MultisetCountClass : public CountClass {
protected:
    const std::vector<int> Reps;

public:
    ~MultisetCountClass() override = default;

    explicit MultisetCountClass(const std::vector<int>& reps)
        : Reps(reps) {}

    double GetCount(
        int n, int m,
        const std::vector<int>& allowed = std::vector<int>(),
        int strtLen = 0
    ) override = 0;

    void GetCount(
        mpz_class& res, int n, int m,
        const std::vector<int>& allowed = std::vector<int>(),
        int strtLen = 0, bool bLiteral = true
    ) override = 0;
};

std::unique_ptr<MultisetCountClass> MakeMultisetCount(
    PartitionType ptype, const std::vector<int>& reps
);

class PartsMultiset : public MultisetCountClass {
public:
    explicit PartsMultiset(const std::vector<int>& reps)
        : MultisetCountClass(reps) {}

    double GetCount(
        int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0
    ) override;

    void GetCount(
        mpz_class &res, int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0, bool bLiteral = true
    ) override;
};

class CompsMultiset : public MultisetCountClass {
public:
    explicit CompsMultiset(const std::vector<int>& reps)
        : MultisetCountClass(reps) {}

    double GetCount(
        int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0
    ) override;

    void GetCount(
        mpz_class &res, int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0, bool bLiteral = true
    ) override;
};

class CompsMultisetZero : public MultisetCountClass {
public:
    explicit CompsMultisetZero(const std::vector<int>& reps)
        : MultisetCountClass(reps) {}

    double GetCount(
        int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0
    ) override;

    void GetCount(
        mpz_class &res, int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0, bool bLiteral = true
    ) override;
};

class CompsMultisetWeak : public MultisetCountClass {
public:
    explicit CompsMultisetWeak(const std::vector<int>& reps)
        : MultisetCountClass(reps) {}

    double GetCount(
        int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0
    ) override;

    void GetCount(
        mpz_class &res, int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0, bool bLiteral = true
    ) override;
};
