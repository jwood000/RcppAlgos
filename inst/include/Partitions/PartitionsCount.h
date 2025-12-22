#pragma once

#include "Partitions/PartitionsTypes.h"
#include <memory>

int PartitionsCount(const std::vector<int> &Reps,
                    PartDesign &part, int lenV, bool bIsCount);

class CountClass {
protected:
    std::vector<mpz_class> p1;
    std::vector<mpz_class> p2;
    std::vector<std::vector<mpz_class>> p2d;

    int size = 0;
    int width = 0;

public:
    virtual ~CountClass() = default;
    void SetArrSize(PartitionType ptype, int n, int m);

    virtual double GetCount(
        int n, int m, const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0
    ) = 0;

    virtual void GetCount(
        mpz_class &res, int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0, bool bLiteral = true
    ) = 0;

    void InitializeMpz();
};

std::unique_ptr<CountClass> MakeCount(PartitionType ptype);

class DistinctAll : public CountClass {
    double GetCount(
        int n, int m, const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0
    );

    void GetCount(
        mpz_class &res, int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0, bool bLiteral = true
    );
};

class DistinctLen : public CountClass {
    double GetCount(
        int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0
    );

    void GetCount(
        mpz_class &res, int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0, bool bLiteral = true
    );
};

class DistinctLenRstrctd : public CountClass {
    double GetCount(
        int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0
    );

    void GetCount(
        mpz_class &res, int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0, bool bLiteral = true
    );
};

class DistinctMZ : public CountClass {
    double GetCount(
        int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0
    );

    void GetCount(
        mpz_class &res, int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0, bool bLiteral = true
    );
};

class DistinctRstrctdMZ : public CountClass {
    double GetCount(
        int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0
    );

    void GetCount(
        mpz_class &res, int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0, bool bLiteral = true
    );
};

class RepAll : public CountClass {
    double GetCount(
        int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0
    );

    void GetCount(
        mpz_class &res, int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0, bool bLiteral = true
    );
};

class RepLen : public CountClass {
    double GetCount(
        int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0
    );

    void GetCount(
        mpz_class &res, int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0, bool bLiteral = true
    );
};

class RepLenRstrctd : public CountClass {
    double GetCount(
        int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0
    );

    void GetCount(
        mpz_class &res, int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0, bool bLiteral = true
    );
};

class PermDstnctRstrctd : public CountClass {
    double GetCount(
        int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0
    );

    void GetCount(
        mpz_class &res, int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0, bool bLiteral = true
    );
};

class PermDstnctRstrctdMZ : public CountClass {
    double GetCount(
        int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0
    );

    void GetCount(
        mpz_class &res, int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0, bool bLiteral = true
    );
};

class CompsRepLen : public CountClass {
    double GetCount(
        int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0
    );

    void GetCount(
        mpz_class &res, int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0, bool bLiteral = true
    );
};

class CompsRepZero : public CountClass {
    double GetCount(
        int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0
    );

    void GetCount(
        mpz_class &res, int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0, bool bLiteral = true
    );
};

class CompsDistinctLen : public CountClass {
    double GetCount(
        int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0
    );

    void GetCount(
        mpz_class &res, int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0, bool bLiteral = true
    );
};

class CompsDistLenMZ : public CountClass {
    double GetCount(
        int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0
    );

    void GetCount(
        mpz_class &res, int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0, bool bLiteral = true
    );
};

class CompsDistLenMZWeak : public CountClass {
    double GetCount(
        int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0
    );

    void GetCount(
        mpz_class &res, int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0, bool bLiteral = true
    );
};

class CompsDstnctRstrctdMZ : public CountClass {
    double GetCount(
        int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0
    );

    void GetCount(
        mpz_class &res, int n, int m,
        const std::vector<int> &allowed = std::vector<int>(),
        int strtLen = 0, bool bLiteral = true
    );
};
