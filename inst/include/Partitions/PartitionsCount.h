#pragma once

#include "Partitions/PartitionsTypes.h"
#include <memory>

void PartitionsCount(const std::vector<int> &Reps,
                     PartDesign &part, int lenV, bool bIsCount);

class CountClass {
protected:
    std::vector<mpz_class> p1;
    std::vector<mpz_class> p2;
    int size = 0;

public:
    virtual ~CountClass() = default;
    void SetArrSize(PartitionType ptype, int n, int m, int cap);

    virtual double GetCount(int n, int m, int cap, int strtLen)  = 0;
    virtual void GetCount(mpz_class &res, int n, int m, int cap,
                          int strtLen, bool bLiteral = true) = 0;
    void InitializeMpz();
};

std::unique_ptr<CountClass> MakeCount(PartitionType ptype,
                                      bool IsComp = false);

class DistinctAll : public CountClass {
    double GetCount(int n, int m, int cap, int strtLen);
    void GetCount(mpz_class &res, int n, int m, int cap,
                  int strtLen, bool bLiteral = true);
};

class DistinctLen : public CountClass {
    double GetCount(int n, int m, int cap, int strtLen);
    void GetCount(mpz_class &res, int n, int m, int cap,
                  int strtLen, bool bLiteral = true);
};

class DistinctLenCap : public CountClass {
    double GetCount(int n, int m, int cap, int strtLen);
    void GetCount(mpz_class &res, int n, int m, int cap,
                  int strtLen, bool bLiteral = true);
};

class DistinctMZ : public CountClass {
    double GetCount(int n, int m, int cap, int strtLen);
    void GetCount(mpz_class &res, int n, int m, int cap,
                  int strtLen, bool bLiteral = true);
};

class DistinctCapMZ : public CountClass {
    double GetCount(int n, int m, int cap, int strtLen);
    void GetCount(mpz_class &res, int n, int m, int cap,
                  int strtLen, bool bLiteral = true);
};

class RepAll : public CountClass {
    double GetCount(int n, int m, int cap, int strtLen);
    void GetCount(mpz_class &res, int n, int m, int cap,
                  int strtLen, bool bLiteral = true);
};

class RepLen : public CountClass {
    double GetCount(int n, int m, int cap, int strtLen);
    void GetCount(mpz_class &res, int n, int m, int cap,
                  int strtLen, bool bLiteral = true);
};

class RepLenCap : public CountClass {
    double GetCount(int n, int m, int cap, int strtLen);
    void GetCount(mpz_class &res, int n, int m, int cap,
                  int strtLen, bool bLiteral = true);
};

class CompsRepLen : public CountClass {
    double GetCount(int n, int m, int cap, int strtLen);
    void GetCount(mpz_class &res, int n, int m, int cap,
                  int strtLen, bool bLiteral = true);
};

class CompsRepZero : public CountClass {
    double GetCount(int n, int m, int cap, int strtLen);
    void GetCount(mpz_class &res, int n, int m, int cap,
                  int strtLen, bool bLiteral = true);
};
