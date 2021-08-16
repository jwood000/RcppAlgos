#ifndef PARTITIONS_COUNT_H
#define PARTITIONS_COUNT_H

#include "Partitions/PartitionsTypes.h"
#include "Cpp14MakeUnique.h"

void PartitionsCount(const std::vector<int> &Reps, PartDesign &part,
                     int lenV, bool bCalcDifficult, bool IsComb);

class CountClass {
protected:
    std::unique_ptr<mpz_t[]> p1;
    std::unique_ptr<mpz_t[]> p2;
    int size = 0;

public:
    virtual ~CountClass() = default;
    void SetArrSize(PartitionType ptype, int n, int m, int cap);
    static std::unique_ptr<CountClass> MakeCount(PartitionType ptype);

    virtual double GetCount(int n, int m, int cap, int strtLen) {
        return 0.0;
    };

    virtual void GetCount(mpz_t res, int n, int m, int cap,
                          int strtLen, bool bLiteral = true) {};

    void InitializeMpz();
    void ClearMpz();
};

class DistinctAll : public CountClass {
public:
    double GetCount(int n, int m, int cap, int strtLen);
    void GetCount(mpz_t res, int n, int m, int cap,
                  int strtLen, bool bLiteral = true);
};

class DistinctLen : public CountClass {
public:
    double GetCount(int n, int m, int cap, int strtLen);
    void GetCount(mpz_t res, int n, int m, int cap,
                  int strtLen, bool bLiteral = true);
};

class DistinctLenCap : public CountClass {
public:
    double GetCount(int n, int m, int cap, int strtLen);
    void GetCount(mpz_t res, int n, int m, int cap,
                  int strtLen, bool bLiteral = true);
};

class DistinctMZ : public CountClass {
public:
    double GetCount(int n, int m, int cap, int strtLen);
    void GetCount(mpz_t res, int n, int m, int cap,
                  int strtLen, bool bLiteral = true);
};

class DistinctCapMZ : public CountClass {
public:
    double GetCount(int n, int m, int cap, int strtLen);
    void GetCount(mpz_t res, int n, int m, int cap,
                  int strtLen, bool bLiteral = true);
};

class RepAll : public CountClass {
public:
    double GetCount(int n, int m, int cap, int strtLen);
    void GetCount(mpz_t res, int n, int m, int cap,
                  int strtLen, bool bLiteral = true);
};

class RepLen : public CountClass {
public:
    double GetCount(int n, int m, int cap, int strtLen);
    void GetCount(mpz_t res, int n, int m, int cap,
                  int strtLen, bool bLiteral = true);
};

class RepLenCap : public CountClass {
public:
    double GetCount(int n, int m, int cap, int strtLen);
    void GetCount(mpz_t res, int n, int m, int cap,
                  int strtLen, bool bLiteral = true);
};

#endif
