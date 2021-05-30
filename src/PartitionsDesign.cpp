#include "Partitions/PartitionsTypes.h"
#include <RcppThread.h>
#include <numeric>
#include <string>

std::string GetPartitionType(const PartDesign &part) {

    std::string res;

    switch (part.ptype) {
        case PartitionType::NotPartition: {
            res = "NotPartition";
            break;
        } case PartitionType::DistCapped: {
            res = "DistCapped";
            break;
        } case PartitionType::DstctNoZero : {
            res = "DstctNoZero";
            break;
        } case PartitionType::DstctOneZero: {
            res = "DstctOneZero";
            break;
        } case PartitionType::DstctShort: {
            res = "DstctShort";
            break;
        } case PartitionType::DstctSpecial : {
            res = "DstctSpecial";
            break;
        } case PartitionType::DstctStdAll: {
            res = "DstctStdAll";
            break;
        } case PartitionType::Multiset: {
            res = "Multiset";
            break;
        } case PartitionType::RepCapped : {
            res = "RepCapped";
            break;
        } case PartitionType::RepNoZero: {
            res = "RepNoZero";
            break;
        } case PartitionType::RepShort : {
            res = "RepShort";
            break;
        } default: {
            res = "RepStdAll";
            break;
        }
    }

    return res;
}

SEXP GetDesign(const PartDesign &part, int lenV, bool verbose) {

    std::vector<int> vMap(lenV);
    const int strt = part.mapZeroFirst ? 0 : 1;
    std::iota(vMap.begin(), vMap.end(), strt);
    const std::string ptype = GetPartitionType(part);

    if (verbose) {
        RcppThread::Rcout << "          Partition Design Overview\n";
        RcppThread::Rcout << "*********************************************\n\n";

        if (part.isMult) {
            RcppThread::Rcout << "Partitions of Multiset of width: " <<
                part.width << "\n\n";
        } else if (part.isRep) {
            RcppThread::Rcout << "Partitions with Repetition of width: " <<
                part.width << "\n\n";
        } else {
            RcppThread::Rcout << "Distinct Partitions of width: " <<
                part.width << "\n\n";
        }

        RcppThread::Rcout << "Partition Type:  " << ptype << "\n\n";

        std::string strBool = (part.mIsNull) ? "TRUE" : "FALSE";
        RcppThread::Rcout << "Is m NULL?: " << strBool << "\n";
        strBool = (part.solnExist) ? "TRUE" : "FALSE";
        RcppThread::Rcout << "Does Soln Exist?: " << strBool << "\n";

        RcppThread::Rcout << "\nThe isomorphic vector:\nv: ";

        for (int i = 0; i < (lenV - 1); ++i)
            RcppThread::Rcout << vMap[i] << ", ";

        RcppThread::Rcout << vMap.back() << "\n\n";
        RcppThread::Rcout << "The first indexing vector is given by:\nstartZ: ";

        for (int i = 0; i < (part.startZ.size() - 1); ++i)
            RcppThread::Rcout << part.startZ[i] << ", ";

        RcppThread::Rcout << part.startZ.back() << "\n\n";
        RcppThread::Rcout << "Number of partitions: " <<  part.count << "\n\n";

        RcppThread::Rcout << "Shift:           " <<  part.shift << "\n";
        RcppThread::Rcout << "Slope:           " <<  part.slope << "\n\n";
        RcppThread::Rcout << "Mapped target:   " <<  part.mapTar << "\n";
        RcppThread::Rcout << "Original target: " <<  part.target << "\n\n";

        RcppThread::Rcout << "Confirm MappedTar = (Target + Width * Shift) / Slope\n";
        const std::string eqn_check_str = std::to_string(part.mapTar) + " == (" +
            std::to_string(part.target) + " + " + std::to_string(part.width) +
            " * " + std::to_string(part.shift) + ") / " + std::to_string(part.slope);
        RcppThread::Rcout << eqn_check_str << std::endl;
    }

    bool eqn_check_val = part.mapTar == (part.target +
                                         part.width * part.shift) /
                                         part.slope;

    SEXP sexp_vec = PROTECT(Rf_allocVector(INTSXP, lenV));
    SEXP sexp_index = PROTECT(Rf_allocVector(INTSXP, part.startZ.size()));

    for (int i = 0; i < lenV; ++i)
        INTEGER(sexp_vec)[i] = vMap[i];

    for (int i = 0; i < part.startZ.size(); ++i)
        INTEGER(sexp_index)[i] = part.startZ[i] + 1;

    const char *names[] = {"num_partitions", "mapped_vector",
                           "mapped_target", "first_index_vector",
                           "eqn_check", "partition_type", ""};

    SEXP res = PROTECT(Rf_mkNamed(VECSXP, names));

    SET_VECTOR_ELT(res, 0, Rf_ScalarReal(part.count));
    SET_VECTOR_ELT(res, 1, sexp_vec);
    SET_VECTOR_ELT(res, 2, Rf_ScalarInteger(part.mapTar));
    SET_VECTOR_ELT(res, 3, sexp_index);
    SET_VECTOR_ELT(res, 4, Rf_ScalarLogical(eqn_check_val));
    SET_VECTOR_ELT(res, 5, Rf_mkString(ptype.c_str()));

    UNPROTECT(3);
    return res;
}
