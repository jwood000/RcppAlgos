#include "Constraints/ConstraintsTypes.h"
#include "Partitions/PartitionsTypes.h"
#include "CleanConvert.h"
#include <numeric>

std::string GetPartitionType(const PartDesign &part) {

    std::string res;

    switch (part.ptype) {
        case PartitionType::NotPartition: {
            res = "NotPartition";
            break;
        } case PartitionType::LengthOne: {
            res = "LengthOne";
            break;
        } case PartitionType::DstctCapped: {
            res = "DistCapped";
            break;
        } case PartitionType::DstctCappedMZ: {
            res = "DstctCappedMZ";
            break;
        } case PartitionType::DstctNoZero : {
            res = "DstctNoZero";
            break;
        } case PartitionType::DstctOneZero: {
            res = "DstctOneZero";
            break;
        } case PartitionType::DstctMultiZero : {
            res = "DstctMultiZero";
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

SEXP GetDesign(const PartDesign &part, ConstraintType ctype,
               int lenV, bool verbose) {

    std::vector<int> vMap(lenV);
    const int strt = part.includeZero ? 0 : 1;
    std::iota(vMap.begin(), vMap.end(), strt);
    const std::string ptype = GetPartitionType(part);

    if (verbose) {
        Rprintf("          Partition Design Overview\n");
        Rprintf("*********************************************\n\n");
        std::string strWidth = std::to_string(part.width);

        if (part.isMult) {
            Rprintf("Partitions of Multiset of width: %s\n",
                    strWidth.c_str());
        } else if (part.isRep) {
            Rprintf("Partitions with Repetition of width: %s\n",
                    strWidth.c_str());
        } else {
            Rprintf("Distinct Partitions of width: %s\n",
                    strWidth.c_str());
        }

        Rprintf("Partition Type: %s\n\n", ptype.c_str());

        std::string strBool = (part.mIsNull) ? "TRUE" : "FALSE";
        Rprintf("Is m NULL?: %s\n", strBool.c_str());
        strBool = (part.solnExist) ? "TRUE" : "FALSE";
        Rprintf("Does Soln Exist?: %s\n", strBool.c_str());

        Rprintf("\nThe isomorphic vector:\nv: ");
        std::string strVMap = "";

        for (auto v_i: vMap) {
            strVMap += (std::to_string(v_i) + " ");
        }

        Rprintf("%s\n\n", strVMap.c_str());
        Rprintf("The first indexing vector is given by:\nstartZ: ");
        std::string strStartZ = "";

        for (auto z_i: part.startZ) {
            strStartZ += (std::to_string(z_i) + " ");
        }

        Rprintf("%s\n\n", strStartZ.c_str());
        Rprintf("Number of partitions: %s\n",
                std::to_string(part.count).c_str());
        Rprintf("Shift:           %s\n",
                std::to_string(part.shift).c_str());
        Rprintf("Slope:           %s\n",
                std::to_string(part.slope).c_str());
        Rprintf("Mapped target:   %s\n",
                std::to_string(part.mapTar).c_str());
        Rprintf("Original target: %s\n\n",
                std::to_string(part.target).c_str());

        Rprintf("Confirm MappedTar = (Target + Width * Shift) / Slope\n");
        const std::string eqn_check_str = std::to_string(part.mapTar) + " == (" +
            std::to_string(part.target) + " + " + std::to_string(part.width) +
            " * " + std::to_string(part.shift) + ") / " + std::to_string(part.slope);
        Rprintf("%s\n\n", eqn_check_str.c_str());
    }

    bool eqn_check_val = part.mapTar == (part.target +
                                         part.width * part.shift) /
                                         part.slope;

    SEXP sexp_vec = PROTECT(Rf_allocVector(INTSXP, lenV));
    SEXP sexp_index = PROTECT(Rf_allocVector(INTSXP, part.startZ.size()));

    for (int i = 0; i < lenV; ++i) {
        INTEGER(sexp_vec)[i] = vMap[i];
    }

    if (ctype == ConstraintType::PartStandard) {
        for (std::size_t i = 0; i < part.startZ.size(); ++i) {
            INTEGER(sexp_index)[i] = part.startZ[i] +
                static_cast<int>(part.includeZero);
        }
    } else {
        for (std::size_t i = 0; i < part.startZ.size(); ++i) {
            INTEGER(sexp_index)[i] = part.startZ[i] + 1;
        }
    }

    const char *names[] = {"num_partitions", "mapped_vector",
                           "mapped_target", "first_index_vector",
                           "eqn_check", "partition_type", ""};

    SEXP res = PROTECT(Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(res, 0, CleanConvert::GetCount(part.isGmp, part.bigCount, part.count));
    SET_VECTOR_ELT(res, 1, sexp_vec);
    SET_VECTOR_ELT(res, 2, Rf_ScalarInteger(part.mapTar));
    SET_VECTOR_ELT(res, 3, sexp_index);
    SET_VECTOR_ELT(res, 4, Rf_ScalarLogical(eqn_check_val));
    SET_VECTOR_ELT(res, 5, Rf_mkString(ptype.c_str()));

    UNPROTECT(3);
    return res;
}
