#include "Constraints/ConstraintsTypes.h"
#include "Partitions/PartitionsUtils.h"
#include "CppConvert.h"
#include <numeric>

static std::string GetCTypeName(ConstraintType ctype) {
    switch (ctype) {
        case ConstraintType::General:
            return "General";
        case ConstraintType::NoConstraint:
            return "NoConstraint";
        case ConstraintType::PartitionEsque:
            return "PartitionEsque";
        case ConstraintType::PartMapping:
            return "PartMapping";
        case ConstraintType::PartStandard:
            return "PartStandard";
        case ConstraintType::SpecialCnstrnt:
            return "SpecialCnstrnt";
        default:
            return "Constraint Not Recognized";
    }
}

SEXP GetDesign(const PartDesign &part, ConstraintType ctype,
               int lenV, bool verbose) {

    std::vector<int> vMap(lenV);
    const int strt = part.includeZero ? 0 : 1;

    std::iota(vMap.begin(), vMap.end(), strt);
    const std::string str_ptype = GetPTypeName(part.ptype);
    const std::string str_ctype = GetCTypeName(ctype);

    cpp11::sexp sexp_count = part.numUnknown
        ? Rf_ScalarReal(NA_REAL)
        : CppConvert::GetCount(part.isGmp, part.bigCount, part.count);

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

        Rprintf("Partition Type: %s\n\n", str_ptype.c_str());

        std::string strBool = (part.mIsNull) ? "TRUE" : "FALSE";
        Rprintf("Is m NULL?: %s\n", strBool.c_str());
        strBool = part.solnExist ? "TRUE" : "FALSE";
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

        if (part.numUnknown) {
            Rprintf("Number of partitions: Unknown\n");
        } else if (part.isGmp) {
            const std::string bigCount = part.bigCount.get_str();
            Rprintf("Number of partitions: %s\n", bigCount.c_str());
        } else {
            Rprintf("Number of partitions: %s\n",
                    std::to_string(part.count).c_str());
        }

        if (part.solnExist) {
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
    }

    const bool canCheck = part.solnExist && part.slope != 0;
    bool eqn_check_val = false;

    if (canCheck) {
        eqn_check_val = part.mapTar ==
            (part.target + part.width * part.shift) / part.slope;
    }

    cpp11::sexp sexp_eqn_check = Rf_ScalarLogical(
        canCheck ? static_cast<int>(eqn_check_val) : NA_LOGICAL
    );

    cpp11::sexp sexp_vec = Rf_allocVector(INTSXP, lenV);
    cpp11::sexp sexp_index = Rf_allocVector(INTSXP, part.startZ.size());

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
                           "eqn_check", "partition_type", "count_unknown",
                           "is_partition", "constraint_type", ""};

    cpp11::sexp res = Rf_mkNamed(VECSXP, names);
    SET_VECTOR_ELT(res, 0, sexp_count);
    SET_VECTOR_ELT(res, 1, sexp_vec);
    SET_VECTOR_ELT(res, 2, Rf_ScalarInteger(part.mapTar));
    SET_VECTOR_ELT(res, 3, sexp_index);
    SET_VECTOR_ELT(res, 4, sexp_eqn_check);
    SET_VECTOR_ELT(res, 5, Rf_mkString(str_ptype.c_str()));
    SET_VECTOR_ELT(res, 6, Rf_ScalarLogical(part.numUnknown));
    SET_VECTOR_ELT(res, 7, Rf_ScalarLogical(part.isPart));
    SET_VECTOR_ELT(res, 8, Rf_mkString(str_ctype.c_str()));

    return res;
}
