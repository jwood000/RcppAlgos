#include "cpp11/R.hpp"

const char* RcppAlgos_version = "2.10.1";

/**
 * This file records the expected package version in the shared
 * library (or DLL) of the package. This is useful to check that users
 * have properly installed your package. Installation issues where the
 * package is updated but the DLL isn't are common on Windows in
 * particular. To automatically check that the native library of the
 * package was properly installed:
 *
 * - Register the function below as a C callable under the name
 *   "RcppAlgos_linked_version".
 *
 * - Call `rlang::check_linked_version(pkg_name)` from your
 *   `.onLoad()` hook. If you don't depend on rlang copy the
 *   standalone-linked-version.R file from the rlang repository to your R
 *   folder. Find it at
 *   <https://github.com/r-lib/rlang/blob/main/R/standalone-linked-version.R>
 */

[[cpp11::register]]
SEXP linked_version(void) {
    return Rf_mkString(RcppAlgos_version);
}
