# Run from the repository root:
# Rscript --vanilla .github/development/native-regressions/run.R [new_output_directory]
main <- function() {
    repository <- normalizePath(".", winslash = "/")
    fixtures <- file.path(repository, ".github/development/native-regressions")
    stopifnot(file.exists(file.path(fixtures, "check.h")),
              requireNamespace("cpp11", quietly = TRUE),
              requireNamespace("gmp", quietly = TRUE))

    arguments <- commandArgs(TRUE)
    if (length(arguments) > 1L) stop("Expected at most one new output directory")
    output <- if (length(arguments)) arguments[[1L]] else
        tempfile("rcppalgos-native-regressions-")
    if (file.exists(output)) stop("Use a new output directory: ", output)
    if (!dir.create(output, recursive = TRUE)) stop("Cannot create output directory")
    output <- normalizePath(output, winslash = "/")
    cat("Native regression output: ", output, "\n", sep = "")

    snapshot <- file.path(output, "source")
    dir.create(snapshot)
    dir.create(file.path(snapshot, "src"))
    dir.create(file.path(snapshot, "checks"))

    # Link the actual production implementations and their dependencies.
    # No installed RcppAlgos DLL, copied algorithms, or counting stubs are used.
    source_names <- c(
        "BigComboCount.cpp", "BigPartsCountDistinct.cpp",
        "BigPartsCountMultiset.cpp", "BigPartsCountRep.cpp",
        "BigPartsCountSection.cpp", "BigPermuteCount.cpp", "ComboCount.cpp",
        "ConvertUtils.cpp", "GmpConvert.cpp", "MultisetCountEngine.cpp",
        "NthCombination.cpp", "NthComposition.cpp", "NthPartition.cpp",
        "NthPartsDispatch.cpp", "NthPermutation.cpp", "PartitionsCount.cpp",
        "PartitionsCountDistinct.cpp", "PartitionsCountMultiset.cpp",
        "PartitionsCountRep.cpp", "PartitionsCountSection.cpp",
        "PartitionsUtils.cpp", "PermuteCount.cpp", "PollardRhoUtils.cpp",
        "RankComposition.cpp", "RankPartition.cpp", "RankPartsDispatch.cpp",
        "SetUpUtils.cpp", "StdConvert.cpp"
    )
    fixture_names <- c("gmp-convert.cpp", "primality.cpp",
                       "multiset-gmp-buffers.cpp", "multiset-gmp-internals.cpp")
    stopifnot(all(file.copy(file.path(repository, "src", source_names),
                           file.path(snapshot, "src"))),
              file.copy(file.path(repository, "inst/include"), snapshot,
                        recursive = TRUE),
              all(file.copy(file.path(fixtures, c(fixture_names, "check.h", "run.R")),
                            file.path(snapshot, "checks"))))

    paths <- list.files(snapshot, recursive = TRUE, full.names = TRUE)
    write.csv(data.frame(file = substring(paths, nchar(snapshot) + 2L),
                         md5 = unname(tools::md5sum(paths))),
              file.path(output, "source-hashes.csv"), row.names = FALSE)
    writeLines(capture.output(sessionInfo()), file.path(output, "session-info.txt"))

    previous <- Sys.getenv(c("CPATH", "PKG_LIBS"), unset = NA_character_)
    on.exit(for (name in names(previous)) {
        if (is.na(previous[[name]])) Sys.unsetenv(name)
        else do.call(Sys.setenv, setNames(list(previous[[name]]), name))
    }, add = TRUE)
    includes <- c(file.path(snapshot, "include"), file.path(snapshot, "checks"))
    if (!is.na(previous[["CPATH"]]) && nzchar(previous[["CPATH"]])) {
        includes <- c(includes, previous[["CPATH"]])
    }
    libraries <- if (is.na(previous[["PKG_LIBS"]])) "" else previous[["PKG_LIBS"]]
    Sys.setenv(CPATH = paste(includes, collapse = .Platform$path.sep),
               PKG_LIBS = paste(libraries, "-lgmp"))

    compiled <- new.env(parent = globalenv())
    build <- file.path(output, "build")
    # Existing production counting factories use std::make_unique (C++14).
    # This affects only the temporary fixture build, not package requirements.
    cpp11::cpp_source(
        c(file.path(snapshot, "checks", fixture_names),
          file.path(snapshot, "src", source_names)),
        env = compiled, clean = FALSE, quiet = FALSE,
        cxx_std = "CXX14", dir = build
    )

    checks <- c("check_signed_gmp_conversions", "check_multiset_buffer_reuse",
                "check_reset_p2d", "check_multiset_changing_reps",
                "check_multiset_gmp_rank_unrank", "check_primality_upper_bound",
                "check_primality_out_of_range")
    timings <- vapply(checks, function(check) {
        elapsed <- system.time(compiled[[check]]())[["elapsed"]]
        cat(sprintf("PASS %s (%.3fs)\n", check, elapsed))
        elapsed
    }, numeric(1))
    write.csv(data.frame(check = checks, elapsed_seconds = unname(timings)),
              file.path(output, "results.csv"), row.names = FALSE)
    cat(sprintf("PASS all %d native regression checks; output=%s\n",
                length(checks), output))
}
main()
