## Run against the freshly installed, isolated production-integration build.
## Usage: Rscript --vanilla r-tests.R ARTIFACT_ROOT REPOSITORY_ROOT
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
artifact <- normalizePath(args[1L], mustWork = TRUE)
repository <- normalizePath(args[2L], mustWork = TRUE)
libraryPath <- file.path(artifact, "library")
library(RcppAlgos, lib.loc = libraryPath)

dll <- normalizePath(getLoadedDLLs()[["RcppAlgos"]][["path"]])
stopifnot(startsWith(dll, paste0(normalizePath(libraryPath), "/")))
metadata <- list(R = R.version.string,
                 testthat = as.character(packageVersion("testthat")),
                 package = as.character(packageVersion("RcppAlgos")),
                 dll = dll, maxThreads = RcppAlgos:::pkgEnv$maxThreads,
                 session = capture.output(sessionInfo()))
saveRDS(metadata, file.path(artifact, "r-session.rds"))
dput(metadata, file = file.path(artifact, "r-session.txt"))

## PrimeCount runs first so its existing parallel-cache tests start in a fresh
## R process. All files use this isolated installed DLL, not load_all().
files <- c("testPrimeCount.R", "testPartialSieve.R", "testErrors.R",
           "testPrimeSieve.R", "testPrimeFactorize.R",
           "testPrimeFactorizeSieve.R", "testIsPrimeRcpp.R",
           "testEulerPhiSieve.R", "testDivisorsRcpp.R",
           "testDivisorsSieve.R", "testNumDivisorSieve.R")
summaries <- vector("list", length(files))
details <- vector("list", length(files))
for (i in seq_along(files)) {
    cat("\nRunning", files[i], "\n")
    result <- testthat::test_file(file.path(repository, "tests/testthat", files[i]),
                                  reporter = "summary", stop_on_failure = TRUE)
    detail <- as.data.frame(result)
    details[[i]] <- detail
    summaries[[i]] <- data.frame(file = files[i], tests = nrow(detail),
                                expectations = sum(detail$nb),
                                passed = sum(detail$passed),
                                failed = sum(detail$failed),
                                errors = sum(detail$error),
                                warnings = sum(detail$warning),
                                skipped = sum(detail$skipped))
    write.csv(do.call(rbind, summaries[seq_len(i)]),
              file.path(artifact, "r-tests.csv"), row.names = FALSE)
    saveRDS(details[seq_len(i)], file.path(artifact, "r-test-details.rds"))
}
summary <- do.call(rbind, summaries)
print(summary, row.names = FALSE)
stopifnot(sum(summary$failed) == 0L, sum(summary$errors) == 0L,
          sum(summary$warnings) == 0L)

## Only this public help file changes; interfaces and registration are intact.
rd <- file.path(repository, "man/semiPrimeCount.Rd")
tools::parse_Rd(rd)
tools::checkRd(rd)
example("semiPrimeCount", package = "RcppAlgos", lib.loc = libraryPath,
        ask = FALSE, echo = FALSE)
cat("\nPASS R expectations:", sum(summary$passed),
    "; skips:", sum(summary$skipped), "; semiPrimeCount Rd/examples checked\n")
