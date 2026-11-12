args <- commandArgs(TRUE)
artifact <- normalizePath(args[[1L]], mustWork = TRUE)
libraryPath <- file.path(artifact, "library")
library(RcppAlgos, lib.loc = libraryPath)
options(warn = 2)
dll <- getLoadedDLLs()[["RcppAlgos"]]
stopifnot(startsWith(normalizePath(dll[["path"]]),
                    paste0(normalizePath(libraryPath), "/")))
stopifnot(identical(formals(primeCount),
                    alist(n = , nThreads = NULL, method = "auto")))
routine <- getDLLRegisteredRoutines(dll)[[".Call"]][["_RcppAlgos_PrimeCountCpp"]]
stopifnot(identical(routine$numParameters, 4L))
rd <- file.path(artifact, "package-source/man/primeCount.Rd")
invisible(tools::parse_Rd(rd))
tools::checkRd(rd)
example("primeCount", package = "RcppAlgos", lib.loc = libraryPath,
        ask = FALSE, echo = FALSE)
cat("PASS primeCount Rd/examples, unchanged formals and four-argument registration\n")
