## One fresh process per method/input/thread/block. Timing includes the public
## R wrapper and the complete count, but excludes package/process setup.
args <- commandArgs(TRUE)
stopifnot(length(args) == 6L)
libraryPath <- normalizePath(args[[1L]], mustWork = TRUE)
n <- as.numeric(args[[2L]])
threads <- as.integer(args[[3L]])
methodName <- args[[4L]]
expected <- as.numeric(args[[5L]])
timerDll <- dyn.load(normalizePath(args[[6L]], mustWork = TRUE))
timerSymbol <- getNativeSymbolInfo("AutoDispatchMilliseconds", PACKAGE = timerDll)$address
library(RcppAlgos, lib.loc = libraryPath)
dll <- normalizePath(getLoadedDLLs()[["RcppAlgos"]][["path"]])
stopifnot(startsWith(dll, paste0(libraryPath, "/")),
          threads >= 1L, threads <= RcppAlgos:::pkgEnv$maxThreads,
          methodName %in% c("legendre", "meissel", "lehmer"))

## Avoid introducing JIT compilation during a measured sequence. The installed
## package is already byte-compiled. This warms method matching/conversion and
## the wrapper via its small sieve path without evaluating or warming Phi.
invisible(compiler::enableJIT(0L))
for (i in seq_len(64L)) invisible(primeCount(9973L, method = "auto"))
tick <- compiler::cmpfun(function() .Call(timerSymbol))
for (i in seq_len(16L)) invisible(tick())
count <- compiler::cmpfun(function() {
    primeCount(n, nThreads = threads, method = methodName)
})
times <- numeric(4L)
results <- numeric(4L)
invisible(gc())
for (i in seq_len(4L)) {
    start <- tick()
    answer <- count()
    times[[i]] <- tick() - start
    results[[i]] <- answer
    stopifnot(identical(answer, if (expected <= .Machine$integer.max)
                                  as.integer(expected) else expected))
}
write.table(data.frame(method = methodName, n = n, threads = threads,
                       repeat_index = 0:3, result = results, total_ms = times,
                       maxThreads = RcppAlgos:::pkgEnv$maxThreads),
            file = stdout(), sep = ",", row.names = FALSE, quote = FALSE)
