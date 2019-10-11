purifyVector <- function(v, freqs) {
    
    vIsNA <- is.na(v)
    
    if (!is.null(freqs) && length(v) != length(freqs)) {
        stop("the length of freqs must equal the length of v")
    }
    
    if (any(vIsNA)) {
        if (!is.null(freqs)) freqs <- freqs[-vIsNA]
        v <- v[-vIsNA]
    }
    
    vOut <- sort(unique(v))
    
    if (length(vOut) < length(v)) {
        vReps <- tabulate(match(v, vOut))
        
        if (is.null(freqs)) {
            return(list(vOut, vReps))
        } else {
            freqs <- vapply(split(freqs, v), sum, 1, USE.NAMES = F)
            return(list(vOut, freqs))
        }
    }
    
    return(list(vOut, NULL))
}

comboGeneral <- function(v, m = NULL, repetition = FALSE, freqs = NULL, 
                         lower = NULL, upper = NULL, constraintFun = NULL,
                         comparisonFun = NULL, limitConstraints = NULL, 
                         keepResults = NULL, FUN = NULL, Parallel = FALSE,
                         nThreads = NULL, tolerance = NULL, cleanVec = FALSE) {
    
    if (cleanVec) {
        myL <- purifyVector(v, freqs)
        v <- myL[[1]]; freqs <- myL[[2]]
    }
    
    isFactor <- is.factor(v)
    CombinatoricsRcpp(v, m, repetition, freqs, lower, upper, constraintFun, 
                      comparisonFun, limitConstraints, TRUE, keepResults, isFactor, FALSE,
                      FUN, new.env(), Parallel, nThreads, pkgEnv$nThreads, tolerance, cleanVec)
}

permuteGeneral <- function(v, m = NULL, repetition = FALSE, freqs = NULL, 
                           lower = NULL, upper = NULL, constraintFun = NULL,
                           comparisonFun = NULL, limitConstraints = NULL,
                           keepResults = NULL, FUN = NULL, Parallel = FALSE,
                           nThreads = NULL, tolerance = NULL, cleanVec = FALSE) {
    
    if (cleanVec) {
        myL <- purifyVector(v, freqs)
        v <- myL[[1]]; freqs <- myL[[2]]
    }
    
    isFactor <- is.factor(v)
    CombinatoricsRcpp(v, m, repetition, freqs, lower, upper, constraintFun, 
                      comparisonFun, limitConstraints, FALSE, keepResults, isFactor, FALSE,
                      FUN, new.env(), Parallel, nThreads, pkgEnv$nThreads, tolerance, cleanVec)
}

comboCount <-  function(v, m = NULL, repetition = FALSE, freqs = NULL, cleanVec = FALSE) {
    
    if (cleanVec) {
        myL <- purifyVector(v, freqs)
        v <- myL[[1]]; freqs <- myL[[2]]
    }
    
    CombinatoricsRcpp(v, m, repetition, freqs, NULL, NULL,
                      NULL, NULL, NULL, TRUE, FALSE, FALSE,
                      TRUE, NULL, NULL, FALSE, NULL, 0, NULL)
}

permuteCount <- function(v, m = NULL, repetition = FALSE, freqs = NULL, cleanVec = FALSE) {
    
    if (cleanVec) {
        myL <- purifyVector(v, freqs)
        v <- myL[[1]]; freqs <- myL[[2]]
    }
    
    CombinatoricsRcpp(v, m, repetition, freqs, NULL, NULL,
                      NULL, NULL, NULL, FALSE, FALSE, FALSE,
                      TRUE, NULL, NULL, FALSE, NULL, 0, NULL, cleanVec)
}

comboSample <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                        n = NULL, sampleVec = NULL, seed = NULL, 
                        FUN = NULL, Parallel = FALSE, nThreads = NULL, 
                        namedSample = FALSE, cleanVec = FALSE) {
    
    if (cleanVec) {
        myL <- purifyVector(v, freqs)
        v <- myL[[1]]; freqs <- myL[[2]]
    }
    
    isFactor <- is.factor(v)
    if (!is.null(seed)) {set.seed(seed)}
    SampleRcpp(v, m, repetition, freqs, sampleVec, TRUE, isFactor,
               seed, n, sample, FUN, new.env(), Parallel, nThreads,
               pkgEnv$nThreads, namedSample)
}

permuteSample <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                          n = NULL, sampleVec = NULL, seed = NULL, 
                          FUN = NULL, Parallel = FALSE, nThreads = NULL, 
                          namedSample = FALSE, cleanVec = FALSE) {
    
    if (cleanVec) {
        myL <- purifyVector(v, freqs)
        v <- myL[[1]]; freqs <- myL[[2]]
    }
    
    isFactor <- is.factor(v)
    if (!is.null(seed)) {set.seed(seed)}
    SampleRcpp(v, m, repetition, freqs, sampleVec, FALSE, isFactor,
               seed, n, sample, FUN, new.env(), Parallel, nThreads,
               pkgEnv$nThreads, namedSample)
}

comboGroups <- function(v, numGroups, retType = "matrix", lower = NULL,
                       upper = NULL, Parallel = FALSE, nThreads = NULL) {
    isFactor <- is.factor(v)
    ComboGroupsRcpp(v, numGroups, retType, lower, upper, 
                    isFactor, FALSE, Parallel, nThreads, pkgEnv$nThreads,
                    FALSE, NULL, NULL, NULL, sample, FALSE)
}

comboGroupsSample <- function(v, numGroups, retType = "matrix", n = NULL,
                              sampleVec = NULL, seed = NULL, Parallel = FALSE, 
                              nThreads = NULL, namedSample = FALSE) {
    isFactor <- is.factor(v)
    if (!is.null(seed)) {set.seed(seed)}
    ComboGroupsRcpp(v, numGroups, retType, NULL, NULL, 
                    isFactor, FALSE, Parallel, nThreads, pkgEnv$nThreads,
                    TRUE, sampleVec, seed, n, sample, namedSample)
}

comboGroupsCount <- function(v, numGroups) {
    ComboGroupsRcpp(v, numGroups, NULL, NULL, NULL, FALSE, TRUE, 
                    FALSE, FALSE, 0, FALSE, NULL, NULL, NULL, sample, FALSE)
}

