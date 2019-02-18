pkgEnv <- new.env(parent = emptyenv())
pkgEnv$nCores <- NULL
pkgEnv$nThreads <- NULL

physicalCoreCount <- function() {
    
    if (grepl("darwin|solaris", R.version$os) || .Platform$OS.type == "windows")
        return(parallel::detectCores(logical = FALSE))
    
    ## According to R News as of # 3.4.4 
    ##    "parallel::detectCores(logical = FALSE) is ignored on Linux systems, since
    ##                     the information is not available with virtualized OSes."
    ##
    ## This implementation tries the old command available from version 3.3.0 until
    ## 3.4.3 (It was removed in 3.4.4), as I have had success with this command on
    ## several Linux flavors. If this fails, we also try the method outlined
    ## by user teambob in an answer to:
    ##       "How to obtain the number of CPUs/cores in Linux from the command line?"
    ## https://stackoverflow.com/a/18051445/4408538
    ##
    ## Otherwise, the implementation below is nearly identical to parallel::detectCores
    ## with the argument logical set to FALSE. You will note that we have excluded 
    ## freevsd, openbsd, & irix as they do not have this option available.
    
    systems <- list(linux = c("cat /proc/cpuinfo | grep 'cpu cores'| uniq | cut -f2 -d:", 
                              "grep '^core id' /proc/cpuinfo |sort -u|wc -l"), 
                    darwin = "/usr/sbin/sysctl -n hw.physicalcpu 2>/dev/null", 
                    solaris = "/bin/kstat -p -m cpu_info | grep :core_id | cut -f2 | uniq | wc -l")
    
    for (i in seq(systems)) if (grepl(paste0("^", names(systems)[i]), R.version$os))
        for (cmd in systems[i]) {
            if (is.null(a <- tryCatch(suppressWarnings(
                system(cmd, TRUE)), error = function(e) NULL))) {next}
            
            a <- gsub("^ +", "", a[1])
            if (grepl("^[1-9]", a)) 
                return(as.integer(a))
        }
    
    ## if the above fails, we try all of them
    for (i in seq(systems))
        for (cmd in systems[i]) {
            if (is.null(a <- tryCatch(suppressWarnings(
                system(cmd, TRUE)), error = function(e) NULL))) {next}
            
            a <- gsub("^ +", "", a[1])
            if (grepl("^[1-9]", a)) 
                return(as.integer(a))
        }
    
    ## If we get here, we assume 1 core
    1L
}

## This will set the maximum number of cores
## and number of threads on a given machine
## when the package is loaded
.onLoad <- function(libname, pkgname) {
    pkgEnv$nCores <- physicalCoreCount()
    tempThreads <- stdThreadMax()
    
    if (is.na(tempThreads)) {
        pkgEnv$nThreads <- 1L
    } else {
        pkgEnv$nThreads <- tempThreads
    }
    
    invisible()
}

