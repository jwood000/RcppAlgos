physicalCoreCount <- function() {
    
    if (.Platform$OS.type == "windows")
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
            if (is.null(a <- tryCatch(suppressWarnings(system(cmd, TRUE)), error = function(e) NULL))) 
                next
            a <- gsub("^ +", "", a[1])
            if (grepl("^[1-9]", a)) 
                return(as.integer(a))
        }
    
    ## if the above fails, we try all of them
    for (i in seq(systems))
        for (cmd in systems[i]) {
            if (is.null(a <- tryCatch(suppressWarnings(system(cmd, TRUE)), error = function(e) NULL))) 
                next
            a <- gsub("^ +", "", a[1])
            if (grepl("^[1-9]", a)) 
                return(as.integer(a))
        }
    
    NA_integer_
}

getL2Cache <- function(nCores, nThreads) {
    
    # Typically, there is 256KiB per core of L2 Cache. This is an estimate and
    # more care should be taken when determining this for critical code.
    #
    # For a proper handling of this subject see here:
    #         https://github.com/kimwalisch/primesieve/blob/master/src/CpuInfo.cpp
    #
    # Our simple treament here assumes hyperthreading with 2 threads per core, as
    # this is common on the market. This means we must divide 256KiB by 2 to obtain
    # the amount of memory per thread if physicalCoreCount fails.
    #
    # For all of these reasons, we provide an argument for total L2 Cache Size
    # (in Kibibytes) for a more fine-tuned user approach.
    
    if (is.na(nCores)) {
        if (is.na(nThreads)) {
            L2Cache = 256
        } else {
            L2Cache = nThreads * 128
        }
    } else {
        L2Cache = nCores * 256
    }
    
    L2Cache
}
