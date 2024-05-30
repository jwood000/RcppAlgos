rperm <- function(m, size=2) { # Obtain m unique permutations of 1:size
    max.failures <- 10

    # Function to index into the upper-level cache.
    prefix <- function(p, k) {    # p is a permutation, k is the prefix size
        sum((p[1:k] - 1) * (size ^ ((1:k)-1))) + 1
    } # Returns a value from 1 through size^k

    # Function to obtain a new permutation.
    newperm <- function() {
        # References cache, k.head, and failures in parent context.
        # Modifies cache and failures.

        count <- 0                # Protects against infinite loops
        repeat {
            # Generate a permutation and check against previous ones.
            p <- sample(1:size)
            k <- prefix(p, k.head)
            ip <- cache[[k]]
            hash.p <- paste(tail(p,-k.head), collapse="")
            if (is.null(ip[[hash.p]])) break

            # Prepare to try again.
            n.failures <<- n.failures + 1
            count <- count+1
            if (count > max.failures) {
                p <- NA           # NA indicates a new permutation wasn't found
                hash.p <- ""
                break
            }
        }
        if (count <= max.failures) {
            ip[[hash.p]] <- TRUE      # Update the list of permutations found
            cache[[k]] <<- ip
        }
        p                         # Return this (new) permutation
    }

    # Initialize the cache.
    k.head <- min(size-1, max(1, floor(log(m / log(m)) / log(size))))
    cache <- as.list(1:(size^k.head))
    for (i in 1:(size^k.head)) cache[[i]] <- list()

    # Count failures (for benchmarking and error checking).
    n.failures <- 0

    # Obtain (up to) m unique permutations.
    s <- replicate(m, newperm())
    s[is.na(s)] <- NULL
    list(failures=n.failures, sample=matrix(unlist(s), ncol=size))
} # Returns an m by size matrix; each row is a permutation of 1:size.