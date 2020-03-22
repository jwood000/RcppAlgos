\name{RcppAlgos-package}
\alias{RcppAlgos-package}
\alias{RcppAlgos}
\docType{package}

\title{
High Performance Tools for Combinatorics and Computational Mathematics
}

\description{
The \strong{RcppAlgos} package attacks age-old problems in combinatorics and computational mathematics.
}

\section{Goals}{
    \enumerate{
        \item{The main goal is to encourage fresh and creative approaches to foundational problems. The question that most appropriately summarizes \code{RcppAlgos} is: \strong{\emph{"Can we do better?"}}.}
        \item{Provide highly optimized functions that facilitates a broader spectrum of researchable cases. \emph{E.g}
            \itemize{
                \item{Investigating the structure of large numbers over wide ranges:
                    \itemize{
                        \item{\code{primeFactorizeSieve(10^13 - 10^7, 10^13 + 10^7)}}
                        \item{\code{primeSieve(2^53 - 10^10, 2^53 - 1, nThreads = 32)}}
                        }
                    }
                \item{Easily explore combinations/permutations that would otherwise be inaccessible due to time of execution/memory constraints:
                    \itemize{
                        \item{\preformatted{iter = comboIter(10000, 100)
bigSamp = gmp::urand.bigz(3, gmp::log2.bigz(comboCount(10000, 100)))
iter[[bigSamp]]}}
                        \item{\code{comboGeneral(150, 5, constraintFun = "sum", Parallel = TRUE)}}
                        \item{\preformatted{parallel::mclapply(seq(...), function(x) {
            temp = permuteGeneral(15, 10, lower = x, upper = y)
            ## analyze permutations
            ## output results
    }, mc.cores = detectCores() - 1))}}
                        \item{\preformatted{comboGeneral(0:80, m = 80, repetition = TRUE,
            constraintFun = "sum", comparisonFun = "==",
            limitConstraints = 80)}}
                        \item{\code{permuteSample(rnorm(100), 10, freqs = rep(1:4, 25), n = 15, seed = 123)}}
                        \item{\preformatted{set.seed(123)
comboGeneral(runif(42, 0, 50), 10,
             constraintFun = "mean",
             comparisonFun = c(">","<"),
             limitConstraints = c(39.876, 42.123))}}
                    }
                }
            }}
        \item{\emph{Speed!!!...}. You will find that the functions in \code{RcppAlgos} are some of the fastest of their type available in \code{R}.}
    }
}

\author{
Joseph Wood
}

\keyword{ package }
