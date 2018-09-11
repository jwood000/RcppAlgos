\name{comboSample}
\alias{comboSample}
\alias{permuteSample}
\title{
Sample combinations/permutations
}
\description{
    \itemize{
        \item Generate a specific (lexicographically) or random sample of combinations/permutations.
        \item Produce results in parallel using the \code{Parallel} argument.
        \item GMP support allows for exploration of combinations/permutations of vectors with many elements.
    }
}
\usage{
comboSample(v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL,
             sampleVec = NULL, seed = NULL, FUN = NULL, Parallel = FALSE)

permuteSample(v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL,
               sampleVec = NULL, seed = NULL, FUN = NULL, Parallel = FALSE)
}
\arguments{
\item{v}{Source vector. If \code{v} is an integer, it will be converted to the sequence \code{1:v}.}
  \item{m}{Number of elements to choose. If \code{repetition = TRUE}, \code{m} can exceed the length of \code{v}. The default is \code{NULL}.}
  \item{repetition}{Logical value indicating whether combinations/permutations should be with or without repetition. The default is \code{FALSE}.}
  \item{freqs}{A vector of frequencies used for producing all combinations/permutations of a multiset of \code{v}. Each element of \code{freqs} represents how many times each element of the source vector, \code{v}, is repeated. It is analogous to the \code{times} argument in \code{\link[base]{rep}}. The default value is \code{NULL}.}
  \item{n}{Number of combinations/permutations to return. The default is \code{NULL}.}
  \item{sampleVec}{A vector of numbers representing the lexicographical combination/permutations to return. Accepts vectors of class \code{bigz} as well as vectors of characters}
  \item{seed}{Random seed initialization. The default is \code{NULL}.}
  \item{FUN}{Function to be applied to each combination/permutation. The default is \code{NULL}.}
  \item{Parallel}{Logical value indicating whether combinations/permutations should be generated in parallel. The default is \code{FALSE}. If \code{TRUE}, the number of cores used is equal to the minimum of the number of cores available on your system minus one and the number of results requested (\emph{E.g.} if user has 16 cores and only needs 5 results, 5 cores will be used). If your system has fewer than 3 cores, \code{Parallel} is forced to be \code{FALSE}.}
}
\details{
These algorithms rely on efficiently generating the \eqn{n^{th}}{nth} lexicographical combination/permutation (sometimes called the \href{https://rosettacode.org/wiki/Permutations/Rank_of_a_permutation}{rank}).
}
\value{
In general, a matrix is returned with each row containing a vector of length \eqn{m}. If \code{m} isn't supplied and \code{freqs} is given, a matrix is returned with each row containing a vector of length \code{sum(freqs)}. If \code{FUN} is utilized, a list is returned.
}
\references{
\href{https://en.wikipedia.org/wiki/Lexicographical_order}{Lexicographical order}
}
\author{
Joseph Wood
}
\note{
    \itemize{
        \item \code{n} and \code{sampleVec} cannot both be \code{NULL}.
        \item Factor vectors are accepted. Class and level attributes are preserved.
        \item {\code{Parallel} will be ignored in the following cases (Research is on-going to enable parallelism in all cases below):
            \itemize{
                \item If the class of the vector passed is \code{character} (N.B. \code{Rcpp::CharacterMatrix} is not thread safe). Alternatively, you can generate an indexing matrix in parallel.
                \item If \code{FUN} is utilized.
            }
        }
    }
}

\examples{
set.seed(11)
comboSample(30, 15, TRUE, n = 10)

set.seed(2)
permuteSample(15, 10, freqs = sample(1:2, 15, TRUE), 
              sampleVec = c(1, 10^2, 10^5, 10^8, 10^11))
              
all.equal(comboSample(10, 5, 
            sampleVec = 1:comboCount(10, 5)),
         comboGeneral(10, 5))
         
## Examples with enormous number of total permutations
num = permuteCount(10000, 20)
gmp::log2.bigz(num)
## [1] 265.7268

first = gmp::urand.bigz(n = 1, size = 265, seed = 123)
mySamp = do.call(c, lapply(0:10, function(x) gmp::add.bigz(first, x)))

class(mySamp)
## [1] "bigz"

## using permuteSample
pSamp = permuteSample(10000, 20, sampleVec = mySamp)

## using permuteGeneral
pGeneral = permuteGeneral(10000, 20, 
                          lower = first,
                          upper = gmp::add.bigz(first, 10))

identical(pSamp, pGeneral)
## [1] TRUE

## Using FUN and Parallel
permuteSample(10000, 50, n = 3, seed = 10, FUN = sd, Parallel = TRUE)
# Seed initialisation
# [[1]]
# [1] 2694.073
# 
# [[2]]
# [1] 2923.153
# 
# [[3]]
# [1] 2971.491

}