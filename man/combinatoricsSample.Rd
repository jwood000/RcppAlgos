\name{comboSample}
\alias{comboSample}
\alias{permuteSample}
\title{
Sample Combinations and Permutations
}
\description{
    \itemize{
        \item Generate a specific (lexicographically) or random sample of combinations/permutations.
        \item Produce results in parallel using the \code{Parallel} or \code{nThreads} arguments.
        \item GMP support allows for exploration of combinations/permutations of vectors with many elements.
    }
}
\usage{
comboSample(v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL,
            sampleVec = NULL, seed = NULL, FUN = NULL, Parallel = FALSE,
            nThreads = NULL, namedSample = FALSE, FUN.VALUE = NULL)

permuteSample(v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL,
              sampleVec = NULL, seed = NULL, FUN = NULL, Parallel = FALSE,
              nThreads = NULL, namedSample = FALSE, FUN.VALUE = NULL)
}
\arguments{
  \item{v}{Source vector. If \code{v} is a positive integer, it will be converted to the sequence \code{1:v}. If \code{v} is a negative integer, it will be converted to the sequence \code{v:-1}. All atomic types are supported (See \code{\link{is.atomic}}).}
  \item{m}{Number of elements to choose.  If \code{repetition = TRUE} or \code{freqs} is utilized, \code{m} can exceed the length of \code{v}. If \code{m = NULL}, the length will default to \code{length(v)} or \code{sum(freqs)}.}
  \item{repetition}{Logical value indicating whether combinations/permutations should be with or without repetition. The default is \code{FALSE}.}
  \item{freqs}{A vector of frequencies used for producing all combinations/permutations of a multiset of \code{v}. Each element of \code{freqs} represents how many times each element of the source vector, \code{v}, is repeated. It is analogous to the \code{times} argument in \code{\link[base]{rep}}. The default value is \code{NULL}.}
  \item{n}{Number of combinations/permutations to return. The default is \code{NULL}.}
  \item{sampleVec}{A vector of indices representing the lexicographical combination/permutations to return. Accepts whole numbers as well as vectors of class \code{bigz} as well as vectors of characters}
  \item{seed}{Random seed initialization. The default is \code{NULL}. N.B. If the gmp library is needed, this parameter must be set in order to have reproducible results (\emph{E.g} \code{set.seed()} has no effect in these cases).}
  \item{FUN}{Function to be applied to each combination/permutation. The default is \code{NULL}.}
  \item{Parallel}{Logical value indicating whether combinations/permutations should be generated in parallel. The default is \code{FALSE}. If \code{TRUE} and \code{nThreads = NULL}, the number of threads used is equal to the minimum of one minus the number of threads available on your system and the number of results requested (\emph{e.g.} if user has 16 threads and only needs 5 results, 5 threads will be used (\emph{i.e.} \code{min(16 - 1, 5) = 5})). If \code{nThreads} is not \code{NULL}, it will be given preference (\emph{e.g.} if user has 8 threads with \code{Parallel = TRUE} and \code{nThreads = 4}, only 4 threads will be spawned). If your system is single-threaded, the arguments \code{Parallel} and \code{nThreads} are ignored.}
  \item{nThreads}{Specific number of threads to be used. The default is \code{NULL}. See \code{Parallel}.}
  \item{namedSample}{Logical flag. If \code{TRUE}, \code{rownames} corresponding to the lexicographical combination/permutation, will be added to the returned matrix. The default is \code{FALSE}.}
  \item{FUN.VALUE}{A template for the return value from \code{FUN}. See 'Details' of \code{\link[base]{vapply}} for more information.}
 }
\details{
These algorithms rely on efficiently generating the \eqn{n^{th}}{nth} lexicographical combination/permutation. This is the process of \href{https://rosettacode.org/wiki/Permutations/Rank_of_a_permutation}{unranking}.
}
\value{
    \itemize{
        \item In general, a matrix with \eqn{m} or \eqn{m + 1} columns, depending on the value of \code{keepResults}
        \item If \code{FUN} is utilized and \code{FUN.VALUE = NULL}, a list is returned
        \item When both \code{FUN} and \code{FUN.VALUE} are not \code{NULL}, the return is modeled after the return of \code{vapply}. See the 'Value' section of \code{\link[base]{vapply}}.
    }
}
\references{
\href{https://en.wikipedia.org/wiki/Lexicographical_order}{Lexicographical order}
}
\seealso{
\code{\link{comboRank}}, \code{\link{permuteRank}}
}
\author{
Joseph Wood
}
\note{
    \itemize{
        \item {\code{Parallel} and \code{nThreads} will be ignored in the following cases:
            \itemize{
                \item If the class of the vector passed is \code{character} (N.B. \code{Rcpp::CharacterMatrix} is not thread safe). Alternatively, you can generate an indexing matrix in parallel.
                \item If \code{FUN} is utilized.
            }
        }
        \item \code{n} and \code{sampleVec} cannot both be \code{NULL}.
        \item Factor vectors are accepted. Class and level attributes are preserved except when \code{FUN} is used.
    }
}

\examples{
## generate 10 random combinations
comboSample(30, 8, TRUE, n = 5, seed = 10)

## using sampleVec to generate specific permutations
permuteSample(15, 10, freqs = c(1,2,2,1,2,2,1,2,1,2,2,1,2,1,1),
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

## Using nThreads
permPar = permuteSample(10000, 50, n = 8, seed = 10, nThreads = 2)

## Using FUN
permuteSample(10000, 50, n = 4, seed = 10, FUN = sd)

\dontrun{
## Using Parallel
permuteSample(10000, 50, n = 80, seed = 10, Parallel = TRUE)
}
}
\keyword{random}
\keyword{sample}
\keyword{combinations}
\keyword{permutations}
\keyword{combinatorics}
