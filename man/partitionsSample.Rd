\name{partitionsSample}
\alias{partitionsSample}
\title{
Sample Partitions
}
\description{
    \itemize{
        \item Generate a specific (lexicographically) or random sample of partitions.
        \item Produce results in parallel using the \code{Parallel} or \code{nThreads} arguments.
        \item GMP support allows for exploration of partitions of vectors with many elements.
    }
}
\usage{
partitionsSample(v, m = NULL, repetition = FALSE, freqs = NULL,
                 target = NULL, n = NULL, sampleVec = NULL,
                 seed = NULL, Parallel = FALSE,
                 nThreads = NULL, namedSample = FALSE)
}
\arguments{
  \item{v}{Source vector. If \code{v} is a positive integer, it will be converted to the sequence \code{1:v}. If \code{v} is a negative integer, it will be converted to the sequence \code{v:-1}. Only integer and numeric vectors are accepted.}
  \item{m}{Number of elements to choose.  If \code{repetition = TRUE} or \code{freqs} is utilized, \code{m} can exceed the length of \code{v}. If \code{m = NULL}, the length will default to \code{length(v)} or \code{sum(freqs)}.}
  \item{repetition}{Logical value indicating whether partitions should be with or without repetition. The default is \code{FALSE}.}
  \item{freqs}{A vector of frequencies used for producing all partitions of a multiset of \code{v}. Each element of \code{freqs} represents how many times each element of the source vector, \code{v}, is repeated. It is analogous to the \code{times} argument in \code{\link[base]{rep}}. The default value is \code{NULL}.}
  \item{target}{Number to be partitioned.}
  \item{n}{Number of partitions to return. The default is \code{NULL}.}
  \item{sampleVec}{A vector of numbers representing the lexicographical partitions to return. Accepts vectors of class \code{bigz} as well as vectors of characters}
  \item{seed}{Random seed initialization. The default is \code{NULL}. N.B. If the gmp library is needed, this parameter must be set in order to have reproducible results (\emph{E.g} \code{set.seed()} has no effect in these cases).}
  \item{Parallel}{Logical value indicating whether partitions should be generated in parallel. The default is \code{FALSE}. If \code{TRUE} and \code{nThreads = NULL}, the number of threads used is equal to the minimum of one minus the number of threads available on your system and the number of results requested (\emph{e.g.} if user has 16 threads and only needs 5 results, 5 threads will be used (\emph{i.e.} \code{min(16 - 1, 5) = 5})). If \code{nThreads} is not \code{NULL}, it will be given preference (\emph{e.g.} if user has 8 threads with \code{Parallel = TRUE} and \code{nThreads = 4}, only 4 threads will be spawned). If your system is single-threaded, the arguments \code{Parallel} and \code{nThreads} are ignored.}
  \item{nThreads}{Specific number of threads to be used. The default is \code{NULL}. See \code{Parallel}.}
  \item{namedSample}{Logical flag. If \code{TRUE}, \code{rownames} corresponding to the lexicographical partition, will be added to the returned matrix. The default is \code{FALSE}.}
 }
\details{
These algorithms rely on efficiently generating the \eqn{n^{th}}{nth} lexicographical partition. This is the process of \href{https://rosettacode.org/wiki/Permutations/Rank_of_a_permutation}{unranking}.
}
\value{
In general, a matrix is returned with each row containing a vector of length \eqn{m}. If \code{FUN} is utilized, a list is returned.
}
\references{
\href{https://en.wikipedia.org/wiki/Lexicographical_order}{Lexicographical order}
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
        \item Factor vectors are accepted. Class and level attributes are preserved.
    }
}

\examples{
partitionsSample(100, 10, n = 10)
}