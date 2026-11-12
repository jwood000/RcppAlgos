\name{comboGroupsSample}
\alias{comboGroupsSample}
\title{
Sample Partitions of a Vector into Groups
}
\description{
    \itemize{
        \item Generate a specific (lexicographically) or random sample of partitions of groups.
        \item Produce results in parallel using the \code{Parallel} or \code{nThreads} arguments.
        \item GMP support allows for exploration where the number of results is large.
    }
}
\usage{
comboGroupsSample(v, numGroups = NULL, grpSizes = NULL, retType = "matrix",
                  n = NULL, sampleVec = NULL, seed = NULL, Parallel = FALSE,
                  nThreads = NULL, namedSample = FALSE)
}
\arguments{
  \item{v}{Source vector. If \code{v} is a positive integer, it will be converted to the sequence \code{1:v}. If \code{v} is a negative integer, it will be converted to the sequence \code{v:-1}. All atomic types are supported (See \code{\link{is.atomic}}).}
  \item{numGroups}{An Integer. The number of groups that the vector will be partitioned into. The default is \code{NULL}. If provided and \code{grpSize} is \code{NULL}, it must divide the length of v (if v is a vector) or v (if v is a scalar).}
  \item{grpSizes}{A vector of whole numbers representing the size of each group. The default is \code{NULL}. If provided, the sum of the elements must total the length of v (if v is a vector) or v (if v is a scalar).}
  \item{retType}{A string, "3Darray" or "matrix", that determines the shape of the output. The default is "matrix". Note, "3Darray" can only be used when the size of each group is uniform. When the size of each group varies, the return output will always be a matrix.}
  \item{n}{Number of results to return. The default is \code{NULL}.}
  \item{sampleVec}{A vector of numbers representing the lexicographical partition of groups to return. Accepts vectors of class \code{bigz} as well as vectors of characters}
  \item{seed}{Random seed initialization. The default is \code{NULL}. N.B. If the gmp library is needed, this parameter must be set in order to have reproducible results (\emph{E.g} \code{set.seed()} has no effect in these cases).}
  \item{Parallel}{Logical value indicating whether results
  should be generated in parallel. The default is \code{FALSE}. If
  \code{TRUE} and \code{nThreads = NULL}, up to \eqn{n - 1} threads
  are used, where \eqn{n} is the number of threads available on your
  system (both threads may be used on a two-thread system).
  This is also capped by the number of results requested (\emph{e.g.}
  with 16 available threads and 5 results requested, up to 5 threads
  may be used, since \code{min(16 - 1, 5) = 5}).
  If \code{nThreads} is not \code{NULL}, it will be given preference
  (\emph{e.g.} with 8 available threads, \code{Parallel = TRUE} and
  \code{nThreads = 4}, up to 4 threads may be used). In either case,
  the actual number of threads depends on the useful work in the
  calculation; small or simple cases may use fewer threads.
  If your system is single-threaded, the arguments \code{Parallel}
  and \code{nThreads} are ignored.}
  \item{nThreads}{Requested number of threads. The default is
  \code{NULL} for serial execution unless \code{Parallel = TRUE}.
  Requests are capped by the available hardware and the useful work in
  the calculation. Small or simple cases may use fewer threads.
  See \code{Parallel}.}
  \item{namedSample}{Logical flag. If \code{TRUE}, \code{rownames} corresponding to the lexicographical result, will be added to the returned matrix. The default is \code{FALSE}.}
 }
\details{
These algorithms rely on efficiently generating the \eqn{n^{th}}{nth} lexicographical result.
}
\value{
By default, a matrix is returned with column names corresponding to the associated group. If \code{retType = "3Darray"}, a 3D array is returned.
}
\references{
\href{https://en.wikipedia.org/wiki/Lexicographical_order}{Lexicographical order}
}
\author{
Joseph Wood
}
\examples{
## generate 10 random partitions of groups of equal size
comboGroupsSample(10, 2, n = 10, seed = 123)

## generate 10 random partitions of groups of varying sizes
comboGroupsSample(10, grpSizes = 1:4, n = 10, seed = 123)

## using sampleVec to generate specific results
comboGroupsSample(15, 5, sampleVec = c(1, 100, 1e3, 1e6))

all.equal(comboGroupsSample(10, 5,
            sampleVec = 1:comboGroupsCount(10, 5)),
         comboGroups(10, 5))

## Examples with enormous number of total results
num = comboGroupsCount(100, 20)
gmp::log2.bigz(num)
## [1] 325.5498

first = gmp::urand.bigz(n = 1, size = 325, seed = 123)
mySamp = do.call(c, lapply(0:10, function(x) gmp::add.bigz(first, x)))

class(mySamp)
## [1] "bigz"

## using the sampling function
cbgSamp = comboGroupsSample(100, 20, sampleVec = mySamp)

## using the standard function
cbgGeneral = comboGroups(100, 20,
                         lower = first,
                         upper = gmp::add.bigz(first, 10))

identical(cbgSamp, cbgGeneral)
## [1] TRUE

\dontrun{
## Using Parallel
system.time(comboGroupsSample(1000, 20, n = 80, seed = 10, Parallel = TRUE))
}
}
\keyword{random}
\keyword{sample}
\keyword{rank}
\keyword{unrank}
\keyword{partitions}
\keyword{comboGroups}
\keyword{combinatorics}
