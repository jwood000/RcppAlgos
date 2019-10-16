\name{comboGroupsSample}
\alias{comboGroupsSample}
\title{
Sample Partitions of a Vector into Groups of Equal Size
}
\description{
    \itemize{
        \item Generate a specific (lexicographically) or random sample of partitions of groups of equal size.
        \item Produce results in parallel using the \code{Parallel} or \code{nThreads} arguments.
        \item GMP support allows for exploration where the number of results is large.
    }
}
\usage{
comboGroupsSample(v, numGroups, retType = "matrix", n = NULL, 
                  sampleVec = NULL, seed = NULL, Parallel = FALSE,
                  nThreads = NULL, namedSample = FALSE)
}
\arguments{
  \item{v}{Source vector. If \code{v} is an integer (including nonpositive integers), it will be converted to the sequence \code{1:v}. All atomic types are supported (See \code{?is.atomic}).}
  \item{numGroups}{An Integer. The number of groups that the vector will be partitioned into. Must divide the length of v (if v is a vector) or v (if v is a scalar).
}
  \item{retType}{A string, "3Darray" or "matrix", that determines the shape of the output. The default is "matrix".
}
  \item{n}{Number of results to return. The default is \code{NULL}.}
  \item{sampleVec}{A vector of numbers representing the lexicographical partition of groups to return. Accepts vectors of class \code{bigz} as well as vectors of characters}
  \item{seed}{Random seed initialization. The default is \code{NULL}. N.B. If the gmp library is needed, this parameter must be set in order to have reproducible results (\emph{E.g} \code{set.seed()} has no effect in these cases).}
  \item{Parallel}{Logical value indicating whether results should be generated in parallel. The default is \code{FALSE}. If \code{TRUE} and \code{nThreads = NULL}, the number of threads used is equal to the minimum of one minus the number of threads available on your system and the number of results requested (\emph{e.g.} if user has 16 threads and only needs 5 results, 5 threads will be used (\emph{i.e.} \code{min(16 - 1, 5) = 5})). If \code{nThreads} is not \code{NULL}, it will be given preference (\emph{e.g.} if user has 8 threads with \code{Parallel = TRUE} and \code{nThreads = 4}, only 4 threads will be spawned). If your system is single-threaded, the arguments \code{Parallel} and \code{nThreads} are ignored.}
  \item{nThreads}{Specific number of threads to be used. The default is \code{NULL}. See \code{Parallel}.}
  \item{namedSample}{Logical flag. If \code{TRUE}, \code{rownames} corresponding to the lexicographical partitions of groups of equal size, will be added to the returned matrix. The default is \code{FALSE}.}
 }
\details{
These algorithms rely on efficiently generating the \eqn{n^{th}}{nth} lexicographical partition of groups of equal size.
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
## generate 10 random partitions of groups
comboGroupsSample(10, 2, n = 10, seed = 123)

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