\name{expandGridSample}
\alias{expandGridSample}
\title{
Sample the Cartesian Product
}
\description{
    \itemize{
        \item Generate a specific (lexicographically) or random sample of the Cartesian product of the input vectors.
        \item Produce results in parallel using the \code{nThreads} arguments.
        \item GMP support allows for exploration where the number of results is large.
    }
}
\usage{
expandGridSample(
    ..., n = NULL, sampleVec = NULL, seed = NULL,
    nThreads = NULL, namedSample = FALSE, return_df = FALSE
)
}
\arguments{
  \item{...}{vectors, factors or a list containing these. (See \code{?expand.grid}).}
  \item{n}{Number of results to return. The default is \code{NULL}.}
  \item{sampleVec}{A vector of numbers representing the lexicographical partition of groups to return. Accepts vectors of class \code{bigz} as well as vectors of characters}
  \item{seed}{Random seed initialization. The default is \code{NULL}. N.B. If the gmp library is needed, this parameter must be set in order to have reproducible results (\emph{E.g} \code{set.seed()} has no effect in these cases).}
  \item{nThreads}{Specific number of threads to be used. The default is \code{NULL}.}
  \item{namedSample}{Logical flag. If \code{TRUE}, \code{rownames} corresponding to the lexicographical result, will be added to the returned matrix. The default is \code{FALSE}.}
  \item{return_df}{Logical flag to force the output to be a \code{data.frame}. The default is \code{FALSE}.}
}
\details{
These algorithms rely on efficiently generating the \eqn{n^{th}}{nth} lexicographical result.
}
\value{
When all of the input is of the same type, by default \code{expandGridSample} produces a \code{matrix} (a \code{data.frame} otherwise). This can be ignored by setting the argument \code{return_df = TRUE}.
}
\references{
\href{https://en.wikipedia.org/wiki/Lexicographical_order}{Lexicographical order}
}
\author{
Joseph Wood
}
\examples{
## input vectors
lst = list(factor(state.abb), euro, islands)

## generate 10 random products
expandGridSample(lst, n = 10, seed = 123)

## using sampleVec to generate specific results
expandGridSample(lst, sampleVec = c(1, 100, 1e3))

all.equal(expandGridSample(lst, sampleVec = 1:expandGridCount(lst)),
         expandGrid(lst))

## Examples with enormous number of total results
big_lst = Map(function(x, y) x:y, 8:33, 15:40)
num = expandGridCount(big_lst)
gmp::log2.bigz(num)
## [1] 78

first = gmp::urand.bigz(n = 1, size = 78, seed = 123)
mySamp = do.call(c, lapply(0:10, function(x) gmp::add.bigz(first, x)))

class(mySamp)
## [1] "bigz"

## using the sampling function
cartSamp = expandGridSample(big_lst, sampleVec = mySamp)

## using the standard function
cartGeneral = expandGrid(big_lst,
                         lower = first,
                         upper = gmp::add.bigz(first, 10))

identical(cartSamp, cartGeneral)
## [1] TRUE
}
\keyword{rank}
\keyword{unrank}
\keyword{random}
\keyword{sample}
\keyword{expand}
\keyword{product}
\keyword{cartesian}
\keyword{expandGrid}
\keyword{expand.grid}
\keyword{combinatorics}
