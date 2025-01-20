\name{expandGridIter}
\alias{expandGridIter}
\title{
expandGrid Iterator
}
\description{
    \itemize{
        \item Returns an iterator for iterating over the Cartesian product of the input vectors.
        \item Supports random access via the \code{[[} method.
        \item GMP support allows for exploration of cases where the number of products is large.
        \item Use the \code{next} methods to obtain results in lexicographical order.
    }
}
\usage{
expandGridIter(..., nThreads = NULL, return_df = FALSE)
}
\arguments{
  \item{...}{vectors, factors or a list containing these. (See \code{?expand.grid}).}
  \item{nThreads}{Specific number of threads to be used. The default is \code{NULL}.}
  \item{return_df}{Logical flag to force the output to be a \code{data.frame}. The default is \code{FALSE}.}
}

\value{
    \itemize{
        \item If \code{nextIter} is called, a named vector is returned if a \code{matrix} can be returned in the general case. Otherwise, a \code{data.frame} is returned.
        \item When \code{nextNIter} and \code{nextRemaining} are called, a named \code{matrix} is returned when all of the input is of the same type and \code{return_df = FALSE}. Otherwise, a \code{data.frame} is returned.
    }
}
\details{
Once you initialize a new iterator, the following methods are available:
  \describe{
    \item{\code{nextIter}}{Retrieve the \bold{next} lexicographical result}
    \item{\code{nextNIter}}{Pass an integer \emph{n} to retrieve the \bold{next} \emph{n} lexicographical results}
    \item{\code{nextRemaining}}{Retrieve all remaining lexicographical results}
    \item{\code{currIter}}{Returns the current iteration}
    \item{\code{startOver}}{Resets the iterator}
    \item{\code{sourceVector}}{View the source input}
    \item{\code{summary}}{Returns a list of summary information about the iterator}
    \item{\code{front}}{Retrieve the \bold{first} lexicographical result}
    \item{\code{back}}{Retrieve the \bold{last} lexicographical result}
    \item{\code{[[}}{Random access method. Pass a single value or a vector of valid indices. If a single value is passed, the internal index of the iterator will be updated, however if a vector is passed the internal state will not change. GMP support allows for flexible indexing.}
  }
}
\note{
    \itemize{
        \item{If \code{nThreads} is utilized, it will only take effect if the number of elements requested is greater than some threshold (determined internally). \emph{E.g}:
            \preformatted{serial   <- expandGridIter(Map(\(x, y) x:y, 1:10, 11:20))
multi    <- expandGridIter(Map(\(x, y) x:y, 1:10, 11:20), nThreads = 4)
fetch1e6 <- multi@nextNIter(1e6)  ## much faster than serial@nextNIter(1e6)
fetch1e3 <- multi@nextNIter(1e3)  ## only one thread used... same as serial@nextNIter(1e3)

library(microbenchmark)
microbenchmark(multi@nextNIter(1e6), serial@nextNIter(1e6), times = 20)
microbenchmark(multi@nextNIter(1e3), serial@nextNIter(1e3), times = 20)}
        }
        \item{The maximum number of expandGrid that can be generated at one time is \eqn{2^{31} - 1}{2^31 - 1}.}
    }
}
\seealso{\link{expandGrid}}
\author{
Joseph Wood
}
\examples{
a = expandGridIter(factor(state.abb), euro, islands)
a@nextIter()
a@nextNIter(3)
a@front()
all_remaining = a@nextRemaining()
dim(all_remaining)
a@summary()
a@back()
a[[5]]
a@summary()
a[[c(1, 17, 3)]]
a@summary()
}
\keyword{expand}
\keyword{product}
\keyword{iterator}
\keyword{cartesian}
\keyword{expandGrid}
\keyword{expand.grid}
\keyword{combinatorics}
