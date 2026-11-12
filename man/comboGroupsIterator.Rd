\name{comboGroupsIter}
\alias{comboGroupsIter}
\title{
comboGroups Iterator
}
\description{
    \itemize{
        \item Returns an iterator for iterating over partitions of a vector into groups.
        \item Supports random access via the \code{[[} method.
        \item GMP support allows for exploration of cases where the number of comboGroups is large.
        \item Use the \code{next} methods to obtain results in lexicographical order.
    }
}
\usage{
comboGroupsIter(v, numGroups = NULL, grpSizes = NULL,
                retType = "matrix", Parallel = FALSE,
                nThreads = NULL)
}
\arguments{
  \item{v}{Source vector. If \code{v} is a positive integer, it will be converted to the sequence \code{1:v}. If \code{v} is a negative integer, it will be converted to the sequence \code{v:-1}. All atomic types are supported (See \code{\link{is.atomic}}).}
  \item{numGroups}{An Integer. The number of groups that the vector will be partitioned into. The default is \code{NULL}. If provided and \code{grpSize} is \code{NULL}, it must divide the length of v (if v is a vector) or v (if v is a scalar).}
  \item{grpSizes}{A vector of whole numbers representing the size of each group. The default is \code{NULL}. If provided, the sum of the elements must total the length of v (if v is a vector) or v (if v is a scalar).}
  \item{retType}{A string, "3Darray" or "matrix", that determines the shape of the output. The default is "matrix". Note, "3Darray" can only be used when the size of each group is uniform. When the size of each group varies, the return output will always be a matrix.}
  \item{Parallel}{Logical value indicating whether results
  should be generated in parallel. The default is \code{FALSE}. If
  \code{TRUE} and \code{nThreads = NULL}, up to \eqn{n - 1} threads
  are used, where \eqn{n} is the number of threads available on your
  system (both threads may be used on a two-thread system).
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
}

\value{
    \itemize{
        \item If \code{nextIter} is called, a named vector is returned if \code{retType = "matrix"}. If \code{retType = "3Darray"}, a named matrix is returned.
        \item Otherwise a named matrix is returned when \code{retType = "matrix"} and a named 3D array is returned when \code{retType = "3Darray"}.
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
    \item{\code{sourceVector}}{View the source vector}
    \item{\code{summary}}{Returns a list of summary information about the iterator}
    \item{\code{front}}{Retrieve the \bold{first} lexicographical result}
    \item{\code{back}}{Retrieve the \bold{last} lexicographical result}
    \item{\code{[[}}{Random access method. Pass a single value or a vector of valid indices. If a single value is passed, the internal index of the iterator will be updated, however if a vector is passed the internal state will not change. GMP support allows for flexible indexing.}
  }
}
\note{
    \itemize{
        \item{If \code{nThreads} is utilized, it will only take effect if the number of elements requested is greater than some threshold (determined internally). \emph{E.g}:
            \preformatted{serial   <- comboGroupsIter(50, 10)
multi    <- comboGroupsIter(50, 10, nThreads = 4)
fetch1e6 <- multi@nextNIter(1e6)  ## much faster than serial@nextNIter(1e6)
fetch1e3 <- multi@nextNIter(1e3)  ## only one thread used... same as serial@nextNIter(1e3)

library(microbenchmark)
microbenchmark(multi@nextNIter(1e6), serial@nextNIter(1e6), times = 20)
microbenchmark(multi@nextNIter(1e3), serial@nextNIter(1e3), times = 20)}
        }
        \item{The maximum number of comboGroups that can be generated at one time is \eqn{2^{31} - 1}{2^31 - 1}.}
        \item If an error occurs while generating results, the iterator's internal state is not guaranteed to remain at its previous position. Call \code{startOver()} before continuing to use the iterator.
    }
}
\seealso{\link{comboGroups}}
\author{
Joseph Wood
}
\examples{
a = comboGroupsIter(12, 3)
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
\keyword{integer}
\keyword{iterator}
\keyword{comboGroups}
\keyword{combinatorics}
