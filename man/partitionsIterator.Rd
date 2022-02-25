\name{partitionsIter}
\alias{partitionsIter}
\title{
Partition Iterator
}
\description{
    \itemize{
        \item Returns an iterator for iterating over combinations or permutations of a vector.
        \item Supports random access via the \code{[[} method.
        \item GMP support allows for exploration of combinations/permutations of vectors with many elements.
        \item The output is in lexicographical order for the \code{next} methods and reverse lexicographical order for the \code{prev} methods.
    }
}
\usage{
partitionsIter(v, m = NULL, repetition = FALSE,
               freqs = NULL, target = NULL, 
               nThreads = NULL, tolerance = NULL)
}
\arguments{
  \item{v}{Source vector. If \code{v} is a positive integer, it will be converted to the sequence \code{1:v}. If \code{v} is a negative integer, it will be converted to the sequence \code{v:-1}. Only integer and numeric vectors are accepted.}
  \item{m}{Number of elements to choose.  If \code{repetition = TRUE} or \code{freqs} is utilized, \code{m} can exceed the length of \code{v}. If \code{m = NULL}, the length will default to \code{length(v)} or \code{sum(freqs)}.}
  \item{repetition}{Logical value indicating whether combinations/permutations should be with or without repetition. The default is \code{FALSE}.}
  \item{freqs}{A vector of frequencies used for producing all combinations/permutations of a multiset of \code{v}. Each element of \code{freqs} represents how many times each element of the source vector, \code{v}, is repeated. It is analogous to the \code{times} argument in \code{\link[base]{rep}}. The default value is \code{NULL}.}
  \item{target}{Number to be partitioned.}
  \item{nThreads}{Specific number of threads to be used. The default is \code{NULL}.}
  \item{tolerance}{A numeric value greater than or equal to zero. This parameter is utilized when a constraint is applied on a numeric vector. The default value is 0 when it can be determined that whole values are being utilized, otherwise it is \code{sqrt(.Machine$double.eps)} which is approximately \eqn{1.5e-8}. N.B. If the input vector is of type integer, this parameter will be ignored and strict equality will be enforced.}
}

\value{
A vector of length \emph{m} is returned if \code{nextIter} or \code{prevIter} is called, otherwise a matrix is returned with each row containing a vector of length \eqn{m}. If \code{FUN} is utilized and either \code{nextIter} or \code{prevIter} is called, the result will be determined by \code{FUN}, otherwise a list is returned.
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
        \item The maximum number of combinations/permutations that can be generated at one time is \eqn{2^{31} - 1}{2^31 - 1}.
        \item Factor vectors are accepted. Class and level attributes are preserved.
        \item Currently, constrained iterators are not supported (Will be available in future releases).
    }
}
\seealso{
\code{\link{partitionsGeneral}}
}
\references{
    \itemize{
        \item{\href{https://en.wikipedia.org/wiki/Lexicographical_order}{Lexicographical Order}}
    }
}
\author{
Joseph Wood
}
\examples{
a = partitionsIter(0:10, repetition = TRUE)
a@nextIter()
a@nextNIter(3)
a@front()
a@nextRemaining()
a@summary()
a@back()
a[[5]]
a@summary()
a[[c(1, 17, 3)]]
a@summary()
}