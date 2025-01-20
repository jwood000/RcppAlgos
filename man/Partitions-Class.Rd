\name{Partitions-class}
\docType{class}
\title{S4-class for Exposing C++ Partitions Class}
\alias{Partitions-class}
\alias{$,Partitions-method}
\description{
  The \code{Partitions} class is an S4-class that exposes C++ classes that provide access to iterators and other useful methods.
}
\section{Slots}{
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
    \item{\code{randomAccess}}{Random access method. Pass a single value or a vector of valid indices. If a single value is passed, the internal index of the iterator will be updated, however if a vector is passed the internal state will not change. GMP support allows for flexible indexing.}
  }
}
\author{Joseph Wood}
\seealso{\code{\link{Combo-class}}, \code{\link{Constraints-class}}}
\examples{
  showClass("Partitions")
}
\keyword{classes}
\keyword{integer}
\keyword{iterator}
\keyword{partitions}
\keyword{combinatorics}
