\name{Constraints-class}
\docType{class}
\title{S4-class for Exposing C++ Constraints Class}
\alias{Constraints-class}
\alias{$,Constraints-method}
\description{
  The \code{Constraints} class is an S4-class that exposes C++ classes that provide access to iterators and other useful methods.
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
  }
}
\author{Joseph Wood}
\seealso{\code{\link{Combo-class}}, \code{\link{Partitions-class}}}
\examples{
  showClass("Constraints")
}
\keyword{classes}
\keyword{iterator}
\keyword{combinations}
\keyword{permutations}
\keyword{combinatorics}
