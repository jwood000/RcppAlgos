\name{Combo-class}
\docType{class}
\title{S4-classes for Exposing C++ Combinatorial Classes}
\alias{Combo-class}
\alias{ComboApply-class}
\alias{ComboRes-class}
%
\alias{$,Combo-method}
\alias{$,ComboApply-method}
\alias{$,ComboRes-method}
%
\description{
  The \code{Combo} class family are S4-classes that expose C++ classes that provide access to iterators and other useful methods.
}
\section{Slots}{
  of \code{"Combo"} and all classes inheriting from it:
  \describe{
    \item{\code{nextIter}}{Retrieve the \bold{next} lexicographical result}
    \item{\code{nextNIter}}{Pass an integer \emph{n} to retrieve the \bold{next} \emph{n} lexicographical results}
    \item{\code{nextRemaining}}{Retrieve all remaining lexicographical results}
    \item{\code{currIter}}{Returns the current iteration}
    \item{\code{prevIter}}{Retrieve the \bold{previous} lexicographical result (the \bold{next} \emph{reverse} lexicographical result)}
    \item{\code{prevNIter}}{Pass an integer \emph{n} to retrieve the \bold{previous} \emph{n} lexicographical results (the \bold{next} \emph{n} \emph{reverse} lexicographical results)}
    \item{\code{prevRemaining}}{Retrieve all remaining \emph{reverse} lexicographical results}
    \item{\code{startOver}}{Resets the iterator}
    \item{\code{sourceVector}}{View the source vector}
    \item{\code{summary}}{Returns a list of summary information about the iterator}
    \item{\code{front}}{Retrieve the \bold{first} lexicographical result}
    \item{\code{back}}{Retrieve the \bold{last} lexicographical result}
    \item{\code{randomAccess}}{Random access method. Pass a single value or a vector of valid indices. If a single value is passed, the internal index of the iterator will be updated, however if a vector is passed the internal state will not change. GMP support allows for flexible indexing.}
  }
}

\author{Joseph Wood}
\seealso{\code{\link{Partitions-class}}, \code{\link{Constraints-class}}}
\examples{
  showClass("Combo")
}
\keyword{classes}
\keyword{combinations}
\keyword{permutations}
\keyword{combinatorics}
