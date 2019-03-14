\name{stdThreadMax}
\alias{stdThreadMax}
\title{
Max Number of Concurrent Threads
}
\description{
Rcpp wrapper of \href{https://en.cppreference.com/w/cpp/thread/thread/hardware_concurrency}{std::thread::hardware_concurrency()}. As stated by \href{https://en.cppreference.com/w/}{cppreference}, the returned value should be considered only a hint.
}
\usage{
stdThreadMax()
}
%- maybe also 'usage' for other objects documented here.
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An integer representing the number of concurrent threads supported by the user implementation. If the value cannot be determined, \code{1L} is returned.
}

\seealso{
\code{\link[parallel]{detectCores}}
}
\examples{
stdThreadMax()
