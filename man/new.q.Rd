\name{new.q}
\alias{new.q}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Helper function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Helper function; should not be called directly by user. This function differs from the quantile function in that it returns an NA when there are NAs in the supplied vector, rather than causing an error. 
}
\usage{
new.q(x, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ vector of numeric values
%%     ~~Describe \code{x} here~~
}
  \item{p}{ sample quantile that is desired
%%     ~~Describe \code{p} here~~
}
}
\value{
returns either NA or the desired sample quantile
%% ...
}
\author{
%%  ~~who you are~~
Layla Parast
}
\keyword{arith}
\keyword{internal}

