\name{delta.estimate.RMST}
\alias{delta.estimate.RMST}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Helper function
}
\description{
Helper function; used by delta.estimate function
}
\usage{
delta.estimate.RMST(xone, xzero, deltaone, deltazero, t, weight = NULL, delta.only = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xone}{
numeric vector, observed event times for the primary outcome in the treatment group.
}
  \item{xzero}{
numeric vector, observed event times for the primary outcome in the control group.
}
  \item{deltaone}{
numeric vector, event/censoring indicators for the primary outcome in the treatment group.
}
  \item{deltazero}{
numeric vector, event/censoring indicators for the primary outcome in the control group.
}
  \item{t}{
time of interest for treatment effect.
}
  \item{weight}{
	weight; optional   
}
  \item{delta.only}{
TRUE or FALSE; if TRUE then only delta is returned, if FALSE then delta and restricted mean survival time for each treatment group is also returned
}
}

\value{
A list is returned:                
 \item{delta }{the difference in the restricted mean survival between the two treatment groups. }
 \item{s1}{only if delta.only = F; the restricted mean survival in group 1.}
 \item{s0}{only if delta.only = F; the restricted mean survival in group 0.}
}

\author{
%%  ~~who you are~~
Layla Parast
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
