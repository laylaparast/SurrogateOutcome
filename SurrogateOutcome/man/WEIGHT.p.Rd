\name{WEIGHT.p}
\alias{WEIGHT.p}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weight function for resampling
}
\description{
Creates a matrix of weights and estimates to be used in estimation of the standard error of the PTE within resam
}
\usage{
WEIGHT.p(xob, deltaob, aob, n, v, t.0, t)
}
\arguments{
   \item{xob}{
numeric vector, the observed event time: X = min(T, C) where T is the time of the primary outcome, C is the censoring time
}
  \item{deltaob}{
numeric vector of 0/1, the censoring indicator: D = I(T<C) where T is the time of the primary outcome, C is the censoring time

}
  \item{aob}{
numeric vector of 0/1, indicating treatment group

}
  \item{n}{
sample size
}
  \item{v}{
random exponential weight
}
  \item{t.0}{
landmark time of interest, \eqn{t_0}
}
 \item{t}{
time of interest for treatment effect.
}
}
\value{
Returns a matrix of weights and estimates to be used in estimation of the PTE within resam

}
\author{
Xuan Wang
}

\keyword{internal}