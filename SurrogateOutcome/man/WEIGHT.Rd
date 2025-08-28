\name{WEIGHT}
\alias{WEIGHT}
\title{
Weight function
}
\description{
Creates a matrix of weights and estimates to be used in estimation of the PTE within R.opt.event
}
\usage{
WEIGHT(xob, deltaob, aob, n, t.0,t)
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
  \item{t.0}{
landmark time of interest, \eqn{t_0}
}
 \item{t}{
time of interest for treatment effect.
}
}

\value{
Returns a matrix of weights and estimates to be used in estimation of the PTE within R.opt.event
}

\author{
Xuan Wang
}

\keyword{internal}
