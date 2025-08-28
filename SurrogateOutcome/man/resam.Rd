\name{resam}
\alias{resam}
\title{
Resampling for standard error estimation
}
\description{
Resamples data to obtain standard error estimates
}
\usage{
resam(vv, t, t.0, tt, data, data1, data2, indexindex)
}
\arguments{
  \item{vv}{
random exponential weight
}
  \item{t}{
number, the time of interest
}
  \item{t.0}{
landmark time of interest, \eqn{t_0}
}
  \item{tt}{
set of time points between t.0 and t
}
  \item{data}{
full dataset
}
  \item{data1}{
a random half from data
}
  \item{data2}{
the remaining random half from data
}
  \item{indexindex}{
the indices of data contained within data1
}
}
\value{
A vector with the first value being the resampled estimate of the PTE, the second value being the resampled estimate of the optimal g2, and the remaining values being the resampled estimate of the optimal g1(s)
}
\author{
Xuan Wang
}
\keyword{internal}
