\name{delta.t.RMST}
\alias{delta.t.RMST}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculates the residual treatment effect (the difference in restricted mean survival time at time t) after accounting for the treatment effect on the primary outcome up to the landmark time

}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates the residual treatment effect (the difference in restricted mean survival time at time t) after accounting for the treatment effect on the primary outcome up to the landmark time
}
\usage{
delta.t.RMST(xone, xzero, deltaone, deltazero, t, weight = NULL, landmark = landmark)
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
  	optional weight.
%%     ~~Describe \code{weight} here~~
}
  \item{landmark}{
landmark time of interest, \eqn{t_0}.
}

}
\details{
See documentation for R.t.estimate for details.
}
\value{
 \item{delta.t}{the estimated residual treatment effect after accounting for the treatment effect on the primary outcome up to the landmark time
}

}
\references{
Parast L, Tian L, and Cai T (2020). Assessing the Value of a Censored Surrogate Outcome. Lifetime Data Analysis, 26(2):245-265.}
\author{
Layla Parast
}


\examples{
data(ExampleData)
names(ExampleData)
delta.t.RMST(xone = ExampleData$x1, xzero = ExampleData$x0, deltaone = ExampleData$delta1, 
deltazero = ExampleData$delta0, t = 5, landmark=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}% use one of  RShowDoc("KEYWORDS")
\keyword{survival}% __ONLY ONE__ keyword per line
\keyword{RMST}% __ONLY ONE__ keyword per line
