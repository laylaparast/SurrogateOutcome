\name{delta.q.event.semi.RMST}
\alias{delta.q.event.semi.RMST}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculates the residual treatment effect (the difference in restricted mean survival time at time t) after accounting for the treatment effect on the surrogate outcome information up to the landmark time
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates the residual treatment effect (the difference in restricted mean survival time at time t) after accounting for the treatment effect on the surrogate outcome information up to the landmark time; uses semi-parametric estimation.
}
\usage{
delta.q.event.semi.RMST(xone, xzero, deltaone, deltazero, sone, szero, t, 
weight = NULL, landmark = landmark, deltaslist = TRUE, number)
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
  \item{sone}{
numeric vector, observed event times for the surrogate outcome in the treatment group.
}
  \item{szero}{
numeric vector, observed event times for the surrogate outcome in the control group.
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
  \item{deltaslist}{
  	TRUE or FALSE; if TRUE, each component of the residual treatment effect is returned along with the residual treatment effect itself, if FALSE, only the residual treatment effect is returned.
%%     ~~Describe \code{deltaslist} here~~
}
  \item{number}{
%%     ~~Describe \code{number} here~~
number of points for RMST calculation, default is 40.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
See documentation for R.q.event for details.
}
\value{
A list is returned:                
 \item{delta.q }{the estimated residual treatment effect}
 \item{first.term}{the first term of the residual treatment effect, if deltaslist = TRUE}
  \item{second.term}{the second term of the residual treatment effect, if deltaslist = TRUE}
   \item{third.term}{the third term of the residual treatment effect, if deltaslist = TRUE}
}

\references{
Parast L, Tian L, and Cai T (2020). Assessing the Value of a Censored Surrogate Outcome. Lifetime Data Analysis, 26(2):245-265.}
\author{
Layla Parast
}
\examples{
data(ExampleData)
names(ExampleData)
\donttest{
delta.q.event.semi.RMST(xone = ExampleData$x1, xzero = ExampleData$x0, deltaone = 
ExampleData$delta1, deltazero = ExampleData$delta0, sone = ExampleData$s1, 
szero = ExampleData$s0, t = 5, landmark=2, number = 40)
delta.q.event.semi.RMST(xone = ExampleData$x1, xzero = ExampleData$x0, deltaone = 
ExampleData$delta1, deltazero = ExampleData$delta0, sone = ExampleData$s1, 
szero = ExampleData$s0, t = 3, landmark=2, number = 40)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}% __ONLY ONE__ keyword per line
\keyword{RMST}% __ONLY ONE__ keyword per line
