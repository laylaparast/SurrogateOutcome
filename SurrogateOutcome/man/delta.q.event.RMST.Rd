\name{delta.q.event.RMST}
\alias{delta.q.event.RMST}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the residual treatment effect (the difference in restricted mean survival time at time t) after accounting for the treatment effect on the surrogate outcome information up to the landmark time
}
\description{
Calculates the residual treatment effect (the difference in restricted mean survival time at time t) after accounting for the treatment effect on the surrogate outcome information up to the landmark time; uses nonparametric estimation.
}
\usage{
delta.q.event.RMST(xone, xzero, deltaone, deltazero, sone, szero, t, weight = NULL,
 landmark = landmark, deltaslist = TRUE, transform = FALSE, extrapolate=TRUE, 
 number, warn.extrapolate=TRUE)
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
  \item{transform}{
%%     ~~Describe \code{transform} here~~
TRUE or FALSE; indicates whether a transformation should be used, default is FALSE.
}
  \item{extrapolate}{
%%     ~~Describe \code{extrapolate} here~~
TRUE or FALSE; indicates whether local constant extrapolation should be used, default is TRUE.
}
  \item{number}{
%%     ~~Describe \code{number} here~~
number of points for RMST calculation, default is 40.
}
  \item{warn.extrapolate}{
%%     ~~Describe \code{extrapolate} here~~
TRUE or FALSE; indicates whether user prefers a warning message when extrapolation is used, default is TRUE.
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
delta.q.event.RMST(xone = ExampleData$x1, xzero = ExampleData$x0, deltaone = ExampleData$delta1, 
deltazero = ExampleData$delta0, sone = ExampleData$s1, szero = ExampleData$s0, t = 5, landmark=2, 
number = 40)
delta.q.event.RMST(xone = ExampleData$x1, xzero = ExampleData$x0, deltaone = ExampleData$delta1, 
deltazero = ExampleData$delta0, sone = ExampleData$s1, szero = ExampleData$s0, t = 4, landmark=2, 
number = 40)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}% use one of  RShowDoc("KEYWORDS")
\keyword{survival}% __ONLY ONE__ keyword per line
\keyword{RMST}% __ONLY ONE__ keyword per line
