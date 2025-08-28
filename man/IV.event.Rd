\name{IV.event}
\alias{IV.event}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the incremental value of the surrogate outcome information 
}
\description{
Calculates the incremental value of the surrogate outcome information 
}
\usage{
IV.event(xone, xzero, deltaone, deltazero, sone, szero, t, landmark, number = 40, 
transform = FALSE, extrapolate = TRUE, std = FALSE, conf.int = FALSE, 
weight.perturb = NULL, type = "np")
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
  \item{landmark}{
landmark time of interest, \eqn{t_0}.
}
  \item{number}{
%%     ~~Describe \code{number} here~~
number of points for RMST calculation, default is 40.
}
  \item{transform}{
%%     ~~Describe \code{transform} here~~
TRUE or FALSE; indicates whether a transformation should be used, default is FALSE.
}
  \item{extrapolate}{
%%     ~~Describe \code{extrapolate} here~~
TRUE or FALSE; indicates whether local constant extrapolation should be used, default is FALSE.
}
  \item{std}{
TRUE or FALSE; indicates whether standard error estimates should be provided, default is FALSE. Estimates are calculated using perturbation-resampling. Two versions are provided: one that takes the standard deviation of the perturbed estimates (denoted as "sd") and one that takes the median absolute deviation (denoted as "mad"). 
}
  \item{conf.int}{
TRUE or FALSE; indicates whether 95\% confidence intervals should be provided. Confidence intervals are calculated using the percentiles of perturbed estimates, default is FALSE. If this is TRUE, standard error estimates are automatically provided.
}
  \item{weight.perturb}{
weights used for perturbation resampling.
}
  \item{type}{
Type of estimate that should be provided; options are "np" for the nonparametric estimate or "semi" for the semiparametric estimate, default is "np".
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The incremental value of the surrogate outcome information only is quantified as \eqn{IV_S(t,t_0) = R_Q(t,t_0) - R_T(t,t_0)} where the definition and estimation procedures for \eqn{R_Q(t,t_0)} and \eqn{R_T(t,t_0)} are described in the documentation for R.q.event and R.t.estimate, respectively. The estimate of the incremental value is  \eqn{\hat{IV}_S(t,t_0) = \hat{R}_Q(t,t_0) - \hat{R}_T(t,t_0)}. 

}
\value{
A list is returned:                
 \item{delta }{the estimate, \eqn{\hat{\Delta}(t)}, described in delta.estimate documentation.}
 \item{delta.q }{the estimate, \eqn{\hat{\Delta}_Q(t,t_0)}, described in R.q.event documention.}
 \item{R.q }{the estimate, \eqn{\hat{R}_Q(t,t_0)}, described in R.q.event documention.}
  \item{delta.t }{the estimate, \eqn{\hat{\Delta}_T(t,t_0)}, described in R.t.estimate documention.}
 \item{R.t }{the estimate, \eqn{\hat{R}_T(t,t_0)}, described in R.t.estimate documention.}
 \item{IV}{the estimated incremental value of the surrogate outcome information, described above.}
 \item{delta.sd }{the standard error estimate of \eqn{\hat{\Delta}(t)}; if std = TRUE or conf.int = TRUE.}
 \item{delta.mad }{the standard error estimate of \eqn{\hat{\Delta}(t)} using the median absolute deviation; if std = TRUE or conf.int = TRUE.} 
  \item{delta.q.sd}{the standard error estimate of \eqn{\hat{\Delta}_Q(t,t_0)}; if std = TRUE or conf.int = TRUE.}
  \item{delta.q.mad}{the standard error estimate of \eqn{\hat{\Delta}_Q(t,t_0)} using the median absolute deviation; if std = TRUE or conf.int = TRUE.}  
 \item{R.q.sd }{the standard error estimate of \eqn{\hat{R}_Q(t,t_0)}; if std = TRUE or conf.int = TRUE.}
  \item{R.q.mad }{the standard error estimate of \eqn{\hat{R}_Q(t,t_0)} using the median absolute deviation; if std = TRUE or conf.int = TRUE.}
\item{delta.t.sd}{the standard error estimate of \eqn{\hat{\Delta}_T(t,t_0)}; if std = TRUE or conf.int = TRUE.}
  \item{delta.t.mad}{the standard error estimate of \eqn{\hat{\Delta}_T(t,t_0)} using the median absolute deviation; if std = TRUE or conf.int = TRUE.}  
 \item{R.t.sd }{the standard error estimate of \eqn{\hat{R}_T(t,t_0)}; if std = TRUE or conf.int = TRUE.}
  \item{R.t.mad }{the standard error estimate of \eqn{\hat{R}_T(t,t_0)} using the median absolute deviation; if std = TRUE or conf.int = TRUE.}
   \item{IV.sd }{the standard error estimate of the incremental value; if std = TRUE or conf.int = TRUE.}
 \item{IV.mad }{the standard error estimate of the incremental value using the median absolute deviation; if std = TRUE or conf.int = TRUE.}
  \item{conf.int.delta}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}(t)} based on sample quantiles of the perturbed values; if conf.int = TRUE.}
 \item{conf.int.delta.q}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}_Q(t,t_0)} based on sample quantiles of the perturbed values; if conf.int = TRUE.} 
 \item{conf.int.R.q}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_Q(t,t_0)} based on sample quantiles of the perturbed values; if conf.int = TRUE.}
  \item{conf.int.delta.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}_T(t,t_0)} based on sample quantiles of the perturbed values; if conf.int = TRUE.} 
 \item{conf.int.R.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_T(t,t_0)} based on sample quantiles of the perturbed values; if conf.int = TRUE.}
 \item{conf.int.IV}{a vector of size 2; the 95\% confidence interval for the incremental value based on sample quantiles of the perturbed values; if conf.int = TRUE.}
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
IV.event(xone = ExampleData$x1, xzero = ExampleData$x0, deltaone = ExampleData$delta1, 
deltazero = ExampleData$delta0, sone = ExampleData$s1, szero = ExampleData$s0, t = 5,
 landmark=2, type = "np")

IV.event(xone = ExampleData$x1, xzero = ExampleData$x0, deltaone = ExampleData$delta1, 
deltazero = ExampleData$delta0, sone = ExampleData$s1, szero = ExampleData$s0, t = 5,
 landmark=2, type = "np", std = TRUE, conf.int = TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}% use one of  RShowDoc("KEYWORDS")
\keyword{survival}% __ONLY ONE__ keyword per line
