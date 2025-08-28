\name{R.t.estimate}
\alias{R.t.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the proportion of the treatment effect (the difference in restricted mean survival time at time t) explained by primary outcome information observed up to the landmark time
}
\description{
Calculates the proportion of the treatment effect (the difference in restricted mean survival time at time t) explained by primary outcome information observed up to the landmark time; also provides standard error estimate and confidence interval.
}
\usage{
R.t.estimate(xone, xzero, deltaone, deltazero, t, landmark, std = FALSE, conf.int
 = FALSE, weight.perturb = NULL)
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
  \item{landmark}{
landmark time of interest, \eqn{t_0}.
}
  \item{std}{
TRUE or FALSE; indicates whether standard error estimates should be provided, default is FALSE. Estimates are calculated using perturbation-resampling. Two versions are provided: one that takes the standard deviation of the perturbed estimates (denoted as "sd") and one that takes the median absolute deviation (denoted as "mad"). 
}
  \item{conf.int}{
TRUE or FALSE; indicates whether 95\% confidence intervals should be provided. Confidence intervals are calculated using the percentiles of perturbed estimates, default is FALSE. If this is TRUE, standard error estimates are automatically provided.
}
  \item{weight.perturb}{
weights used for perturbation resampling.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Let \eqn{G \in \{1,0\}} be the randomized treatment indicator, \eqn{T} denote the time of the primary outcome of interest, and \eqn{S} denote the time of the surrogate outcome. We use potential outcomes notation such that \eqn{T^{(G)}} and \eqn{S^{(G)}} denote the respective times of the primary and surrogate outcomes under treatment G, for \eqn{G \in \{1, 0\}}.  In the absence of censoring, we only observe \eqn{(T, S)=(T^{(1)}, S^{(1)})} or \eqn{(T^{(0)}, S^{(0)})} for each individual depending on whether \eqn{G=1} or \eqn{0}. Due to censoring, data consist of \eqn{n = n_1 + n_0} independent observations \eqn{\{X_{gi}, \delta_{gi}, I(S_{gi}< t_0)I(X_{gi} > t_0),  S_{gi}\wedge t_0 I(X_{gi} > t_0), i=1,...,n_g, g = 1,0\}}, where \eqn{X_{gi} = T_{gi}\wedge C_{ gi}}, \eqn{\delta_{gi} = I(T_{gi} < C_{gi})}, \eqn{C_{gi}} denotes the censoring time, \eqn{T_{gi}} denotes the time of the primary outcome, \eqn{S_{gi}} denotes the time of the surrogate outcome, \eqn{\{(T_{gi}, C_{gi}, S_{gi}), i = 1, ..., n_g\}} are identically distributed within treatment group, and \eqn{t_0} is the landmark time of interest.

The proportion of treatment effect explained by primary outcome information observed up to the landmark time, \eqn{t_0}, is defined as \eqn{R_T(t, t_0) = 1-\Delta_T(t,t_0) / \Delta(t)} where 
\deqn{\Delta_T(t,t_0)   =  P(T^{(0)} > t_0) E\{T^{(1)} \wedge t - T^{(0)} \wedge t \mid T > t_0 \}} and \eqn{\Delta(t)} is the treatment effect on the primary outcome, defined in the documentation for delta.estimate. The quantity \eqn{\Delta_T(t,t_0)} is estimated using \deqn{\hat{\Delta}_T(t,t_0) = n_0^{-1} \sum_{i=1}^{n_0} I(X_{0i} > t_0) /\hat{W}_0^C(t_0)  \{\hat{\nu}_{1}(t|t_0) - \hat{\nu}_{0}(t|t_0) \}} where \eqn{\hat{W}^C_0(t)} is the Kaplan-Meier estimator of \eqn{P(C_{gi} \ge t)}, \eqn{\hat{\nu}_{0}(t|t_0)} is defined in the documentation for R.q.event and \eqn{\hat{\nu}_{1}(t|t_0)} is obtained by replacing \eqn{0} with \eqn{1}.
}
\value{
A list is returned:                
 \item{delta }{the estimate, \eqn{\hat{\Delta}(t)}, described in delta.estimate documentation.}
 \item{delta.t }{the estimate, \eqn{\hat{\Delta}_T(t,t_0)}, described above.}
 \item{R.t }{the estimate, \eqn{\hat{R}_T(t,t_0)}, described above.}
 \item{delta.sd }{the standard error estimate of \eqn{\hat{\Delta}(t)}; if std = TRUE or conf.int = TRUE.}
 \item{delta.mad }{the standard error estimate of \eqn{\hat{\Delta}(t)} using the median absolute deviation; if std = TRUE or conf.int = TRUE.} 
  \item{delta.t.sd}{the standard error estimate of \eqn{\hat{\Delta}_T(t,t_0)}; if std = TRUE or conf.int = TRUE.}
  \item{delta.t.mad}{the standard error estimate of \eqn{\hat{\Delta}_T(t,t_0)} using the median absolute deviation; if std = TRUE or conf.int = TRUE.}  
 \item{R.t.sd }{the standard error estimate of \eqn{\hat{R}_T(t,t_0)}; if std = TRUE or conf.int = TRUE.}
  \item{R.t.mad }{the standard error estimate of \eqn{\hat{R}_T(t,t_0)} using the median absolute deviation; if std = TRUE or conf.int = TRUE.}
  \item{conf.int.delta}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}(t)} based on sample quantiles of the perturbed values; if conf.int = TRUE.}
 \item{conf.int.delta.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}_T(t,t_0)} based on sample quantiles of the perturbed values; if conf.int = TRUE.} 
 \item{conf.int.R.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_T(t,t_0)} based on sample quantiles of the perturbed values; if conf.int = TRUE.}

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
R.t.estimate(xone = ExampleData$x1, xzero = ExampleData$x0, deltaone =ExampleData$delta1, 
deltazero = ExampleData$delta0, t = 5, landmark=2)

R.t.estimate(xone = ExampleData$x1, xzero = ExampleData$x0, deltaone = ExampleData$delta1, 
deltazero = ExampleData$delta0, t = 5, landmark=2, std = TRUE, conf.int = TRUE)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}% use one of  RShowDoc("KEYWORDS")
\keyword{survival}% __ONLY ONE__ keyword per line
\keyword{RMST}% __ONLY ONE__ keyword per line
