\name{delta.estimate}
\alias{delta.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates the treatment effect at time t, defined as the difference in the restricted mean survival time
}
\description{
Estimates the treatment effect at time t, defined as the difference in the restricted mean survival time.}
\usage{
delta.estimate(xone, xzero, deltaone, deltazero, t, std = FALSE, conf.int = FALSE, 
weight.perturb = NULL)
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
 \item{std}{
TRUE or FALSE; indicates whether standard error estimates should be provided, default is FALSE. Estimates are calculated using perturbation-resampling.  Two versions are provided: one that takes the standard deviation of the perturbed estimates (denoted as "sd") and one that takes the median absolute deviation (denoted as "mad").
}
  \item{conf.int}{
TRUE or FALSE; indicates whether 95\% confidence intervals should be provided. Confidence intervals are calculated using the percentiles of perturbed estimates, default is FALSE. If this is TRUE, standard error estimates are automatically provided.
}
  \item{weight.perturb}{
weights used for perturbation resampling.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
Let \eqn{G \in \{1,0\}} be the randomized treatment indicator and \eqn{T} denote the time of the primary outcome of interest. We use potential outcomes notation such that \eqn{T^{(G)}} denotes the time of the primary outcome under treatment G, for \eqn{G \in \{1, 0\}}. We define the treatment effect as the difference in restricted mean survival time up to a fixed time \eqn{t} under treatment 1 versus under treatment 0,
\deqn{ \Delta(t)=E\{T^{(1)}\wedge t\} - E\{T^{(0)}\wedge t \}}
  where \eqn{\wedge} indicates the minimum. Due to censoring, data consist of \eqn{n = n_1 + n_0} independent observations \eqn{\{X_{gi}, \delta_{gi}, i=1,...,n_g, g = 1,0\}}, where \eqn{X_{gi} = T_{gi}\wedge C_{ gi}}, \eqn{\delta_{gi} = I(T_{gi} < C_{gi})}, \eqn{C_{gi}} denotes the censoring time, \eqn{T_{gi}} denotes the time of the primary outcome, and \eqn{\{(T_{gi}, C_{gi}), i = 1, ..., n_g\}} are identically distributed within treatment group. The quantity \eqn{ \Delta(t)} is estimated using inverse probability of censoring weights:
\deqn{\hat{\Delta}(t)  = n_1^{-1} \sum_{i=1}^{n_1} \hat{M}_{1i}(t)- n_0^{-1} \sum_{i=1}^{n_0} \hat{M}_{0i}(t)}
where  \eqn{\hat{M}_{gi}(t) = I(X_{gi} > t)t/\hat{W}^C_g(t) +  I(X_{gi} < t)X_{gi}\delta_{gi}/\hat{W}^C_g(X_{gi})} and \eqn{\hat{W}^C_g(t)} is the Kaplan-Meier estimator of \eqn{P(C_{gi} \ge t).}

}
\value{
A list is returned:                
 \item{delta }{the estimate, \eqn{\hat{\Delta}(t)}, described above.}
  \item{rmst.1 }{the estimated restricted mean survival time in group 1, described above.}
   \item{rmst.0 }{the estimated restricted mean survival time in group 0, described above.}
 \item{delta.sd }{the standard error estimate of \eqn{\hat{\Delta}(t)}; if std = TRUE or conf.int = TRUE.}
 \item{delta.mad }{the standard error estimate of \eqn{\hat{\Delta}(t)} using the median absolute deviation; if std = TRUE or conf.int = TRUE.} 
    \item{conf.int.delta}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}(t)} based on sample quantiles of the perturbed values; if conf.int = TRUE.}
}
\references{
Parast L, Tian L, and Cai T (2020). Assessing the Value of a Censored Surrogate Outcome. Lifetime Data Analysis, 26(2):245-265.

Tian, L, Zhao, L, & Wei, LJ (2013). Predicting the restricted mean event time with the subject's baseline covariates in survival analysis. Biostatistics, 15(2), 222-233.

Royston, P, & Parmar, MK (2011). The use of restricted mean survival time to estimate the treatment effect in randomized clinical trials when the proportional hazards assumption is in doubt. Statistics in Medicine, 30(19), 2409-2421.
}
\author{
Layla Parast
}
\examples{
data(ExampleData)
names(ExampleData)

delta.estimate(xone = ExampleData$x1, xzero = ExampleData$x0, deltaone = ExampleData$delta1, 
deltazero = ExampleData$delta0, t = 5)

\donttest{
delta.estimate(xone = ExampleData$x1, xzero = ExampleData$x0, deltaone = ExampleData$delta1, 
deltazero = ExampleData$delta0, t = 5, std = TRUE, conf.int = TRUE)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}% use one of  RShowDoc("KEYWORDS")
\keyword{survival}% __ONLY ONE__ keyword per line
\keyword{RMST}% __ONLY ONE__ keyword per line
