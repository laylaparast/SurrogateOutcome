\name{R.q.event}
\alias{R.q.event}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the proportion of the treatment effect (the difference in restriced mean survival time at time t) explained by surrogate outcome information observed up to the landmark time
}
\description{
Calculates the proportion of the treatment effect (the difference in restriced mean survival time at time t) explained by surrogate outcome information observed up to the landmark time; also provides standard error estimate and confidence interval.
}
\usage{
R.q.event(xone, xzero, deltaone, deltazero, sone, szero, t, landmark, number = 40, 
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
Let \eqn{G \in \{1,0\}} be the randomized treatment indicator, \eqn{T} denote the time of the primary outcome of interest, and \eqn{S} denote the time of the surrogate outcome. We use potential outcomes notation such that \eqn{T^{(G)}} and \eqn{S^{(G)}} denote the respective times of the primary and surrogate outcomes under treatment G, for \eqn{G \in \{1, 0\}}.  In the absence of censoring, we only observe \eqn{(T, S)=(T^{(1)}, S^{(1)})} or \eqn{(T^{(0)}, S^{(0)})} for each individual depending on whether \eqn{G=1} or \eqn{0}. Due to censoring, data consist of \eqn{n = n_1 + n_0} independent observations \eqn{\{X_{gi}, \delta_{gi}, I(S_{gi}< t_0)I(X_{gi} > t_0),  S_{gi}\wedge t_0 I(X_{gi} > t_0), i=1,...,n_g, g = 1,0\}}, where \eqn{X_{gi} = T_{gi}\wedge C_{ gi}}, \eqn{\delta_{gi} = I(T_{gi} < C_{gi})}, \eqn{C_{gi}} denotes the censoring time, \eqn{T_{gi}} denotes the time of the primary outcome, \eqn{S_{gi}} denotes the time of the surrogate outcome, \eqn{\{(T_{gi}, C_{gi}, S_{gi}), i = 1, ..., n_g\}} are identically distributed within treatment group, and \eqn{t_0} is the landmark time of interest.

 We define the treatment effect as the difference in restricted mean survival time up to a fixed time \eqn{t} under treatment 1 versus under treatment 0,
\deqn{ \Delta(t)=E\{T^{(1)}\wedge t\} - E\{T^{(0)}\wedge t \}}
  where \eqn{\wedge} indicates the minimum. To define the proportion of treatment effect explained by the surrogate outcome information, let \deqn{Q_{t_0} ^{(g)} =  (Q_{t_01}, Q_{t_02})'=\{S ^{(g)} \wedge t_0I(T ^{(g)} > t_0), T^{(g)} I(T^{(g)} \le t_0)\}',  g=1, 0}
and define the residual treatment effect after accounting for the treatment effect on the surrogate outcome information as:
\deqn{ \Delta_Q(t,t_0) = P ^{(0)}_{t_0,2}\int_0^{t_0} \phi_1(t|t_0,s)dF_0(s) + P^{(0)}_{t_0,3}\psi_1(t|t_0) - P(T ^{(0)}> t_0) \nu_0(t|t_0)}  where  \eqn{P^{(0)}_{t_0,2} = P(T^{(0)} > {t_0}, S ^{(0)} < t_0)} and \eqn{ P^{(0)}_{t_0,3} = P(T^{(0)} > {t_0}, S ^{(0)} > t_0)},  \eqn{ \psi_1(t \mid t_0) = E(T^{(1)}\wedge t \mid  T^{(1)}> t_0, S^{(1)} > t_0)}, 
\eqn{\phi_1(t\mid t_0,s) = E(T^{(1)}\wedge t \mid  T ^{(1)}> t_0, S ^{(1)} = s), \quad \nu_0(t|t_0) = E(T ^{(0)} \wedge t | T ^{(0)}> t_0)}, and  \eqn{F_0(\cdot\mid t_0)} is the cumulative distribution function of \eqn{S^{(0)}} conditional on \eqn{T ^{(0)}> t_0} and \eqn{S ^{(0)} < t_0}. Then, the proportion of treatment effect on the primary outcome that is explained by surrogate information up to \eqn{t_0}, \eqn{Q_{t_0}}, can be expressed as a contrast between  \eqn{\Delta(t)} and \eqn{\Delta_Q(t,t_0)}:
\deqn{R_Q(t,t_0) = \{\Delta(t) - \Delta_Q(t,t_0) \} / \Delta(t) = 1-  \Delta_Q(t,t_0) / \Delta(t).} 
 
 
 The quantity \eqn{ \Delta(t)} is estimated using inverse probability of censoring weights:
\deqn{\hat{\Delta}(t)  = n_1^{-1} \sum_{i=1}^{n_1} \hat{M}_{1i}(t)- n_0^{-1} \sum_{i=1}^{n_0} \hat{M}_{0i}(t)}
where  \eqn{\hat{M}_{gi}(t) = I(X_{gi} > t)t/\hat{W}^C_g(t) +  I(X_{gi} < t)X_{gi}\delta_{gi}/\hat{W}^C_g(X_{gi})} and \eqn{\hat{W}^C_g(t)} is the Kaplan-Meier estimator of \eqn{P(C_{gi} \ge t).} The residual treatment effect \eqn{ \Delta_Q(t,t_0)} can be estimated nonparametrically or semi-parametrically. For nonparametric estimation, \eqn{ \psi_{1}(t|t_0)} is estimated by \eqn{\hat{\psi}_{1}(t|t_0) = \sum_{i=1}^{n_1}\frac{  { \hat{W}^C_1(t_0)} I(S_{1i}>t_0, X_{1i} > t_0) }{ \sum_{i=1}^{n_1}I(S_{1i}>t_0, X_{1i} > t_0)} \hat{M}_{1i}(t)}, and \eqn{\phi_1(t \mid t_0,s) = E(T^{(1)}\wedge t\mid  X^{(1)}> t_0, S ^{(1)} = s)} is estimated using a
nonparametric kernel Nelson-Aalen estimator for \eqn{\Lambda_1(t\mid t_0,s ),} the cumulative hazard function of \eqn{T^{(1)}} conditional on \eqn{S^{(1)}=s} and \eqn{T^{(1)}>t_0,} as
\deqn{\hat \phi_1(t \mid t_0,s)  = t_0+\int_ {t_0}^t \exp\{-\hat{\Lambda}_1(t\mid t_0,s) \}dt,}
where
\deqn{\hat{\Lambda}_1(t\mid t_0,s) = \int_{t_0}^t \frac{\sum_{i=1}^{n_1} I(X_{1i}>t_0, {S_{1i} < t_0})  K_h\{\gamma(S_{1i}) - \gamma(s)\}dN_{1i}(z)}{\sum_{i=1}^{n_1}  I(X_{1i}>t_0,  {S_{1i} < t_0})  K_h\{\gamma(S_{1i}) - \gamma(s)\} Y_{1i}(z)},}
is a consistent estimate of \eqn{\Lambda_1(t\mid t_0,s ),} \eqn{Y_{1i}(t) = I(X_{1i} \geq t)},  \eqn{N_{1i}(t) = I(X_{1i} \leq t) \delta_i,  K(\cdot)} is a smooth symmetric density function, \eqn{K_h(x) = K(x/h)/h}, \eqn{\gamma(\cdot)} is a given monotone transformation function, and \eqn{h=O(n_1^{-\eta})} is a specified bandwidth with \eqn{\eta \in (1/2,1/4)}.  Finally, we let \deqn{
\hat{\nu}_{0}(t|t_0) = \sum_{i=1}^{n_0}\frac{  {\hat{W}^C_0(t_0)}I(X_{0i} > t_0) }{ \sum_{i=1}^{n_0}I(X_{0i} > t_0)} \hat{M}_{0i}(t).}
We then estimate  \eqn{\Delta_{Q}(t,t_0)} as \eqn{\hat{\Delta}_{Q}(t,t_0)} defined as
\deqn{
n_0^{-1} \sum_{i=1}^{n_0} \left \{ \frac{I_{t_0,2}(X_{0i}, S_{0i})\hat{\phi}_1(t\mid t_0, S_{0i}) + I_{t_0,3}(X_{0i}, S_{0i})\hat{\psi}_1(t\mid t_0) - I_{t_0}(X_{0i})\hat{\nu}(t|t_0) }{\hat{W}^C_0(t_0)} \right \}}
where \eqn{ I_{t_0,2}(x, s) = I(x > {t_0}, s < t_0)} and \eqn{ I_{t_0,3}(x, s) = I(x > {t_0}, s > t_0)} and \eqn{ I_{t_0}(x)=I(x > {t_0})} and thus, \eqn{\hat{R}_Q(t,t_0) =1- \hat{\Delta}_Q(t,t_0)/\hat{\Delta}(t).}

For the semi-parametric estimate,  \eqn{\hat \phi_1(t| t_0,s)} is replaced with an estimate obtained using a landmark Cox proportional hazards model \deqn{P(T^{(1)}> t\mid  T^{(1)}> t_0, S^{(1)} < t_0, S ^{(1)}) = \exp \{ -\Lambda_0(t|t_0)\exp(\beta_0S ^{(1)})\}}
 where \eqn{\Lambda_0(t|t_0)} is the unspecified baseline cumulative hazard among \eqn{\Omega_{t_0} = \{T^{(1)}> t_0, S^{(1)} < t_0\}} and \eqn{\beta_0} is unknown. That is, let \eqn{\tilde \phi_1(t| t_0,s) =  t_0+\int_{t_0}^{t}\exp \{ -\hat{\Lambda}_0(t|t_0)\exp(\hat{\beta}s)\} dt,}  where \eqn{\hat{\beta}} is estimated by fitting a Cox model to the subpopulation \eqn{\Omega_{t_0}} with a single predictor \eqn{S} and \eqn{\hat{\Lambda}_0(\cdot|t_0)} is the corresponding Breslow estimator. Then the semiparametric estimator for
\eqn{\Delta_{Q}(t,t_0)} is \eqn{\tilde{\Delta}_{Q}(t,t_0)} defined as
\deqn{ n_0^{-1} \sum_{i=1}^{n_0} \left \{ \frac{I_{t_0,2}(X_{0i}, S_{0i})\tilde{\phi}_1(t\mid t_0, S_{0i}) + I_{t_0,3}(X_{0i}, S_{0i})\hat{\psi}_1(t\mid t_0) - I_{t_0}(X_{0i})\hat{\nu}(t|t_0) }{\hat{W}^C_0(t_0)} \right \} }
 and \eqn{\tilde{R}_Q(t,t_0) =1- \tilde{\Delta}_Q(t,t_0)/\hat{\Delta}(t).}

}
\value{
A list is returned:                
 \item{delta }{the estimate, \eqn{\hat{\Delta}(t)}, described in delta.estimate documentation.}
 \item{delta.q }{the estimate, \eqn{\hat{\Delta}_Q(t,t_0)}, described above.}
 \item{R.q }{the estimate, \eqn{\hat{R}_Q(t,t_0)}, described above.}
 \item{delta.sd }{the standard error estimate of \eqn{\hat{\Delta}(t)}; if std = TRUE or conf.int = TRUE.}
 \item{delta.mad }{the standard error estimate of \eqn{\hat{\Delta}(t)} using the median absolute deviation; if std = TRUE or conf.int = TRUE.} 
  \item{delta.q.sd}{the standard error estimate of \eqn{\hat{\Delta}_Q(t,t_0)}; if std = TRUE or conf.int = TRUE.}
  \item{delta.q.mad}{the standard error estimate of \eqn{\hat{\Delta}_Q(t,t_0)} using the median absolute deviation; if std = TRUE or conf.int = TRUE.}  
 \item{R.q.sd }{the standard error estimate of \eqn{\hat{R}_Q(t,t_0)}; if std = TRUE or conf.int = TRUE.}
  \item{R.q.mad }{the standard error estimate of \eqn{\hat{R}_Q(t,t_0)} using the median absolute deviation; if std = TRUE or conf.int = TRUE.}
  \item{conf.int.delta}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}(t)} based on sample quantiles of the perturbed values; if conf.int = TRUE.}
 \item{conf.int.delta.q}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}_Q(t,t_0)} based on sample quantiles of the perturbed values; if conf.int = TRUE.} 
 \item{conf.int.R.q}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_Q(t,t_0)} based on sample quantiles of the perturbed values; if conf.int = TRUE.}
}
\references{
Parast L, Tian L, and Cai T (2020). Assessing the Value of a Censored Surrogate Outcome. Lifetime Data Analysis, 26(2):245-265.

Parast, L and Cai, T (2013). Landmark risk prediction of residual life for breast cancer survival. Statistics in Medicine, 32(20), 3459-3471.

}
\author{
Layla Parast
}
\examples{
data(ExampleData)
names(ExampleData)

\donttest{
R.q.event(xone = ExampleData$x1, xzero = ExampleData$x0, deltaone = ExampleData$delta1, 
deltazero = ExampleData$delta0, sone = ExampleData$s1, szero = ExampleData$s0, t = 5, 
landmark=2, type = "np")
R.q.event(xone = ExampleData$x1, xzero = ExampleData$x0, deltaone = ExampleData$delta1, 
deltazero = ExampleData$delta0, sone = ExampleData$s1, szero = ExampleData$s0, t = 5, 
landmark=2, type = "semi")
R.q.event(xone = ExampleData$x1, xzero = ExampleData$x0, deltaone = ExampleData$delta1, 
deltazero = ExampleData$delta0, sone = ExampleData$s1, szero = ExampleData$s0, t = 5, 
landmark=2, type = "np", std = TRUE, conf.int = TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}% use one of  RShowDoc("KEYWORDS")
\keyword{survival}% __ONLY ONE__ keyword per line
\keyword{RMST}% __ONLY ONE__ keyword per line
