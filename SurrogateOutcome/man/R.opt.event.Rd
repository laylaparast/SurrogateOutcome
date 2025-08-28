\name{R.opt.event}
\alias{R.opt.event}
\title{
Calculates the proportion of the treatment effect (the difference in survival at time t) explained by surrogate outcome information observed up to the landmark time, using the optimal transformation of the surrogate}
\description{
Calculates the proportion of the treatment effect (the difference in survival at time t) explained by surrogate outcome information observed up to the landmark time, using the optimal transformation of the surrogate; also provides standard error estimate and confidence interval. Details are provided in: Wang X, Cai T, Tian L, Parast L (2025). Model-free Approach to Evaluate a Censored Intermediate Outcome as a Surrogate for Overall Survival. arXiv preprint arXiv:2412.14129.
}
\usage{
R.opt.event(xone, xzero, deltaone, deltazero, sone, szero, t, landmark, std = FALSE, 
conf.int = FALSE, gopt = FALSE, ind = FALSE)
}
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
  \item{std}{
TRUE or FALSE; indicates whether standard error estimates should be provided, default is FALSE. Estimates are calculated using perturbation-resampling. 
}
  \item{conf.int}{
TRUE or FALSE; indicates whether a 95\% confidence interval for the PTE should be provided. Confidence intervals are calculated using the percentiles of perturbed estimates, default is FALSE. 
}
  \item{gopt}{
TRUE or FALSE; indicates whether the estimates of the optimal g1(s) and g2 should be provided, default is FALSE.
}
  \item{ind}{
TRUE or FALSE; indicates whether the estimate of PTE for primary outcome information only up to the landmark time should be provided, default is FALSE.
}
}
\value{
A list is returned:                
 \item{R.opt }{the estimate of the PTE of the surrogate outcome}
 \item{R.opt.ind }{the estimate of the PTE for primary outcome information only up to the landmark time; if ind=TRUE}
\item{g1.opt}{the estimate of the optimal transformation g1(s); if gopt=TRUE}
\item{g2.opt}{the estimate of the optimal transformation g2; if gopt=TRUE}
\item{R.opt.std }{the standard error estimate of the PTE of the surrogate outcome; if std=TRUE}
\item{conf.int.R }{the confidence interval for the PTE of the surrogate outcome; if cont.int=TRUE}
\item{g1.opt.sd}{the standard error estimate of the optimal transformation g1(s); if gopt=TRUE and std=TRUE}
\item{g2.opt.sd}{the standard error estimate of the optimal transformation g2; if gopt=TRUE and std=TRUE}
}
\references{
Wang X, Cai T, Tian L, Parast L (2025). Model-free Approach to Evaluate a Censored Intermediate Outcome as a Surrogate for Overall Survival. arXiv preprint arXiv:2412.14129.}
\author{
Xuan Wang
}
\examples{
data(ExampleData)
names(ExampleData)
\donttest{
R.opt.event(xone = ExampleData$x1, xzero = ExampleData$x0, deltaone = ExampleData$delta1, 
deltazero = ExampleData$delta0, sone = ExampleData$s1, szero = ExampleData$s0, t = 5, 
landmark=2)

#with all options as TRUE
R.opt.event(xone = ExampleData$x1, xzero = ExampleData$x0, deltaone = ExampleData$delta1, 
deltazero = ExampleData$delta0, sone = ExampleData$s1, szero = ExampleData$s0, t = 5, 
landmark=2, std = TRUE, conf.int = TRUE, gopt = TRUE, ind = TRUE)
}
}