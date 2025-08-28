\name{ExampleData}
\alias{ExampleData}
\docType{data}
\title{
Hypothetical data
}
\description{
Hypothetical data to be used in examples; t=5 and the landmark time = 2.

}
\usage{data(ExampleData)}
\format{
 A list with 6 elements representing 1000 observations from a control group and 1000 observations from a treatment group:
  \describe{
    \item{\code{s1}}{Time of the occurrence of the surrogate outcome for treated observations.}
    \item{\code{x1}}{The observed event or censoring time for treated observations; X  = min(T, C) where T is the time of the primary outcome and C is the censoring time.}
    \item{\code{delta1}}{The indicator identifying whether the treated observation was observed to have the event or was censored; D =1*(T<C) where T is the time of the primary outcome and C is the censoring time.}
   \item{\code{s0}}{Time of the occurrence of the surrogate outcome for control observations.}
    \item{\code{x0}}{The observed event or censoring time for control observations; X  = min(T, C) where T is the time of the primary outcome and C is the censoring time.}
    \item{\code{delta0}}{The indicator identifying whether the control observation was observed to have the event or was censored; D =1*(T<C) where T is the time of the primary outcome and C is the censoring time.}
}}
\details{
Note that the time of the surrogate outcome is used in all functions only if the surrogate outcome occurs before the minimum of the event time and censoring time.
}
\examples{
data(ExampleData)
names(ExampleData)
}
\keyword{datasets}
