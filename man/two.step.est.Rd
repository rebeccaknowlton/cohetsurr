\name{two.step.est}
\alias{two.step.est}
\title{
Estimates the proportion of treatment effect explained as a function of multiple baseline covariates, W, using a two step, semiparametric model. 
}
\description{
Estimates the proportion of treatment effect explained as a function of multiple baseline covariates, W, using a two step, semiparametric model in a randomized treatment setting.
}
\usage{
two.step.est(data.control, data.treat, W.grid.expand.function)
}
\arguments{
  \item{data.control}{
dataframe containing data from the control group, specifically the outcome Y, the surrogate S, and the baseline covariates W
}
  \item{data.treat}{
dataframe containing data from the treamtent group, specifically the outcome Y, the surrogate S, and the baseline covariates W
}
  \item{W.grid.expand.function}{
expanded version of the W grid of baseline covariates, where each row is a specific combination of the covariates for which the estimates should be provided
}
}
\value{
A grid of estimates is returned of the proportion of treatment effect explained, the overall treatment effect, and the residual treatment effect for the given baseline covariate combinations.
}
\keyword{internal}