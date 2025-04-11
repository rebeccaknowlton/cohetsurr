\name{parametric.est}
\alias{parametric.est}
\title{
Estimates the proportion of treatment effect explained as a function of multiple baseline covariates, W, using a parametric model. 
}
\description{
Estimates the proportion of treatment effect explained as a function of multiple baseline covariates, W, using a parametric model in a randomized treatment setting.
}
\usage{
parametric.est(data.control, data.treat, W.grid.expand)
}
\arguments{
  \item{data.control}{
dataframe containing data from the control group, specifically the outcome Y, the surrogate S, and the baseline covariates W
}
  \item{data.treat}{
dataframe containing data from the treamtent group, specifically the outcome Y, the surrogate S, and the baseline covariates W
}
  \item{W.grid.expand}{
expanded version of the W grid of baseline covariates, where each row is a specific combination of the covariates for which the estimates should be provided
}
}
\value{
A grid of estimates is returned of the proportion of treatment effect explained, the overall treatment effect, and the residual treatment effect for the given baseline covariate combinations.
}
\keyword{internal}