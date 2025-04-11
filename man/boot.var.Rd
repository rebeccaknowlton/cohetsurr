\name{boot.var}
\alias{boot.var}
\title{
Performs bootstrap estimation procedures for the variance of the proportion of treatment effect explained, the omnibus test, and identifying a region above a treshold.
}
\description{
Performs bootstrap estimation procedures for the variance of the proportion of treatment effect explained, the omnibus test, and identifying a region above a treshold in a randomized treatment setting.
}
\usage{
boot.var(data.control, data.treat, W.grid.expand, type, test = FALSE, 
data.all = NULL, num.cov = NULL, results.for.test = NULL, threshold = NULL)
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
  \item{type}{
options are "model", "two step", or "both"; specifies the estimation method that should be used for the proportion of treatment effect explained
}
  \item{test}{
TRUE or FALSE, if test for heterogeneity is wanted
}
  \item{data.all}{
dataframe containing data from the control and treamtent groups, specifically the outcome Y, the surrogate S, and the baseline covariates W
}
  \item{num.cov}{
number of baseline covariates in the matrix W
}
  \item{results.for.test}{
the grid of estimates for the proportion of treatment effect explained obtained prior to bootstrapping, needed for the omnibus test
}
  \item{threshold}{
threshold to flag regions where the estimated proportion of the treatment effect explained is at least that high
}
}
\value{
A list is returned:
\item{return.grid}{grid of variance estimates for the overall treatment effect, the residual treatment effect, and the proportion of treatment effect explained as a function of the baseline covariates, W. If requested by user, includes regions flagged above the threshold.}
\item{pval}{p-value(s) from the F test and the two step omnibus test for heterogeneity, depending on type argument.}
}
\keyword{internal}