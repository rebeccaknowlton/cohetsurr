\name{complex.heterogeneity}
\alias{complex.heterogeneity}
\title{
Estimates the proportion of treatment effect explained by the surrogate marker as a function of multiple baseline covariates.
}
\description{
Assesses complex heterogeneity in the utility of a surrogate marker by estimating the proportion of treatment effect explained by the surrogate marker as a function of multiple baseline covariates. Optionally, tests for evidence of heterogeneity overall and flags regions where the proportion of treatment effect explained is above a given threshold.
}
\usage{
complex.heterogeneity(y, s, a, W.mat, type = "model", variance = FALSE, 
test = FALSE, W.grid = NULL, grid.size = 4, threshold = NULL)
}
\arguments{
  \item{y}{
y, the outcome
}
  \item{s}{
s, the surrogate marker
}
  \item{a}{
a, the treatment assignment with 1 indicating the treatment group and 0 indicating the control group
}
  \item{W.mat}{
matrix of baseline covariate observations, where the first column is W1, second columns is W2, and so on.
}
  \item{type}{
options are "model", "two step", or "both"; specifies the estimation method that should be used for the proportion of treatment effect explained
}
  \item{variance}{
TRUE or FALSE, if variance/standard error estimates are wanted
}
  \item{test}{
TRUE or FALSE, if test for heterogeneity is wanted
}
  \item{W.grid}{
grid for the baseline covariates W where estimation will be provided
}
  \item{grid.size}{
number of measures for each baseline covariate to include in the estimation grid, if one is not provided by the user directly
}
  \item{threshold}{
threshold to flag regions where the estimated proportion of the treatment effect explained is at least that high
}
}
\value{
A list is returned:
\item{return.grid}{grid of estimates for the overall treatment effect, the residual treatment effect, and the proportion of treatment effect explained as a function of the baseline covariates, W. Includes variance estimates and regions flagged above the threshold, if specified by the user.}
\item{pval}{p-value(s) from the F test and the two step omnibus test for heterogeneity, depending on type argument}
}
\references{
Knowlton, R., Tian, L., Parast, L. (2024). "A General Framework to Assess Complex Heterogeneity in the Utility of a Surrogate Marker." Under Review.
}
\author{
Rebecca Knowlton
}
\examples{
  data(exampledata)
  names(exampledata)
  complex.heterogeneity(y = exampledata$y,
                        s = exampledata$s,
                        a = exampledata$a,
                        W.mat = matrix(cbind(exampledata$w1, exampledata$w2), ncol = 2),
                        type = "model",
                        W.grid = matrix(cbind(exampledata$w1.grid, exampledata$w2.grid),ncol=2))
  #computationally intensive
  \donttest{
    complex.heterogeneity(y = exampledata$y,
                          s = exampledata$s,
                          a = exampledata$a,
                          W.mat = matrix(cbind(exampledata$w1, 
                          exampledata$w2), ncol = 2),
                          type = "both",
                          variance = TRUE,
                          test = TRUE,
                          W.grid = matrix(cbind(exampledata$w1.grid, 
                          exampledata$w2.grid), ncol = 2),
                          threshold = 0.75)
  }
}
