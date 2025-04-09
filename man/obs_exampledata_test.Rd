\name{obs_exampledata_test}
\alias{obs_exampledata_test}
\docType{data}
\title{
Example testing data for observational setting
}
\description{
Example testing data for observational setting
}
\usage{data("obs_exampledata_test")}
\format{
  A data frame with 200 observations on the following 9 variables.
  \describe{
    \item{\code{X1}}{a numeric baseline covariate of interest}
    \item{\code{X2}}{a numeric baseline covariate of interest}
    \item{\code{X3}}{a numeric baseline covariate of interest}
    \item{\code{X4}}{a numeric baseline covariate of interest}
    \item{\code{X5}}{a numeric baseline covariate of interest}
    \item{\code{X6}}{a numeric baseline covariate of interest}
    \item{\code{G}}{the treatment assignment, where 1 indicates treated and 0 indicates control}
    \item{\code{S}}{the surrogate marker}
    \item{\code{Y}}{the primary outcome}
  }
}
\examples{
data(obs_exampledata_test)
names(obs_exampledata_test)
}
