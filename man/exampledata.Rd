\name{exampledata}
\alias{exampledata}
\docType{data}
\title{
Example data
}
\description{
Example data
}
\usage{data("exampledata")}
\format{
A list with 7 elements representing 1000 observations from a treatment group and 1000 observations from a control group, and a grid of baseline covariate values at which to calculate estimates:
\describe{
  \item{\code{y}}{the outcome}
  \item{\code{s}}{the surrogate marker}
  \item{\code{a}}{the treatment assignment, where 1 indicates treatment and 0 indicates control}
  \item{\code{w1}}{the first baseline covariate of interest}
  \item{\code{w2}}{the second baseline covariate of interest}
  \item{\code{w1.grid}}{the grid of first baseline covariate values to provide estimates for}
  \item{\code{w2.grid}}{the grid of second baseline covariate values to provide estimates for}
}
}
\examples{
data(exampledata)
names(exampledata)
}

