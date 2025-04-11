\name{obs.het.surr}
\alias{obs.het.surr}
\title{
Estimate the proportion of the treatment effect explained by the surrogate marker as a function of multiple baseline covariates in an observational setting.
}
\description{
Assesses surrogate heterogeneity in real world data by estimating the proportion of the treatment effect explained as a function of baseline covariates. Optionally tests individuals for strong surrogacy based on a threshold.
}
\usage{
obs.het.surr(df.train, df.test, type, var.want = FALSE, threshold = NULL, 
  use.actual.control.S = FALSE)
}
\arguments{
  \item{df.train}{
dataframe containing training data; must have columns G (treatment assignment), S (surrogate marker), and Y (primary outcome), in addition to the baseline covariates of interest
}
  \item{df.test}{
dataframe containing testing data; must contain the same baseline covariate columns as the training data
}
  \item{type}{
options are "linear", "gam", "trees", or "all"; type of base learners to use
}
  \item{var.want}{
TRUE or FALSE, if variance estimates are wanted
}
  \item{threshold}{
optional threshold to test individuals for the null hypothesis that PTE is greater than the threshold; must have var.want = TRUE to return p-values
}
  \item{use.actual.control.S}{
TRUE or FALSE, if user prefers to use the actual observed values for the surrogate in the control group instead of predicting values from the base learners
}
}
\value{
A dataframe is returned, which is the df.test argument with new columns appended for the estimates and corresponding variances of delta, delta.s, and R.s. If a threshold is specified, returns a p-value for the null hypothesis that PTE > threshold.
}
\references{
Knowlton, R. and Parast, L. (2025) ``Assessing Surrogate Heterogeneity in Real World Data Using Meta-Learners." Under Review.
}
\author{
Rebecca Knowlton
}

\examples{
data(obs_exampledata_train)
data(obs_exampledata_test)
obs.het.surr(df.train = obs_exampledata_train, df.test = obs_exampledata_test,
type = "linear", var.want = FALSE)
}
