\name{obs.boot.var}
\alias{obs.boot.var}
\title{
Calculate bootstrapped variance estimates in an observational setting.
}
\description{
Calculates bootstrapped variance estimates of delta, delta.s, and R.s, and optionally calculates p-values for for identifying individual for whom the surrogate is strong.
}
\usage{
obs.boot.var(df.train, df.test, type, numeric_predictors, categorical_predictors, 
  threshold, use.actual.control.S, gam.smoothers, tree.tuners)
}
\arguments{
  \item{df.train}{
A dataframe containing training data.
}
  \item{df.test}{
A dataframe containing testing data.
}
  \item{type}{
Options are "linear", "gam", "trees", or "all"; type of base learners to use.
}
  \item{numeric_predictors}{
The column names in the dataframes that represent numeric baseline covariates.
}
  \item{categorical_predictors}{
The column names in the dataframes that represent categorical baseline covariates.
}
  \item{threshold}{
An optional threshold to test individuals for the null hypothesis that PTE > threshold.
}
  \item{use.actual.control.S}{
TRUE or FALSE, if user prefers to use the actual observed values for the surrogate in the control group instead of predicting values from the base learners.
}
  \item{gam.smoothers}{
A list of smoothing parameters to use for GAM base learners, so they are not retuned with bootstrapping iterations ("m1sp", "m0sp", "m1ssp", "m0ssp", "s0")
}
  \item{tree.tuners}{
A list of tuning parameters to use for tree base learners, so they are not retuned with bootstrapping iterations ("m1sp", "m0sp", "m1ssp", "m0ssp", "s0") 
}
}
\value{
A dataframe is returned, which is the df.test argument with new columns appended for the estimated variances of delta, delta.s, and R.s, as well as p-values if a threshold is provided.
}
\keyword{internal}