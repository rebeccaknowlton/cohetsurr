\name{obs.estimate.PTE}
\alias{obs.estimate.PTE}
\title{
Estimate the proportion of the treatment effect explained in an observational setting.
}
\description{
Fits base learners using the specified type of model on the training data, and uses those models to calculate delta, delta.s, and R.s on the testing set.
}
\usage{
obs.estimate.PTE(df.train, df.test, type, numeric_predictors, categorical_predictors, 
  use.actual.control.S, gam.smoothers, tree.tuners, want.smooth, want.tune)
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
  \item{use.actual.control.S}{
TRUE or FALSE, if user prefers to use the actual observed values for the surrogate in the control group instead of predicting values from the base learners.
}
  \item{gam.smoothers}{
A list of smoothing parameters to use for GAM base learners, so they are not retuned with bootstrapping iterations ("m1sp", "m0sp", "m1ssp", "m0ssp", "s0")
}
  \item{tree.tuners}{
A list of tuning parameters to use for tree base learners, so they are not retuned with bootstrapping iterations ("m1sp", "m0sp", "m1ssp", "m0ssp", "s0") 
}
  \item{want.smooth}{
TRUE or FALSE, if smoothing parameters for GAM should be saved  
}
  \item{want.tune}{
TRUE or FALSE, if tuning parameters for trees should be saved
}
}
\value{
A list is returned:
  \item{df.test}{df.test argument with new columns appended for the estimates of delta, delta.s, and R.s}
  \item{smooth_params}{A list of smoothing parameters used for GAM base learners ("m1sp", "m0sp", "m1ssp", "m0ssp", "s0")}
  \item{tuner_params}{A list of tuning parameters used for tree base learners ("m1sp", "m0sp", "m1ssp", "m0ssp", "s0")}
}
\keyword{internal}