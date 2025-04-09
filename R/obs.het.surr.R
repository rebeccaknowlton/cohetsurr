obs.het.surr <- function(df.train, df.test, type, var.want = FALSE, threshold = NULL, use.actual.control.S = FALSE) {
  # Check that inputs valid
  if (!is.data.frame(df.train) || !is.data.frame(df.test)) {
    stop("Both df.train and df.test must be data frames.")
  }
  if (!(type %in% c("linear", "gam", "trees", "all"))) {
    stop('The "type" argument must be one of "linear", "gam", "trees", or "all".')
  }
  
  required_columns <- c("G", "S", "Y")
  missing_columns_train <- setdiff(required_columns, colnames(df.train))
  
  if (length(missing_columns_train) > 0) {
    stop(paste("df.train is missing the following required columns:", paste(missing_columns_train, collapse = ", ")))
  }
  
  predictor_columns <- setdiff(colnames(df.train), required_columns)
  if (length(predictor_columns) == 0) {
    stop("No predictor columns (e.g., X1, X2, ...) found in the data frames.")
  }
  
  missing_predictors_test <- setdiff(predictor_columns, colnames(df.test))
  if (length(missing_predictors_test) > 0) {
    stop(paste(
      "df.test is missing the following predictor columns from df.train:",
      paste(missing_predictors_test, collapse = ", ")
    ))
  }
  
  # Separate predictors into numeric and categorical
  numeric_predictors <- predictor_columns[sapply(df.train[predictor_columns], is.numeric)]
  categorical_predictors <- predictor_columns[sapply(df.train[predictor_columns], function(col) is.factor(col) || is.character(col))]
  
  # calculate point estimates
  hold = obs.estimate.PTE(df.train = df.train, df.test = df.test, type = type, 
                          numeric_predictors = numeric_predictors, 
                          categorical_predictors = categorical_predictors,
                          use.actual.control.S = use.actual.control.S, want.smooth = TRUE, want.tune = TRUE)
  # calculate variance and flag region if desired
  df.test <- hold$df.test
  if (var.want == TRUE) {
    df.test <- obs.boot.var(df.train = df.train, df.test = df.test, type = type, 
                        numeric_predictors = numeric_predictors, 
                        categorical_predictors = categorical_predictors,
                        threshold = threshold, use.actual.control.S = use.actual.control.S, gam.smoothers = hold$smooth_params, tree.tuners = hold$tuner_params)
  }
  return(df.test)
}
