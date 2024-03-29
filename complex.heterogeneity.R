#' complex.heterogeneity
#'
#' @param y true outcome of interest
#' @param s potential surrogate marker
#' @param a binary treatment assignment indicator
#' @param W matrix of covariates of interest
#' @param type either "model" or "two step"
#' @param variance TRUE if you want to bootstrap variance estimates
#' @param W.grid optional grid of W values to get point estimates for
#' @param grid.size optional size of the auto-created W.grid
#'
#' @return a grid containing the W.grid values, the point estimates of delta, delta.s, and R.s, and their associated variance estimates if requested
#' @export
#'
#' @examples
complex.heterogeneity <- function(y, s, a, W.mat, type = "model", variance = FALSE, test = FALSE, W.grid = NULL, grid.size = 4, threshold, h.0,h.1,h.3) {
  # create dataframe for control and treat
  W.mat.control <- W.mat[a==0,]
  W.mat.treat <- W.mat[a==1,]
  covariates.control <- split(W.mat.control, rep(1:ncol(W.mat.control), each = nrow(W.mat.control)))
  covariates.treat <- split(W.mat.treat, rep(1:ncol(W.mat.treat), each = nrow(W.mat.treat)))
  data.control <- cbind(data.frame(Y = y[a==0], S = s[a==0]), covariates.control)
  data.treat <- cbind(data.frame(Y = y[a==1], S = s[a==1]), covariates.treat)

  #create dataframe for combined data
  covariates.all <- split(W.mat, rep(1:ncol(W.mat), each = nrow(W.mat)))
  num.cov <- length(covariates.all)
  data.all <- cbind(data.frame(Y = y, S = s, A = a), covariates.all)
  for (i in 1:num.cov) {
    names(data.all)[3+i] <- paste0("W",i)
  }

  # create W.grid if user doesn't specify
  if (is.null(W.grid)) {
    W.grid <- seq(quantile(W.mat[,1], 0.2), quantile(W.mat[,1], 0.8), length = grid.size)
    for (i in 2:ncol(W.mat)) {
      W.grid <- cbind(W.grid, seq(quantile(W.mat[,i], 0.2), quantile(W.mat[,i], 0.8), length = grid.size))
    }
  }

  # expand W.grid
  W.grid.expand <- expand.grid(split(W.grid, rep(1:ncol(W.grid), each = nrow(W.grid))))

  if (type == "model") {
    return.grid.p <- parametric.est(data.control, data.treat, W.grid.expand)
    return.grid = return.grid.p
  }
  if (type == "two step") {
    return.grid.t <- two.step.est(data.control, data.treat, W.grid.expand,h.0=h.0,h.1=h.1,h.3=h.3)
    return.grid = return.grid.t
  }
  if (type == "both") {
    return.grid.p <- parametric.est(data.control, data.treat, W.grid.expand)
    return.grid.t <- two.step.est(data.control, data.treat, W.grid.expand,h.0=h.0,h.1=h.1,h.3=h.3)
    return.grid = cbind(return.grid.p,return.grid.t)
  }
  if (variance == TRUE | test == TRUE) {
  	boot.object = boot.var(data.control, data.treat, W.grid.expand, type,test=test, data.all = data.all, num.cov = num.cov, results.for.test = return.grid, threshold = threshold,h.0=h.0,h.1=h.1,h.3=h.3)
    return.grid <- cbind(return.grid, boot.object$my.grid)
    return.grid <- list(return.grid = return.grid, pval = boot.object$pval)
  }
  else {return.grid <- list(return.grid = return.grid)}
  return(return.grid)
}
