complex.heterogeneity <-
function(y, s, a, W.mat, type = "model", variance = FALSE, test = FALSE, W.grid = NULL, grid.size = 4, threshold = NULL) {
  if(ncol(W.mat) < 2) {stop("If only using one baseline covariate, see the package hetsurr instead")}
  w.numeric.check <- 0
  w.idx.drop <- c()
  for (w in 1:ncol(W.mat)) {
    if(inherits(W.mat[1,w], "character") | inherits(W.mat[1,w], "factor")) {w.numeric.check <- 1}
    if(is.null(W.grid) & all(W.mat[,w] %in% 0:1)) {w.idx.drop <- c(w.idx.drop, w)}
  }
  if(w.numeric.check == 1) {stop("W.mat should only contain numeric or binary variables. For categorical variables, user must input data as binary indicators.")}
  
  # create dataframe for control and treat
  W.mat.control <- W.mat[a==0,]
  W.mat.treat <- W.mat[a==1,]
  covariates.control <- split(W.mat.control, rep(1:ncol(W.mat.control), each = nrow(W.mat.control)))
  covariates.treat <- split(W.mat.treat, rep(1:ncol(W.mat.treat), each = nrow(W.mat.treat)))
  data.control <- cbind(data.frame(Y = y[a==0], S = s[a==0]), covariates.control)
  data.treat <- cbind(data.frame(Y = y[a==1], S = s[a==1]), covariates.treat)
  if(ncol(W.mat) > min(nrow(data.control), nrow(data.treat))) {warning("The dimension of W is greater than the sample size; method may not perform well")}
  if((min(nrow(data.control), nrow(data.treat)) < 200) & (type != "model")) {warning("Small sample size; kernel smoothing may not perform well")}
  
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
  W.grid.expand <- na.omit(W.grid.expand)
  for (w in w.idx.drop) {
    W.grid.expand <- W.grid.expand[which((round(W.grid.expand[,w],2) == 0.00) | (round(W.grid.expand[,w],2) == 1.00)),]
  }
  rownames(W.grid.expand) <- NULL

  if (type == "model") {
    return.grid.p <- parametric.est(data.control, data.treat, W.grid.expand)
    return.grid = return.grid.p
  }
  if (type == "two step") {
    return.grid.t <- two.step.est(data.control, data.treat, W.grid.expand)
    if(return.grid.t$extrapol.warn == 1) {warning("Extrapolation used in kernel smoothing")}
    return.grid = return.grid.t$my.grid
  }
  if (type == "both") {
    return.grid.p <- parametric.est(data.control, data.treat, W.grid.expand)
    return.grid.t <- two.step.est(data.control, data.treat, W.grid.expand)
    if(return.grid.t$extrapol.warn == 1) {warning("Extrapolation used in kernel smoothing")}
    return.grid = cbind(return.grid.p,return.grid.t$my.grid)
  }
  if (variance == TRUE | test == TRUE) {
  	boot.object = boot.var(data.control, data.treat, W.grid.expand, type,test=test, data.all = data.all, num.cov = num.cov, results.for.test = return.grid, threshold = threshold)
    return.grid <- cbind(return.grid, boot.object$my.grid)
    if(!is.null(colnames(W.mat))) {colnames(return.grid)[1:num.cov] <- colnames(W.mat)}
    else if (!is.null(colnames(W.grid))) {colnames(return.grid)[1:num.cov] <- colnames(W.grid)}
    return.grid <- list(return.grid = return.grid, pval = boot.object$pval)
  }
  else {
    if(!is.null(colnames(W.mat))) {colnames(return.grid)[1:num.cov] <- colnames(W.mat)}
    else if (!is.null(colnames(W.grid))) {colnames(return.grid)[1:num.cov] <- colnames(W.grid)}
    return.grid <- list(return.grid = return.grid)
  }
  return(return.grid)
}
