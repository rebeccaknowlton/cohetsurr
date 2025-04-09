obs.boot.var <- function(df.train, df.test, type, numeric_predictors, categorical_predictors, threshold, use.actual.control.S, gam.smoothers = NULL, tree.tuners =NULL) {
  num.boot <- 200
  
  if (type %in% c("linear", "gam", "trees")) {
    boot.delta <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))
    boot.delta.s <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))
    boot.R.s <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))  
  } else {
    # dataframes for all types
    boot.delta.linear <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))
    boot.delta.s.linear <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))
    boot.R.s.linear <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))  
    boot.delta.gam <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))
    boot.delta.s.gam <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))
    boot.R.s.gam <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))  
    boot.delta.trees <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))
    boot.delta.s.trees <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))
    boot.R.s.trees <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))  
  }
  
  for (j in 1:num.boot) {
    boot.data.train <- df.train[sample(1:nrow(df.train), nrow(df.train), replace = TRUE),]
    
    boot.results <- obs.estimate.PTE(df.train = boot.data.train, df.test = df.test, type = type, numeric_predictors, categorical_predictors, use.actual.control.S, gam.smoothers = gam.smoothers, tree.tuners = tree.tuners)$df.test
    if (type %in% c("linear", "gam", "trees")) {
      boot.delta[,j] <- boot.results$delta
      boot.delta.s[,j] <- boot.results$delta.s
      boot.R.s[,j] <- boot.results$R.s
    } else {
      boot.delta.linear[,j] <- boot.results$delta.linear
      boot.delta.s.linear[,j] <- boot.results$delta.s.linear
      boot.R.s.linear[,j] <- boot.results$R.s.linear
      boot.delta.gam[,j] <- boot.results$delta.gam
      boot.delta.s.gam[,j] <- boot.results$delta.s.gam
      boot.R.s.gam[,j] <- boot.results$R.s.gam
      boot.delta.trees[,j] <- boot.results$delta.trees
      boot.delta.s.trees[,j] <- boot.results$delta.s.trees
      boot.R.s.trees[,j] <- boot.results$R.s.trees
    }
  }
  
  if (type %in% c("linear", "gam", "trees")) {
    df.test$delta.var <- apply(as.matrix(boot.delta), 1, mad)^2
    df.test$delta.s.var <- apply(as.matrix(boot.delta.s), 1, mad)^2
    df.test$R.s.var <- apply(as.matrix(boot.R.s), 1, mad)^2
    df.test$delta.lower <- matrixStats::rowQuantiles(as.matrix(boot.delta), probs = 0.025)
    df.test$delta.upper <- matrixStats::rowQuantiles(as.matrix(boot.delta), probs = 0.975)
    df.test$delta.s.lower <- matrixStats::rowQuantiles(as.matrix(boot.delta.s), probs = 0.025)
    df.test$delta.s.upper <- matrixStats::rowQuantiles(as.matrix(boot.delta.s), probs = 0.975)
    df.test$R.s.lower <- matrixStats::rowQuantiles(as.matrix(boot.R.s), probs = 0.025)
    df.test$R.s.upper <- matrixStats::rowQuantiles(as.matrix(boot.R.s), probs = 0.975)
  } else {
    df.test$delta.var.linear <- apply(as.matrix(boot.delta.linear), 1, mad)^2
    df.test$delta.s.var.linear <- apply(as.matrix(boot.delta.s.linear), 1, mad)^2
    df.test$R.s.var.linear <- apply(as.matrix(boot.R.s.linear), 1, mad)^2
    df.test$delta.lower.linear <- matrixStats::rowQuantiles(as.matrix(boot.delta.linear), probs = 0.025)
    df.test$delta.upper.linear <- matrixStats::rowQuantiles(as.matrix(boot.delta.linear), probs = 0.975)
    df.test$delta.s.lower.linear <- matrixStats::rowQuantiles(as.matrix(boot.delta.s.linear), probs = 0.025)
    df.test$delta.s.upper.linear <- matrixStats::rowQuantiles(as.matrix(boot.delta.s.linear), probs = 0.975)
    df.test$R.s.lower.linear <- matrixStats::rowQuantiles(as.matrix(boot.R.s.linear), probs = 0.025)
    df.test$R.s.upper.linear <- matrixStats::rowQuantiles(as.matrix(boot.R.s.linear), probs = 0.975)
    df.test$delta.var.gam <- apply(as.matrix(boot.delta.gam), 1, mad)^2
    df.test$delta.s.var.gam <- apply(as.matrix(boot.delta.s.gam), 1, mad)^2
    df.test$R.s.var.gam <- apply(as.matrix(boot.R.s.gam), 1, mad)^2
    df.test$delta.lower.gam <- matrixStats::rowQuantiles(as.matrix(boot.delta.gam), probs = 0.025)
    df.test$delta.upper.gam <- matrixStats::rowQuantiles(as.matrix(boot.delta.gam), probs = 0.975)
    df.test$delta.s.lower.gam <- matrixStats::rowQuantiles(as.matrix(boot.delta.s.gam), probs = 0.025)
    df.test$delta.s.upper.gam <- matrixStats::rowQuantiles(as.matrix(boot.delta.s.gam), probs = 0.975)
    df.test$R.s.lower.gam <- matrixStats::rowQuantiles(as.matrix(boot.R.s.gam), probs = 0.025)
    df.test$R.s.upper.gam <- matrixStats::rowQuantiles(as.matrix(boot.R.s.gam), probs = 0.975)
    df.test$delta.var.trees <- apply(as.matrix(boot.delta.trees), 1, mad)^2
    df.test$delta.s.var.trees <- apply(as.matrix(boot.delta.s.trees), 1, mad)^2
    df.test$R.s.var.trees <- apply(as.matrix(boot.R.s.trees), 1, mad)^2
    df.test$delta.lower.trees <- matrixStats::rowQuantiles(as.matrix(boot.delta.trees), probs = 0.025)
    df.test$delta.upper.trees <- matrixStats::rowQuantiles(as.matrix(boot.delta.trees), probs = 0.975)
    df.test$delta.s.lower.trees <- matrixStats::rowQuantiles(as.matrix(boot.delta.s.trees), probs = 0.025)
    df.test$delta.s.upper.trees <- matrixStats::rowQuantiles(as.matrix(boot.delta.s.trees), probs = 0.975)
    df.test$R.s.lower.trees <- matrixStats::rowQuantiles(as.matrix(boot.R.s.trees), probs = 0.025)
    df.test$R.s.upper.trees <- matrixStats::rowQuantiles(as.matrix(boot.R.s.trees), probs = 0.975)
  }
  
  if (!is.null(threshold)) {
    if (type %in% c("linear", "gam", "trees")) {
      df.test$p.val <- p.adjust(rowMeans(boot.R.s <= threshold), method = "BH")
    } else {
      df.test$p.val.linear <- p.adjust(rowMeans(boot.R.s.linear <= threshold), method = "BH")
      df.test$p.val.gam <- p.adjust(rowMeans(boot.R.s.gam <= threshold), method = "BH")
      df.test$p.val.trees <- p.adjust(rowMeans(boot.R.s.trees <= threshold), method = "BH")
    }
  }
  
  return(df.test)
}
