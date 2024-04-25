parametric.est <-
function(data.control, data.treat, W.grid.expand) {
   control.model <- lm(Y ~ ., data = data.control)
  treatment.model <- lm(Y ~ ., data = data.treat)

  alpha_1 <- mean(data.treat$S)
  alpha_0 <- mean(data.control$S)

  coef.diffs <- treatment.model$coef - control.model$coef

  beta_1 <- as.numeric(coef.diffs[1])
  beta_2 <- as.numeric(control.model$coef[2])
  beta_3 <- as.numeric(coef.diffs[2])
  beta_5 <- coef.diffs[3:length(coef.diffs)]

  my.grid <- W.grid.expand
  my.grid$delta <- beta_1 + (beta_2 + beta_3) * alpha_1 +
    as.matrix(W.grid.expand) %*% matrix(beta_5) - beta_2 * alpha_0

  my.grid$delta.s <- beta_1 + beta_3 * alpha_0 +
    as.matrix(W.grid.expand) %*% matrix(beta_5)

  my.grid$R.s <- 1 - my.grid$delta.s / my.grid$delta

  return(my.grid)
}
