boot.var <-
function(data.control, data.treat, W.grid.expand, type, test=FALSE, data.all = NULL, num.cov=NULL, results.for.test = NULL, threshold = NULL) {
  num.boot <- 200
  boot.delta.est <- data.frame(matrix(nrow = nrow(W.grid.expand), ncol = num.boot))
  boot.delta.s.est <- data.frame(matrix(nrow = nrow(W.grid.expand), ncol = num.boot))
  boot.R.s.est <- data.frame(matrix(nrow = nrow(W.grid.expand), ncol = num.boot))
	
  boot.delta.est.two.step <- data.frame(matrix(nrow = nrow(W.grid.expand), ncol = num.boot))
  boot.delta.s.est.two.step <- data.frame(matrix(nrow = nrow(W.grid.expand), ncol = num.boot))
  boot.R.s.est.two.step <- data.frame(matrix(nrow = nrow(W.grid.expand), ncol = num.boot))

  boot.R.d.est.two.step <- data.frame(matrix(nrow = nrow(W.grid.expand), ncol = num.boot))

  F.test <-  function(data.all, num.cov) {
    model1.form <- "Y ~ . + A*S"
    for (i in 1:num.cov) {
      model1.form <- paste0(model1.form, " + A*W", i)
    }
    model1 <- lm(model1.form, data.all)
    model2 <- lm(Y ~ . + A*S, data.all)
    
    pval = anova(model1, model2)[2,6]
    
    return(pval)
  }
  
  for (j in 1:num.boot) {
  		boot.data.control <- data.control[sample(1:nrow(data.control), nrow(data.control), replace = TRUE),]
      boot.data.treat <- data.treat[sample(1:nrow(data.treat), nrow(data.treat), replace = TRUE),]
      if (type == "model") {
        boot.estimate <- parametric.est(boot.data.control, boot.data.treat, W.grid.expand)
        boot.delta.est[,j] <- boot.estimate$delta
      	boot.delta.s.est[,j] <- boot.estimate$delta.s
      	boot.R.s.est[,j] <- boot.estimate$R.s
         }
      if (type == "two step") {
        boot.estimate.two.step <- two.step.est(boot.data.control, boot.data.treat, W.grid.expand)$my.grid
        boot.delta.est.two.step[,j] <- boot.estimate.two.step$delta.two.step
      	boot.delta.s.est.two.step[,j] <- boot.estimate.two.step$delta.s.two.step
      	boot.R.s.est.two.step[,j] <- boot.estimate.two.step$R.s.two.step
         }
      if (type == "both"){
      	boot.estimate <- parametric.est(boot.data.control, boot.data.treat, W.grid.expand)
      	boot.delta.est[,j] <- boot.estimate$delta
      	boot.delta.s.est[,j] <- boot.estimate$delta.s
      	boot.R.s.est[,j] <- boot.estimate$R.s

      	boot.estimate.two.step <- two.step.est(boot.data.control, boot.data.treat, W.grid.expand)$my.grid
      	  	boot.delta.est.two.step[,j] <- boot.estimate.two.step$delta.two.step
      	boot.delta.s.est.two.step[,j] <- boot.estimate.two.step$delta.s.two.step
      	boot.R.s.est.two.step[,j] <- boot.estimate.two.step$R.s.two.step

      }
        if(test) {
        	if(type == "two step" | type == "both"){
   	  tau.boot.two.step = mean(boot.estimate.two.step$R.s.two.step)
	  boot.R.d.est.two.step[,j]  = boot.estimate.two.step$R.s.two.step - tau.boot.two.step
          }
        }
  }
  
  my.grid = c()
  if(type == "model" | type == "both"){
      my.grid <- cbind(apply(as.matrix(boot.delta.est), 1, mad)^2,
                       apply(as.matrix(boot.delta.s.est), 1, mad)^2,
                       apply(as.matrix(boot.R.s.est), 1, mad)^2, 
                       matrixStats::rowQuantiles(as.matrix(boot.delta.est), probs=0.025),
                       matrixStats::rowQuantiles(as.matrix(boot.delta.est), probs=0.975),
                       matrixStats::rowQuantiles(as.matrix(boot.delta.s.est), probs=0.025),
                       matrixStats::rowQuantiles(as.matrix(boot.delta.s.est), probs=0.975),									   
                       matrixStats::rowQuantiles(as.matrix(boot.R.s.est), probs=0.025),
                       matrixStats::rowQuantiles(as.matrix(boot.R.s.est), probs=0.975),
                       apply(as.matrix(boot.R.s.est), 1, function(x) mean(x<=threshold)))
      my.grid[,10] <- as.numeric(p.adjust(my.grid[,10], method = "BH") < 0.05)
  		colnames(my.grid) <- c("delta.var", "delta.s.var", "R.s.var", "delta.lower", "delta.upper", "delta.s.lower", "delta.s.upper", "R.s.lower", "R.s.upper","threshold.flag")
  		if(is.null(threshold)) {my.grid <- my.grid[,-which(colnames(my.grid) == "threshold.flag")]}
  }

	if(type == "two step" | type == "both") {
			my.grid.two.step <- cbind(apply(as.matrix(boot.delta.est.two.step), 1, mad)^2,
                                apply(as.matrix(boot.delta.s.est.two.step), 1, mad)^2,
                                apply(as.matrix(boot.R.s.est.two.step), 1, mad)^2,
                                matrixStats::rowQuantiles(as.matrix(boot.delta.est.two.step), probs=0.025),
                                matrixStats::rowQuantiles(as.matrix(boot.delta.est.two.step), probs=0.975),
                                matrixStats::rowQuantiles(as.matrix(boot.delta.s.est.two.step), probs=0.025),
                                matrixStats::rowQuantiles(as.matrix(boot.delta.s.est.two.step), probs=0.975),
                                matrixStats::rowQuantiles(as.matrix(boot.R.s.est.two.step), probs=0.025),
                                matrixStats::rowQuantiles(as.matrix(boot.R.s.est.two.step), probs=0.975),
                                apply(as.matrix(boot.R.s.est.two.step), 1, function(x) mean(x<threshold)))
			my.grid.two.step[,10] <- as.numeric(p.adjust(my.grid.two.step[,10], method = "BH") < 0.05)
  		colnames(my.grid.two.step) <- c("delta.var.two.step", "delta.s.var.two.step", "R.s.var.two.step","delta.lower.two.step", "delta.upper.two.step","delta.s.lower.two.step", "delta.s.upper.two.step","R.s.lower.two.step", "R.s.upper.two.step","threshold.flag.two.step")
  		if(is.null(threshold)) {my.grid.two.step <- my.grid.two.step[,-which(colnames(my.grid.two.step) == "threshold.flag.two.step")]}
		my.grid = cbind(my.grid, my.grid.two.step)
	}

  if(!test) {return(list("my.grid"=my.grid))}
  if(test) {
  	pval=c()
  	if(type == "model" | type == "both") {
  		pval = F.test(data.all, num.cov)
	  }
   	if(type == "two step" | type == "both") {
  		R.d.var.two.step = apply(as.matrix(boot.R.d.est.two.step), 1, mad)^2
  		std.R = (results.for.test$R.s.two.step - mean(results.for.test$R.s.two.step))/sqrt(R.d.var.two.step)
		  t.statistic = max(abs(std.R))
		  
		  aa.two.step = cor(t(boot.R.d.est.two.step))
		  undernull = mvtnorm::rmvnorm(1000, rep(0, length(results.for.test$R.s.two.step)), aa.two.step)
		  dist.two.step = apply(undernull, 1, function(x) max(abs(x)))
		  pval = c(pval,mean(dist.two.step >= t.statistic))
   	}
  return(list("my.grid" = my.grid, "pval" = pval))
	}
}
