two.step.est <-
function(data.control, data.treat, W.grid.expand.function) {
  extrapol.warn <- 0
  names(W.grid.expand.function) <- names(data.control)[3:ncol(data.control)]

   working.model.control <- lm(Y ~ ., data = data.control)
  working.model.treat <- lm(Y ~ ., data = data.treat)

  
  data.control.pred = data.control
  data.control.pred$S = mean(data.control$S)
  data.treat.pred = data.treat
  data.treat.pred$S = mean(data.control$S)
  control.u1 <- predict(working.model.treat, newdata= data.control.pred) - predict(working.model.control, newdata= data.control.pred) 
  treat.u1 <- predict(working.model.treat, newdata= data.treat.pred) - predict(working.model.control, newdata= data.treat.pred)
  W.grid.expand.function$S = mean(data.control$S)
 grid.u1 <- predict(working.model.treat, newdata= W.grid.expand.function) - predict(working.model.control, newdata= W.grid.expand.function)
 
 working.model.control <- lm(Y ~ ., data = data.control[,-which(names(data.control) == "S")])
 working.model.treat <- lm(Y ~ ., data = data.treat[,-which(names(data.control) == "S")])
 
 
control.u2 <- predict(working.model.treat, newdata= data.control) - predict(working.model.control, newdata= data.control) 
treat.u2 <- predict(working.model.treat, newdata= data.treat) - predict(working.model.control, newdata= data.treat)
  grid.u2 <- predict(working.model.treat, newdata= W.grid.expand.function) - predict(working.model.control, newdata= W.grid.expand.function)

all.u1 = c(control.u1, treat.u1)
all.u2 = c(control.u2, treat.u2)
u.pca = prcomp(cbind(all.u1,all.u2), scale = TRUE)
#Uhat is the first PCA
newdata = data.frame(cbind(control.u1, control.u2))
names(newdata) = c("all.u1", "all.u2")
data.control$U.hat = predict(u.pca, newdata = newdata)[,1]
newdata = data.frame(cbind(treat.u1, treat.u2))
names(newdata) = c("all.u1", "all.u2")
data.treat$U.hat = predict(u.pca, newdata = newdata)[,1]
newdata = data.frame(cbind(grid.u1, grid.u2))
names(newdata) = c("all.u1", "all.u2")
W.grid.expand.function$U.hat = predict(u.pca, newdata = newdata)[,1]
W.grid.expand.function = W.grid.expand.function[,-which((names(W.grid.expand.function) == "S"))]

  # kernel smoothing
  kernel <- function(x, h) { return(dnorm(x / h)) }
    
    h.1 = bw.nrd(data.treat$U.hat)*length(data.treat$U.hat)^(-0.14)
    h.0 = bw.nrd(data.control$U.hat)*length(data.control$U.hat)^(-0.14)
    h.3 = bw.nrd(data.treat$S)*length(data.treat$S)^(-0.14)
    
   get.delta <- function(u) {
    m.1 <- sum(kernel(data.treat$U.hat - u, h.1) * data.treat$Y) / sum(kernel(data.treat$U.hat - u, h.1))
    m.0 <- sum(kernel(data.control$U.hat - u, h.0) * data.control$Y) / sum(kernel(data.control$U.hat - u, h.0))
    return(m.1 - m.0)
  }

  get.delta.s <- function(u) {
    mus <- rep(NA, nrow(data.control))
    kernels <- rep(NA, nrow(data.control))
    for (i in 1:nrow(data.control)) {
      mus[i] <- sum(kernel(data.treat$U.hat - u, h.1) * kernel(data.treat$S - data.control$S[i], h.3) * data.treat$Y) / sum(kernel(data.treat$U.hat - u, h.1) * kernel(data.treat$S - data.control$S[i], h.3))
      kernels[i] <- kernel(data.control$U.hat[i] - u, h.0)
    }
    if (sum(is.na(mus))!=0) {
      ind = which(is.na(mus))
      sub.nona = cbind(data.control$S, mus)[-ind,]
      for (i in ind) {
        mat <- cbind(abs(sub.nona[,1] - data.control$S[i]), sub.nona[,2])
        mmm <- which(mat[, 1] == min(mat[, 1]))[1]
        mus[i] <- sub.nona[mmm,2]
      }
      extrapol.warn <- 1
    }
    m.10 <- sum(mus * kernels) / sum(kernels)
    m.0 <- sum(kernel(data.control$U.hat - u, h.0) * data.control$Y )/ sum(kernel(data.control$U.hat - u, h.0))
    return(m.10-m.0)
  }

  my.grid <- W.grid.expand.function

  my.grid$delta.two.step <- unlist(lapply(my.grid$U.hat, get.delta))
  my.grid$delta.s.two.step <- unlist(lapply(my.grid$U.hat, get.delta.s))
  my.grid$R.s.two.step <- 1 - (my.grid$delta.s.two.step / my.grid$delta.two.step)

  return(list("my.grid" = my.grid, "extrapol.warn" = extrapol.warn))

}
