
###########################
#### Method Comparison ####
###########################

####### OLS #######

ols<-function(X, y, beta.true){
  X <- cbind(1, X)
  lm <- lm(data = data.frame(X, y), "y ~ x1 + x2 + x3 + D")
  su <- summary(lm)
  beta.hat <- su$coefficients[,"Estimate"]
  y.hat <- lm$fitted.values
  se <- su$coefficients[,"Std. Error"]
  RMSE <- sqrt(sum(su$residuals ^ 2) / length(y))
  if(nrow(su$coefficients) < 5){
    print("coefficients is empty")
  }
  coverage <- su$coefficients["D", 4]
  ATE <- beta.hat["D"]
  bias <- abs(ATE - beta.true) / beta.true * 100
  
  #95% CI of coefficient of D
  upper <- beta.hat["D"] + qnorm(0.975) * se["D"]
  lower <- beta.hat["D"] - qnorm(0.975) * se["D"]
  
  return(list(ATE = ATE, bias = bias, CI = c(lower, upper), RMSE = RMSE, coverage = coverage))
}

###### Causal Forest ######

forest<-function(X, y, beta.true){
  Y <- y
  W <- X[, "D"]
  X <- subset(X, select = c(x1, x2, x3))
  
  Y.forest <- regression_forest(X, Y)
  Y.hat <- predict(Y.forest)$predictions
  W.forest <- regression_forest(X, W)
  W.hat <- predict(W.forest)$predictions
  
  cf <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = W.hat)
  tau.hat <- predict(cf)$predictions
  ATE.forest <- average_treatment_effect(cf)
  tree <- get_tree(cf, 1)
  
  #95% CI of ATE
  upper <- ATE.forest[1] + qnorm(0.975) * ATE.forest[2]
  lower <- ATE.forest[1] - qnorm(0.975) * ATE.forest[2]
  
  RMSE <- sqrt(sum((y - Y.hat) ^ 2) / length(y))
  
  ATE <- ATE.forest[1]
  bias <- abs(ATE - beta.true) / beta.true * 100
  
  test <- test_calibration(cf)
  coverage <- test[1,4]
  return(list(ATE = ATE, bias = bias, CI = c(lower, upper), RMSE = RMSE, coverage = coverage, tree=tree, tau.hat=tau.hat))
}

###########################
####    Simulation     ####
###########################

simulation <- function(sim.times, data.size, dgp_func){
  plot.data<-data.frame(ATE = rep(NA, sim.times * length(data.size) * 2), 
                        RMSE = NA,
                        CI_lower = NA,
                        CI_upper = NA,
                        bias = NA,
                        coverage = NA,
                        method = factor(x=rep(c("ols", "forest"), sim.times * length(data.size)), levels = c("ols", "forest")), 
                        sim.time = rep(1:sim.times, each = 2, times = length(data.size)),
                        N = rep(data.size, each = 2 * sim.times))
  index <- 0
  for(j in 1:length(data.size)){
    N = data.size[j]
    for(i in 1:sim.times){
      sim.data <- dgp_func(N)
      
      index <- index + 1
      result.ols<-ols(sim.data$X, sim.data$y, sim.data$beta.true)
      plot.data$ATE[index] <- result.ols$ATE
      plot.data$RMSE[index] <- result.ols$RMSE
      plot.data$CI_lower[index] <- result.ols$CI[1]
      plot.data$CI_upper[index] <- result.ols$CI[2]
      plot.data$bias[index] <- result.ols$bias
      plot.data$coverage[index] <- result.ols$coverage
      
      index <- index + 1
      result.cf<- forest(sim.data$X, sim.data$y, sim.data$beta.true)
      plot.data$ATE[index] <- result.cf$ATE
      plot.data$RMSE[index] <- result.cf$RMSE
      plot.data$CI_lower[index] <- result.cf$CI[1]
      plot.data$CI_upper[index] <- result.cf$CI[2]
      plot.data$bias[index] <- result.cf$bias
      plot.data$coverage[index] <- result.cf$coverage
    }
  }
  return(plot.data)
}

###### CATE for model 3 & 4 ######
CATE <- function(dgp.data, x1_groups, x2_groups, model_name){
  
  x1_group_levels <- sort(unique(x1_groups))
  x1_value_levels <- map_dbl(x1_group_levels, function(x1_group) mean(dgp.data$X[x1_groups==x1_group,"x1"])) # value equals mean of x1 in its group
  x2_group_levels <- sort(unique(x2_groups))
  x2_value_levels <- map_dbl(x2_group_levels, function(x2_group) mean(dgp.data$X[x2_groups==x2_group,"x2"]))
  
  CATE.data <- data.frame(ATE = rep(NA, length(x1_group_levels) * length(x2_group_levels) * 2),
                        method = factor(x=rep(c("ols", "forest"), length(x1_group_levels) * length(x2_group_levels)), levels = c("ols", "forest")),
                        x2_group = rep(x2_group_levels, each = 2, times = length(x1_group_levels)),
                        x1_group = rep(x1_group_levels, each = 2 * length(x2_group_levels)),
                        x2_value = rep(x2_value_levels, each = 2, times = length(x1_value_levels)),
                        x1_value = rep(x1_value_levels, each = 2 * length(x2_value_levels)))
  CATE.matrix.forest <- matrix(NA, nrow = length(x1_value_levels), ncol = length(x2_value_levels), dimnames = list(x1_value_levels,x2_value_levels))
  CATE.matrix.ols <- matrix(NA, nrow = length(x1_value_levels), ncol = length(x2_value_levels), dimnames = list(x1_value_levels,x2_value_levels))
  CATE.matrix.beta_true <- matrix(NA, nrow = length(x1_value_levels), ncol = length(x2_value_levels), dimnames = list(x1_value_levels,x2_value_levels))
  index <- 1
  for (i in 1:length(x1_group_levels)) {
    for(j in 1:length(x2_group_levels)) {
      x1_group <- x1_group_levels[i]
      x2_group <- x2_group_levels[j]
      
      sub.X <- dgp.data$X[(x1_groups == x1_group & x2_groups == x2_group),]
      sub.y <- dgp.data$y[(x1_groups == x1_group & x2_groups == x2_group)]
      
      if(length(sub.y) > 4){
        if(sum(sub.X[,"D"]) != 0 && sum(sub.X[,"D"]) != nrow(sub.X)){
          result.ols <- ols(sub.X, sub.y, dgp.data$beta.true)
          result.cf <- forest(sub.X, sub.y, dgp.data$beta.true)
          
          CATE.data$ATE[index] <- result.ols$ATE
          CATE.data$ATE[index+1] <- result.cf$ATE
          
          CATE.matrix.ols[i,j] <- result.ols$ATE
          CATE.matrix.forest[i,j] <- result.cf$ATE
          
          if(model_name == "model 3"){
            CATE.matrix.beta_true[i,j] <- 0.02 + c(0.01, 0.001, 0.002) %*% c(mean(sub.X[,"x1"]), mean(sub.X[,"x2"]), mean(sub.X[,"x1"]*sub.X[,"x2"]))
          } else if (model_name == "model 4"){
            CATE.matrix.beta_true[i,j] <- 0.02 + c(0.1, 0.1, -0.02, -0.01, 0.01) %*% 
              c(mean(sub.X[,"x1"]), mean(sub.X[,"x2"]), mean(sub.X[,"x1"]^2), mean(sub.X[,"x2"]^2), mean(sub.X[,"x1"]*sub.X[,"x2"]))
          }
        }
      }
      index <- index + 2
    }
  }
  return(list(CATE.data = CATE.data, CATE.matrix.forest = CATE.matrix.forest, CATE.matrix.ols = CATE.matrix.ols, CATE.matrix.beta_true = CATE.matrix.beta_true, x1_value_levels = x1_value_levels, x2_value_levels = x2_value_levels))
}
