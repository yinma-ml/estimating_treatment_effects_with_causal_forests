# Tips: simulation1 and simulation2 takes hours to run. #
# All the resultes are saved in Pics File #
# You can change data.size and sim.times into smaller numbers for quick look. #

rm(list = ls())

if (!require("grf")) install.packages("grf")
if (!require("DiagrammeR")) install.packages("DiagrammeR")
if (!require("akima")) install.packages("akima")
if (!require("plotly")) install.packages("plotly")
library(MASS)
library(ggplot2)
library(dplyr)
library(purrr)
library(cowplot)

source("functions.R")

set.seed(100)


#######################
#####     DGP     #####
#######################

######## model 1 ########
# y = 0.05x1 - 0.005x2 + 0.01x3 + 0.02D + eps #

dgp_model1 <- function(N){
  # mean, standard deviation, correlation matrix for generate (x1,x2,x3,d) using mvrnorm
  mu <- c(1.40, 5.50, 0.10, 0.20)
  sd <- c(1.00, 1.50, 0.30, 0.30)
  cor_matrix <- matrix(c( 1.00,-0.05,-0.30, 0.15,
                          -0.05, 1.00, 0.20, 0.35,
                          -0.30, 0.20, 1.00, 0.50,
                          0.15, 0.35, 0.50, 1.00), nrow = 4, ncol = 4)
  beta <- c(0.05, -0.005, 0.01, 0.02) # for x1,x2,x3,D
  # gerate covariation matrix using standard error and correlation matrix
  cov_matrix <- diag(sd) %*% cor_matrix %*% diag(sd)
  # generate x1,x2,x3,d
  mvnorm_data <- mvrnorm(n = N, mu = mu, Sigma = cov_matrix)
  # generate X : (x1, x2, x3)
  poly_max <- 1
  X <- reduce( c(1,1:poly_max), 
               function(X, poly)
                 cbind(X, mvnorm_data[,1]^poly, mvnorm_data[,2]^poly, mvnorm_data[,3]^poly)
  ) %>% subset(select = -c(1))
  if(poly_max == 1) {
    colnames(X) <- c("x1","x2","x3")
  } else {
    colnames(X) <- rep(c("x1","x2","x3"), poly_max - 1) %>%
      paste0(c(rep("", 3),rep(paste0("^", 2:poly_max), each = 3)))
  }
  
  D <- as.integer(mvnorm_data[,4] > 0)
  X <- cbind(X, D)
  eps <- rnorm(N, 0, 0.065)
  y <- X %*% beta + eps
  X <- subset(X, select = c(x1, x2, x3, D))
  
  beta.true <- 0.02
  return(list(X=X, y=y, d=mvnorm_data[,4], beta.true=beta.true))
}

######## model 2 ########
# y = 0.05x1 + 0.005x2 + 0.01x3 + 0.025x1^2 + 0.01x2^2 + 0.015x3^2 + 0.02D + eps 

dgp_model2 <- function(N){
  # mean, standard deviation, correlation matrix for generate (x1,x2,x3,d) using mvrnorm
  mu <- c(1.40, 5.50, 0.10, 0.20)
  sd <- c(1.00, 1.50, 0.30, 0.30)
  cor_matrix <- matrix(c( 1.00,-0.05,-0.30, 0.15,
                          -0.05, 1.00, 0.20, 0.35,
                          -0.30, 0.20, 1.00, 0.50,
                          0.15, 0.35, 0.50, 1.00), nrow = 4, ncol = 4)
  beta <- c(0.05, -0.005, 0.01, 0.025, -0.01, 0.015, 0.02) # for x1,x2,x3,x1^2,x2^2,x3^2,D
  # gerate covariation matrix using standard error and correlation matrix
  cov_matrix <- diag(sd) %*% cor_matrix %*% diag(sd)
  # generate x1,x2,x3,d
  mvnorm_data <- mvrnorm(n = N, mu = mu, Sigma = cov_matrix)
  # generate X with polynomial: (x1, x2, x3, x1^2, x2^2, x3^2)
  poly_max <- 2
  X <- reduce( c(1,1:poly_max), 
               function(X, poly)
                 cbind(X, mvnorm_data[,1]^poly, mvnorm_data[,2]^poly, mvnorm_data[,3]^poly)
  ) %>% subset(select = -c(1))
  if(poly_max == 1) {
    colnames(X) <- c("x1","x2","x3")
  } else {
    colnames(X) <- rep(c("x1","x2","x3"), poly_max - 1) %>%
      paste0(c(rep("", 3),rep(paste0("^", 2:poly_max), each = 3)))
  }
  
  D <- as.integer(mvnorm_data[,4] > 0)
  X <- cbind(X, D)
  eps <- rnorm(N, 0, 0.065)
  y <- X %*% beta + eps
  X <- subset(X, select = c(x1, x2, x3, D))
  
  beta.true <- 0.02
  return(list(X=X, y=y, d=mvnorm_data[,4], beta.true=beta.true))
}

######## model 3 ########
# y = 0.05x1 + 0.005x2 + 0.01x3 + 0.02D + tau1*x1*D + tau2*x2*D + tau3*x1*x2*D + eps

dgp_model3 <- function(N){
  # mean, standard deviation, correlation matrix for generate (x1,x2,x3,d) using mvrnorm
  mu <- c(1.40, 5.50, 0.10, 0.20)
  sd <- c(1.00, 1.50, 0.30, 0.30)
  cor_matrix <- matrix(c( 1.00,-0.05,-0.30, 0.15,
                          -0.05, 1.00, 0.20, 0.35,
                          -0.30, 0.20, 1.00, 0.50,
                          0.15, 0.35, 0.50, 1.00), nrow = 4, ncol = 4)
  beta <- c(0.05, -0.005, 0.01, 0.02) # for x1,x2,x3,D
  tau <- c(0.01, 0.001, 0.002)
  # gerate covariation matrix using standard error and correlation matrix
  cov_matrix <- diag(sd) %*% cor_matrix %*% diag(sd)
  # generate x1,x2,x3,d
  mvnorm_data <- mvrnorm(n = N, mu = mu, Sigma = cov_matrix)
  # generate X with polynomial: (x1, x2, x3)
  poly_max <- 1
  X <- reduce( c(1,1:poly_max), 
               function(X, poly)
                 cbind(X, mvnorm_data[,1]^poly, mvnorm_data[,2]^poly, mvnorm_data[,3]^poly)
  ) %>% subset(select = -c(1))
  if(poly_max == 1) {
    colnames(X) <- c("x1","x2","x3")
  } else {
    colnames(X) <- rep(c("x1","x2","x3"), poly_max - 1) %>%
      paste0(c(rep("", 3),rep(paste0("^", 2:poly_max), each = 3)))
  }
  
  D <- as.integer(mvnorm_data[,4] > 0)
  X <- cbind(X, D)
  eps <- rnorm(N, 0, 0.065)
  y <- X %*% beta + eps
  
  # add intersection
  X_D <- cbind(x1_D = X[,"x1"] * D, x2_D = X[,"x2"] * D, x1_x2_D = X[,"x1"]*X[,"x2"]*D)
  y <- y + X_D %*% tau
  
  X <- subset(X, select = c(x1, x2, x3, D))
  
  beta.true <- 0.02 + tau %*% c(mu[1:2], (mu[1]*mu[2] + cov_matrix[1,2])) # E(X1)*E(X2)+Cov(X1,X2)=E(X1X2)
  return(list(X=X, y=y, d=mvnorm_data[,4], beta.true=beta.true))
}

######## model 4 ########
# y = 0.05x1 + 0.005x2 + 0.01x3 + 0.02D + tau1*x1*D + tau2*x2*D + tau3*x1^2*D + tau4*x2^2*D + tau5*x1*x2*D + eps

dgp_model4 <- function(N){
  # mean, standard deviation, correlation matrix for generate (x1,x2,x3,d) using mvrnorm
  mu <- c(1.40, 5.50, 0.10, 0.20)
  sd <- c(1.00, 1.50, 0.30, 0.30)
  cor_matrix <- matrix(c( 1.00,-0.05,-0.30, 0.15,
                          -0.05, 1.00, 0.20, 0.35,
                          -0.30, 0.20, 1.00, 0.50,
                          0.15, 0.35, 0.50, 1.00), nrow = 4, ncol = 4)
  beta <- c(0.05, -0.005, 0.01, 0.02) # for x1,x2,x3,D
  tau <- c(0.1, 0.1, -0.02, -0.01, 0.01)
  # gerate covariation matrix using standard error and correlation matrix
  cov_matrix <- diag(sd) %*% cor_matrix %*% diag(sd)
  # generate x1,x2,x3,d
  mvnorm_data <- mvrnorm(n = N, mu = mu, Sigma = cov_matrix)
  # generate X: (x1, x2, x3)
  poly_max <- 1
  X <- reduce( c(1,1:poly_max), 
               function(X, poly)
                 cbind(X, mvnorm_data[,1]^poly, mvnorm_data[,2]^poly, mvnorm_data[,3]^poly)
  ) %>% subset(select = -c(1))
  if(poly_max == 1) {
    colnames(X) <- c("x1","x2","x3")
  } else {
    colnames(X) <- rep(c("x1","x2","x3"), poly_max - 1) %>%
      paste0(c(rep("", 3),rep(paste0("^", 2:poly_max), each = 3)))
  }
  
  D <- as.integer(mvnorm_data[,4] > 0)
  X <- cbind(X, D)
  eps <- rnorm(N, 0, 0.065)
  y <- X %*% beta + eps
  
  # add intersection
  X_D <- cbind(x1_D = X[,"x1"] * D, x2_D = X[,"x2"] * D, x12_D = X[,"x1"]^2 * D, x22_D = X[,"x2"]^2 * D, x1_x2_D = X[,"x1"] * X[,"x2"] * D)
  y <- y + X_D %*% tau
  
  X <- subset(X, select = c(x1, x2, x3, D))
  
  beta.true <- 0.02 + tau %*% c(mu[1:2], (mu[1:2]^2 + sd[1:2]^2), (mu[1]*mu[2] + cov_matrix[1,2]))
  return(list(X=X, y=y, d=mvnorm_data[,4], beta.true=beta.true))
}


###########################
####   Visualization   ####
###########################

sim_model <- list(dgp_model1, dgp_model2, dgp_model3, dgp_model4)
model_names <- c("model 1", "model 2", "model 3", "model 4")

#### plot tree and histogram ####
for(i in c(1:length(sim_model))){
  # plot tree
  dgp_func <- sim_model[[i]]
  name <- model_names[i]
  sim.data <- dgp_func(1000)
  result.cf<- forest(sim.data$X, sim.data$y, sim.data$beta.true)
  print(plot(result.cf$tree))

  # Evaluate the estimate using a histogram
  png(paste0(project_name, "_", name, "_histogram.png"), width=600*4, height=600*4, res=72*4)
  print(ggplot(data.frame(effect=result.cf$tau.hat, y=sim.data$y, treated=sim.data$X[,"D"]), aes(x=effect)) +
    geom_histogram(bins=60) +
    geom_vline(xintercept=result.cf$ATE, linetype="dotted", color="blue", size=1.5) +
    geom_vline(xintercept=sim.data$beta.true, color = "red", size=1.5))
    labs(title = paste0(name, " tau.hat histogram"))
  dev.off()
}


#### simulation1: modify data size ####
data.size1 <- seq(500, 15000, 100)
sim.data.list1 <- list()
for(i in c(1:length(sim_model))){
  dgp_func <- sim_model[[i]]
  name <- model_names[i]
  print(paste0("sim1 ", name, " start"))
  
  sim.data <- simulation(sim.times = 1, data.size = data.size1, dgp_func)
  sim.data.list1[[i]] <- sim.data
  write.csv(sim.data, paste0(project_name, "_", name, "_sim1.csv"), row.names = FALSE)
  # sim.data <- sim.data.list1[[i]] # for replay
  # sim.data <- read.csv(paste0(project_name, "_", name, "_sim1.csv")) # for replay
  beta.true <- dgp_func(2)$beta.true

  png(paste0(project_name, "_", name, "_bias_line.png"), width=600*4, height=600*4, res=72*4)
  plot_bias <- ggplot(sim.data, aes(N)) +
    geom_line(aes(y = bias, color = method)) +
    ylim(min(sim.data$bias), max(sim.data$bias)) +
    labs(title = paste0(name, " bias"))
  print(plot_bias)
  dev.off()

  png(paste0(project_name, "_", name, "_ATE_line.png"), width=600*4, height=600*4, res=72*4)
  plot_ATE <- ggplot(sim.data, aes(N)) +
    geom_line(aes(y = ATE, color = method), cex = 1) +
    geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper, fill = method), alpha = 0.1) +
    geom_hline(yintercept = beta.true, color = I("black")) +
    ylim(min(sim.data$CI_lower), max(sim.data$CI_upper)) +
    labs(title = paste0(name, " ATE"))
  print(plot_ATE)
  dev.off()
  
  print(paste0("sim1 ", name, " end"))
}


#### simulation2: iteration for every data size ####
data.size2 <- c(500, 1000, 5000, 10000)
sim.data.list2 <- list()
mean_sim2 <- data.frame()
for(i in c(1:length(sim_model))){
  dgp_func <- sim_model[[i]]
  name <- model_names[i]
  
  print(paste0("sim2 ", name, " start"))
  sim.data <- simulation(sim.times = 1000, data.size = data.size2, dgp_func)
  sim.data.list2[[i]] <- sim.data
  write.csv(sim.data, file = paste0(project_name, "_", name, "_sim2.csv"), row.names = FALSE)
  # sim.data <- sim.data.list2[[i]] # for replay
  # sim.data <- read.csv(paste0(project_name, "_", name, "_sim2.csv")) # for replay
  beta.true <- dgp_func(2)$beta.true
  
  # generate point and box plot for every data size
  ATE_point_plot.list <- list()
  ATE_box_plot.list <- list()
  for(j in 1:length(data.size2)){
    N <- data.size2[j]
    plot.heigth <- max(beta.true - min(sim.data$ATE), max(sim.data$ATE) - beta.true)
    ATE_point_plot.list[[j]] <- ggplot(sim.data[sim.data$N==N,], aes(sim.time)) +
      geom_point(aes(y = ATE, color = method, shape = method)) +
      geom_hline(yintercept = beta.true, color = I("black")) +
      ylim(beta.true - plot.heigth, beta.true + plot.heigth) +
      labs(subtitle = paste0("sample size = ", N), x = "iteration", y = "ATE")
    
    ATE_box_plot.list[[j]] <- ggplot(sim.data[sim.data$N==N,], aes(method)) +
      geom_boxplot(aes(y = ATE, color = method)) +
      ylim(min(sim.data$ATE), max(sim.data$ATE)) +
      labs(subtitle = paste0("sample size = ", N), x = NULL)
  }
  title <- ggdraw() + 
    draw_label(name, fontface = 'bold', x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 20))
  ATE_point_plot <- plot_grid(plotlist = ATE_point_plot.list)
  ATE_box_plot <- plot_grid(plotlist = ATE_box_plot.list)
  png(paste0(project_name, "_", name, "_ATE_point.png"), width=600*4, height=600*4, res=72*4)
  print(plot_grid(title, ATE_point_plot, ncol = 1, rel_heights = c(0.1, 1)))
  dev.off()
  png(paste0(project_name, "_", name, "_ATE_box.png"), width=600*4, height=600*4, res=72*4)
  print(plot_grid(title, ATE_box_plot, ncol = 1, rel_heights = c(0.1, 1)))
  dev.off()
  
  # generate mean of simulation result data
  sim.data_grouped <- group_by(sim.data, method, N)
  mean_by_method_N <- summarize(sim.data_grouped, ATE = mean(ATE,na.rm=TRUE), RMSE = mean(RMSE,na.rm=TRUE), bias = mean(bias,na.rm=TRUE),
                                       coverage = mean(coverage,na.rm=TRUE)) %>% mutate(model = name)
  mean_sim2 <- rbind(mean_sim2, mean_by_method_N)
  
  print(paste0("sim2 ", name, " end"))
}
write.csv(mean_sim2, paste0(project_name, "_sim2_mean.csv"), row.names = FALSE)

# box plot all in one for sim2
select.model <- model_names
select.data.size <- c(1000, 5000, 10000)
model.box_plot.list <- list()
for(model_name in select.model){
  sim.data <- sim.data.list2[[which(model_name == model_names)]]
  sim.data <- sim.data[sim.data$N %in% select.data.size,]
  model.box_plot.list[[I(model_name)]] <- ggplot(sim.data, aes(method)) +
    geom_boxplot(aes(y = ATE, color = method)) +
    ylim(min(sim.data$ATE), max(sim.data$ATE)) +
    facet_wrap(~N, nrow = 1) +
    labs(title = model_name, x = NULL)
}
pdf(paste0(project_name, "_all_box.pdf"), width=16, height=9)
print(plot_grid(plotlist = model.box_plot.list, align = "hv"))
dev.off()

#### CATE simulation ####
set.seed(200)
CATE_sim_model <- list(dgp_model3, dgp_model4)
CATE_model_names <- c("model 3", "model 4")
for(i in c(1:length(CATE_sim_model))){
  dgp_func <- CATE_sim_model[[i]]
  name <- CATE_model_names[i]
  
  dgp.data <- dgp_func(500000)
  CATE.result <- CATE(dgp.data, ntile(dgp.data$X[,"x1"],50), ntile(dgp.data$X[,"x2"],50), name)
  
  # rainbow plot
  png(paste0(project_name, "_", name, "_CATE_rainbow_cf.png"), width=600*4, height=600*4, res=72*4)
  print(ggplot(subset(CATE.result$CATE.data, method=="forest"), aes(x=x1_group, y=x2_group, color = ATE)) +
          geom_point(size = 4) +
          scale_color_gradientn(colours = rainbow(10)) +
          labs(title = paste0(name, " CATE rainbow of forest")))
  dev.off()
  
  # box plot
  box.plot.data <- CATE.result$CATE.data %>% subset(ATE < 2 & method=="forest") %>% mutate(x1_value = as.factor(.$x1_value))
  png(paste0(project_name, "_", name, "_CATE_box_cf.png"), width=600*4, height=600*4, res=72*4)
  ggplot(box.plot.data, aes(x1_value)) +
    geom_boxplot(aes(y = ATE)) +
    labs(title = paste0(name, " CATE boxplot of forest"))
  dev.off()
  
  # Gridded Bivariate Interpolation plot by package "akima"
  png(paste0(project_name, "_", name, "_CATE_interp_cf.png"), width=600*4, height=600*4, res=72*4)
  im <- with(CATE.result$CATE.data %>% subset(method == "forest" & !is.na(ATE)), interp(x1_value,x2_value,ATE))
  with(im,image(x,y,z, xlab = "x1 value", ylab = "x2 value"))
  dev.off()
  
  # 3D plot by package "plotly"
  plotly_forest <- plot_ly(x=CATE.result$x1_value_levels,y=CATE.result$x2_value_levels,z=CATE.result$CATE.matrix.forest, type="surface")
  htmlwidgets::saveWidget(plotly_forest, file = paste0(project_name, "_", name, "_CATE_3D_cf.html"))
  plotly_ols <- plot_ly(x=CATE.result$x1_value_levels,y=CATE.result$x2_value_levels,z=CATE.result$CATE.matrix.ols, type="surface")
  htmlwidgets::saveWidget(plotly_ols, file = paste0(project_name, "_", name, "_CATE_3D_ols.html"))
  plotly_true <- plot_ly(x=CATE.result$x1_value_levels,y=CATE.result$x2_value_levels,z=CATE.result$CATE.matrix.beta_true, type="surface")
  htmlwidgets::saveWidget(plotly_true, file = paste0(project_name, "_", name, "_CATE_3D_true.html"))
}
