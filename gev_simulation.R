rm(list = ls())
setwd("~/GitHub/GEVreg")
library("Deriv")
library("evd")
library("quantreg")
source("./lib/onestep_library.R")
source("./lib/fgevlibrary.R")

n = 2000
true_theta=c(100,40,0.25)
p= 10
true_beta = rep(0,p)
# optim control
optim_controlList = list()
optim_controlList$maxit = 1e+3

iter.num = 50
loc_mat  = scale_mat = shape_mat = matrix(0, iter.num, 3)
iter = 1
for (iter in 1:iter.num)
{
  cat('iter')
  seed=iter
  set.seed(seed)
  eps = rgev(n,loc=true_theta[1], scale=true_theta[2], shape=true_theta[3])
  Z = matrix(rnorm(n*p), n, p)
  Y = Z%*%true_beta + eps 
  fit_try0 = try(M <- fgev(x = drop(Y), nsloc = as.data.frame(Z))$estimate)
  ##### beta OLS, theta newton
  rdata = data.frame(Y,Z)
  est_beta_lm <- lm( Y ~ ., data = rdata)$coefficient[-1]
  Y_tilde_lm <- Y - as.matrix(Z) %*% est_beta_lm
  est_theta_lm <- fgev(Y_tilde_lm)$estimate
  result_lm <- c(est_theta_lm,est_beta_lm)
  ##### beta QR, theta newton
  est_beta_rq <- rq( Y ~ ., data = rdata)$coefficient[-1]
  Y_tilde_rq <- Y - Z %*% est_beta_rq
  est_theta_rq <- fgev(Y_tilde_rq)$estimate
  result_rq <- c(est_theta_rq,est_beta_rq)
  ##### second approximation after lm&qr
  fit_try1 = try(L2 <- 
                   GEV_regfull_step1(x=Y,z=Z, theta0=est_theta_lm,
                                        beta0=est_beta_lm, alpha=0.01)$root)
  fit_try2 = try(Q2 <- 
                   GEV_regfull_step1(x=Y,z=Z,theta0=est_theta_rq,
                                        beta0=est_beta_rq, alpha=1)$root)
  loc_mat[iter,]  = c(M[1], L2[1], Q2[1])
  scale_mat[iter,] = c(M[p+2], L2[2], Q2[2])
  shape_mat[iter,] = c(M[p+3], L2[3], Q2[3])
}

boxplot(loc_mat)
boxplot(scale_mat)
boxplot(shape_mat)

