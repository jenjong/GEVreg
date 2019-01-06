rm( list = ls()); gc()
if (Sys.info()[1] == "Linux") {
  setwd("/home/jeon/Documents/GitHub/RankConsistency")
} else {
  setwd('C:/Users/Jeon/Documents/GitHub/RankConsistency')
}

library("Deriv")
library("evd")
library("quantreg")
source("fgevlibrary.R")
source("testlibrary.R")

n = 1000
p=10
true_theta=c(100,40,0.1)
true_beta = 10*c(c(1,1),rep(0,8))

result <- list()
for ( s in c(1:100) ) {
set.seed(s)
eps = rgev(n,loc=true_theta[1], scale=true_theta[2], shape=true_theta[3])
Z = matrix(rnorm(p*n),n,p)
Y = Z%*%true_beta + eps 

##### theta, beta MLE
# result_newton <- tryCatch(GEV_regfull(x=Y, z=Z, theta0=true_theta, beta0=rep(0,dim(Z)[2]), expr=expr_reg, alpha=1, maxiter=100),
#                           error=function(e){}, warning=function(e){})$root
result_newton <- GEV_regfull(x=Y, z=Z, theta0=true_theta, beta0=rep(0,dim(Z)[2]), expr=expr_reg, alpha=0.5, maxiter=100)$root

##### beta OLS, theta newton
est_beta1 <- lm( Y ~ Z )$coefficient
Y_tilda1 <- Y - Z %*% est_beta1[-1]
# est_beta1 <-  lm( Y-mean(eps) ~ Z )$coefficient

# Y_tilda1 <-   Y - cbind(rep(1,n),Z) %*% est_beta1
est_theta1 <- GEVnewtonRaphson(x=Y_tilda1,theta0=true_theta,expr=expr_mle,maxiter=1,step_theta=1)
result_lm <- c(est_theta1$root,est_beta1[-1])

##### beta QR, theta newton
est_beta2 <- rq( Y ~ Z )$coefficient
Y_tilda2 <- Y - Z %*% est_beta2[-1]
# est_beta2 <- rq( Y-median(eps) ~ Z )$coefficient
# Y_tilda2 <- Y - cbind(rep(1,n),Z) %*% est_beta2
est_theta2 <- GEVnewtonRaphson(x=Y_tilda2,theta0=true_theta,expr=expr_mle,maxiter=1,step_theta=1)
result_rq <- c(est_theta2$root,est_beta2[-1])

##### second approximation after lm&qr
result_lm_second <- GEV_regfull(x=Y,z=Z,theta0=est_theta1$root,beta0=est_beta1[-1],alpha=1,maxiter=1)$root
result_rq_second <- GEV_regfull(x=Y,z=Z,theta0=est_theta2$root,beta0=est_beta2[-1],alpha=1,maxiter=1)$root

result[[s]] <- rbind(result_newton,result_lm,result_lm_second,result_rq,result_rq_second)
}


# save(result,file="result_10.RData")
# load("./result_10.RData")

nrow_vec <- c()
for ( i in c(1:100)){
  nrow_vec[i] <- nrow(result[[i]])
}
all(nrow_vec==5)


#############################################
############### for boxplot #################
# result_z_all <- list()
# for (i in c(1:10)){
# z_num <- i+3
# result_z <- list()
# result_z[["newton"]] <- unlist(lapply(result,function(x) x[1,z_num]))
# result_z[["lm"]] <- unlist(lapply(result,function(x) x[2,z_num]))
# result_z[["lm+second"]] <- unlist(lapply(result,function(x) x[3,z_num]))
# result_z[["qr"]] <- unlist(lapply(result,function(x) x[4,z_num]))
# result_z[["qr+second"]] <- unlist(lapply(result,function(x) x[5,z_num]))
# result_z_all[[i]] <- result_z
# }
# par(mfrow=c(2,5))
# for ( i in c(1:10)){
# boxplot(result_z_all[[i]],main=paste0("beta",i))
# }


result_x_all <- list()
for (i in c(1:3)){
  result_x <- list()
  result_x[["newton"]] <- unlist(lapply(result,function(x) x[1,i]))
  result_x[["lm"]] <- unlist(lapply(result,function(x) x[2,i]))
  result_x[["lm+second"]] <- unlist(lapply(result,function(x) x[3,i]))
  result_x[["qr"]] <- unlist(lapply(result,function(x) x[4,i]))
  result_x[["qr+second"]] <- unlist(lapply(result,function(x) x[5,i]))
  result_x_all[[i]] <- result_x
}

par(mfrow=c(1,3))
boxplot(result_x_all[[1]],main="mu")
boxplot(result_x_all[[2]],main="alpha")
boxplot(result_x_all[[3]],main="kappa")


# for sigma, (hat-40)^2
par(mfrow=c(1,5))
for ( i in c(1:5)){
  boxplot( (unlist(result_x_all[[2]][i])-40)^2, main="alpha MSE" , ylim=c(0,10))
}

# for kappa, (hat-0.1)^2
par(mfrow=c(1,5))
for ( i in c(1:5)){
  boxplot( (unlist(result_x_all[[3]][i])-0.1)^2 ,main="kappa MSE" , ylim=c(0,0.005))
}

