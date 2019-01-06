rm( list = ls()); gc()
if (Sys.info()[1] == "Linux") {
  setwd("/home/jeon/Documents/GitHub/GEVreg")
} else {
  setwd('C:/Users/Jeon/Documents/GitHub/GEVreg')
}

library("Deriv")
library("evd")
library("quantreg")
source("./test/fgevlibrary.R")
source("./test/testlibrary.R")

n = 100
p=2
true_theta=c(100,40,0.1)
true_beta = 10*rep(1,p)

result <- list()
s = 1
sim.iter = 1000
for ( s in c(1:sim.iter) ) {
set.seed(s)
eps = rgev(n,loc=true_theta[1], scale=true_theta[2], shape=true_theta[3])
Z = matrix(runif(p*n,-1,1),n,p)
#Z[,1] = (1:n)/2
#Z = sweep(Z,2,colMeans(Z))
Y = Z%*%true_beta + eps 

##### theta, beta MLE
##### beta OLS, theta newton
est_beta1 <- lm( Y ~ Z )$coefficient[-1]
Y_tilda1 <- Y - Z %*% est_beta1
est_theta1 <- fgev(Y_tilda1)$estimate
result_lm <- c(est_theta1,est_beta1)

##### beta QR, theta newton
est_beta2 <- rq( Y ~ Z )$coefficient[-1]
#est_beta2 = lad.fit(x = Z, y = Y)
Y_tilda2 <- Y - Z %*% est_beta2
est_theta2 <- fgev(Y_tilda2)$estimate
result_rq <- c(est_theta2,est_beta2)


# result_newton <- tryCatch(GEV_regfull(x=Y, z=Z, theta0=true_theta, beta0=rep(0,dim(Z)[2]), expr=expr_reg, alpha=1, maxiter=100),
#                           error=function(e){}, warning=function(e){})$root
result_newton <- try(GEV_regfull(x=Y, z=Z, theta0=true_theta, 
                                      beta0=est_beta2, expr=expr_reg, 
                                      alpha=0.5, maxiter=100)$root)
fgev(Y - Z %*% result_newton[-(1:3)])
if (class(result_newton) == 'try-error') next

##### second approximation after lm&qr
sfit_1 = GEV_regfull(x=Y,z=Z,theta0=est_theta1,beta0=est_beta1, alpha=1, maxiter=1)
cat("sfit1::", tail(eigen(sfit_1$hess)$values,1), '\n')
sfit_2 = GEV_regfull(x=Y,z=Z,theta0=est_theta2,beta0=est_beta2,alpha=1,maxiter=1)
cat("sfit2::", tail(eigen(sfit_2$hess)$values,1), '\n') 
result_lm_second <- GEV_regfull(x=Y,z=Z,theta0=est_theta1,beta0=est_beta1, alpha=1, maxiter=1)$root
# v1 = sum(dgev(Y, est_theta1[1] + Z%*%est_beta1, est_theta1[2], est_theta1[3], log = T) )
# v2 = sum(dgev(Y, result_lm_second[1] + Z%*%result_lm_second[-(1:3)], 
#          result_lm_second[2], result_lm_second[3], log = T) )

result_rq_second <- GEV_regfull(x=Y,z=Z,theta0=est_theta2,beta0=est_beta2,alpha=1,maxiter=1)$root
# v1 = sum(dgev(Y, est_theta2[1] + Z%*%est_beta2, est_theta2[2], est_theta2[3], log = T) )
# v2 = sum(dgev(Y, result_rq_second[1] + Z%*%result_rq_second[-(1:3)], 
#          result_rq_second[2], result_rq_second[3], log = T) )
# if (v1>v2) next
#Z = matrix(runif(p*n*10,-3,3),n*10,p)
#est = result_rq_second

result[[s]] <- rbind(result_newton,result_lm,result_lm_second,result_rq,result_rq_second)
}


# save(result,file="result_10.RData")
# load("./result_10.RData")

nrow_vec <- c()
for ( i in c(1:sim.iter)){
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
boxplot(result_x_all[[1]][-3],main="mu", ylim = c(80,120))
abline(h = 100, col = 'red')
boxplot(result_x_all[[2]][-3],main="alpha", ylim = c(30,50))
abline(h = 40, col = 'red')
boxplot(result_x_all[[3]][-3],main="kappa")
abline(h = 0.1, col = 'red')



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




