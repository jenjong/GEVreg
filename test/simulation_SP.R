rm( list = ls()); gc()
if (Sys.info()[1] == "Linux") {
  setwd("/home/jeon/Documents/GitHub/GEVreg")
} else {
  setwd('C:/Users/Jeon/Documents/GitHub/GEVreg')
}
library("fda")
library("Deriv")
library("evd")
library("quantreg")
source("./test/fgevlibrary.R")
source("./test/testlibrary.R")

n = 100
true_theta=c(100,30,0.1)
true_beta = 20*sin(seq(0, 2*pi, length = n))

# compute the 11 by 8 matrix of basis function values
x <- 1:n
norder <- 5
breaks <- seq(0,100,20)
Z <- bsplineS(x, breaks, norder)
p = ncol(Z)
basisobj = create.bspline.basis(rangeval=c(1, n), nbasis=NULL,
                                norder=norder, breaks= breaks, dropind=NULL,
                                quadvals=NULL, values=NULL)
Om = eval.penalty(basisobj, Lfdobj=int2Lfd(2), rng=c(1, n))
s = 3
sim.iter = 100
lambda_vec = seq(0, 1, length  = 50)
#for ( s in c(1:sim.iter) ) {
  set.seed(s)
  eps = rgev(n,loc=true_theta[1], scale=true_theta[2], shape=true_theta[3])
  Y = true_beta + eps 
  ##### theta, beta MLE
    lambda = 0.01
    fit <- GEV_regfull(x=Y, z=Z, theta0=true_theta, 
                                     beta0=rep(0,p), expr=expr_reg, 
                                     method = 'B-spline',                                   
                                     Om = Om, lambda = lambda, 
                                     alpha=0.1, maxiter=1000)
  
  #fit$
  plot(Z%*%fit$root[-(1:3)])
  
  plot(Y)
  lines(fit$root[1]+Z%*%fit$root[-(1:3)])
  
  
  # 
  lambda = 10
  plot(Z%*%solve(t(Z)%*%Z + n*lambda*Om)%*%t(Z)%*%Y)
  Y
    
  plot(v)
  lambda = lambda_vec[  which.max(v) ]
  est_beta1 = solve(t(Z)%*%Z + n*lambda*Om)%*%t(Z)%*%Y
  Y_tilda1 <- Y - Z %*% est_beta1
  est_theta1 <- fgev(Y_tilda1)$estimate

  j = 1
  for (j in 1:length(lambda_vec))
  {
    lambda = lambda_vec[j]
    
    
  }
  
  plot(Z %*% result_newton[-(1:3)])
  
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



install.packages("fda")
library(fda)
x <- seq(0,100,1)
# the order willl be 4, and the number of basis functions
# is equal to the number of interior break values (4 here)
# plus the order, for a total here of 8.
norder <- 5
breaks <- seq(0,100,20)
# compute the 11 by 8 matrix of basis function values
basismat <- bsplineS(x, breaks, norder)
plot(basismat[,4])
bsplinepen {fda}
basisobj <- create.bspline.basis(c(0,1),13)
basisobj = create.bspline.basis(rangeval=c(0, 100), nbasis=NULL,
                     norder=4, breaks= breaks, dropind=NULL,
                     quadvals=NULL, values=NULL)
eval.penalty(basisobj, Lfdobj=int2Lfd(0), rng=c(0, 100))
