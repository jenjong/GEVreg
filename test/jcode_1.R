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

sim.iter = 200
n = 100
true_theta=c(100,30,0.3)
p = 2
true_beta = c(1,-1, rep(0,p-2))

# compute the 11 by 8 matrix of basis function values
s = 1
vmat = NULL
for ( s in c(1:sim.iter) ) 
{
    set.seed(s)
    eps = rgev(n,loc=true_theta[1], scale=true_theta[2], shape=true_theta[3])
    Z <- matrix(rnorm(n*2),n,2)
    Y = Z%*%true_beta + eps 
    ##### theta, beta MLE
    fit <- GEV_regfull(x=Y, z=Z, theta0=true_theta, 
                       beta0=rep(0,p), expr=expr_reg, 
                       method = 'linear',                                   
                       Om = Om, lambda = lambda, 
                       alpha=0.1, maxiter=1000)
    
    v1 = norm(fit$root[-(1:3)]-true_beta, type = "2")
    fit = lm(Y~Z)
    v2 = norm(fit$coefficients[-1]-true_beta, "2")
    fit = rq(Y~Z)
    v3 = norm(fit$coefficients[-1] -true_beta, "2")
    vmat = rbind(vmat, c(v1,v2,v3))
}


boxplot(vmat, col = 'royalblue2', names = 
          c("MLE","LSE", "LAD"))

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
