# 
nrow(z)

lad.fit = function(x,y)
{
  n = length(y)
  p = ncol(x) 
  if (is.null(p)) p = 1
  yvec = rep(0,n*(n-1)/2)
  xmat = matrix(0,n*(n-1)/2, p)
  ix = 1
  for(i in 1:(n-1))
  {
    for (j in (i+1):n)
    {
      yvec[ix] = y[i]-y[j]
      xmat[ix,] = x[i,] - x[j,]
      ix = ix + 1
    }
  }
  return(rq.fit(xmat,yvec)$coef)
}

lad.fit(x =Z, y = Y) 
?rq.fit
fit<-rq.fit(zmat,yvec)
names(fit)
fit$coefficients
rq.fit(Z,Y)$coef

result_newton <- GEV_regfull(x=Y, z=Z, theta0=true_theta, beta0=rep(0,dim(Z)[2]), expr=expr_reg, alpha=0.5, maxiter=100)$root
# Finding MLE for nonstationary - fullhessian method
x=Y; z=Z; theta0=c(0,30,0.1);beta0=est_beta1; expr=expr_reg; alpha=0.5; maxiter=100

old_theta <- c(theta0,beta0)
Jaco <- expr$Jaco
Hmat <- expr$Hmat
mu=c(old_theta[1]+z%*%old_theta[-(1:3)]); s=old_theta[2]; k=old_theta[3]
  
  
evalh = eval(expr_reg$Hmat)
h1=matrix(apply(evalh,2,sum),3,3)
h2=rbind(apply(evalh[,1]*z,2,sum), apply(evalh[,2]*z,2,sum), apply(evalh[,3]*z,2,sum))
h3=t(h2)
# h4=apply((eval(expr_reg$Hmat)[,1]*z%*%t(z)),2,sum)
h4 = matrix(0,dim(z)[2],dim(z)[2])
for (i in 1:dim(z)[1]){
  h4row=evalh[,1][i] * z[i,]%*%t(z[i,])
  h4=h4+h4row
}
rbind(cbind(h1,h2),cbind(h3,h4)) 

#GEV_regfull <- function (x, z, theta0, beta0, expr=expr_reg, alpha=1, maxiter = 1000, tol = 1e-05) {
  old_theta <- c(theta0,
                 beta0)
  niter <- 0
  Jaco <- expr$Jaco
  Hmat <- expr$Hmat
  
  for (i in 1:maxiter) {
    niter = niter + 1
    mu=c(old_theta[1]+z%*%old_theta[-c(1,2,3)]); s=old_theta[2]; k=old_theta[3]
    
    grad=apply(cbind(eval(Jaco),eval(Jaco)[,1]*z),2,sum)
    if (max(abs(grad))<1e-07){
      cat("Grad in ML::", grad,'\n')
      cat("   ", i,'th step\n')
      break
    }
    hess = GEVhessian(x,z,mu,s,k)
    new_theta = old_theta - alpha*solve(hess)%*%grad
    
    old_theta = new_theta
  }
  return(list(initial = theta0, root = c(old_theta), step = niter))
}
