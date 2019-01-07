# Finding MLE for nonstationary - fullhessian method
x=Y; z=Z; theta0=c(100,30,0.1);beta0=rep(0,ncol(Z)); expr=expr_reg; alpha=1e-2; maxiter=10000
method = 'B-spline' ; lambda = 0.01

n = length(x)
old_theta <- c(theta0,beta0)
#if (method == 'B-spline' ) old_theta[1] = 0
new_theta <- old_theta
niter <- 0
Jaco <- expr$Jaco
Hmat <- expr$Hmat
for (i in 1:maxiter)
{
    # if the spline is used, the intercept is not updated for identifiability
    niter = niter + 1
    mu = c(old_theta[1]+z%*%old_theta[-(1:3)])
    s  = old_theta[2]
    k  = old_theta[3]
    
    grad = apply(cbind(eval(Jaco),eval(Jaco)[,1]*z),2,mean)
    if (method == 'B-spline' )  
    {
      grad[-(1:3)] = grad[-(1:3)] + lambda*Om%*%old_theta[-(1:3)]*n
      grad = grad[-1]
    }
    if (max(abs(grad))<1e-07){
      break
    }
    
    hess = GEVhessian(x,z,mu,s,k) 
    hess = hess + diag(1e-8, ncol(hess))
    min.h = 1
    if (method == 'B-spline' ) 
    {
      hess[-(1:3),-(1:3)] = hess[-(1:3),-(1:3)] + lambda*Om*n
      hess = hess[-1,-1]
      gd = -solve(hess)%*%grad
      new_theta[-1] = old_theta[-1] + alpha*gd
      if (new_theta[3]>0)
      {
        # check alpha 
        if (new_theta[2] < 0 )
        {
          min.h = (1 - old_theta[2] )/gd[1]
#          stop ("alpha is unstable")
        }  
        # check mu and kappa 
        mu.i = new_theta[1] +Z%*%new_theta[-(1:3)]
        lx = (min(x) <= mu.i - new_theta[2]/new_theta[3])
        if (any(lx))  
        {
          hvec = 2^((-20):0)
          
          new_theta[1] +Z%*%new_theta[-(1:3)]
          stop ("mu and kappa is unstable")
        }
      }
      
      if (new_theta[3]<0)
      {
        # check alpha 
        if (new_theta[2] < 0 )  stop ("alpha is unstable")
        # check mu and kappa 
        mu.i = new_theta[1] +Z%*%new_theta[-(1:3)]
        lx = (max(x) >= mu.i - new_theta[2]/new_theta[3])
        if (any(lx))  stop ("mu and kappa is unstable")
      }
    }
    
    if (method == 'linear' ) 
    {
      new_theta = old_theta - alpha*solve(hess)%*%grad  
    }
  old_theta = new_theta
}



plot(Z%*%new_theta[-(1:3)])
