CI.delta <- function(x,GP,alpha){
  # CI for GPD survival function using delta method, u is the threshold, V is the covariance matrix
  u <- GP$threshold
  sig <- exp(GP$results$par[1])
  pu <- sum(GP$x>GP$threshold)/GP$n
  if (GP$type == "GP"){
    gamma <- GP$results$par[2]
    aa <- 1 + gamma * (x-u)/sig
    prob <- Tail.prob.GP(x,model = GP)
    d.sig <- pu* (x-u)/sig^2 * aa^(-1/gamma-1)
    d.gamma <-  prob * (log(aa)/gamma^2 - x/(gamma*sig*aa) )
    grad <- c(d.sig,d.gamma)
    V <- solve(GP$results$hessian,diag(c(1,1)))
    return(c(prob-qnorm(1-alpha/2)* sqrt(t(grad) %*% V %*% grad),prob,prob+qnorm(1-alpha/2)* sqrt(t(grad) %*% V %*% grad)) )}
  else if (GP$type == "Exponential") {
    # For the Exponential distribution
    d.sig <- -exp(-(x - u) / sig) * (x - u) / sig^2
    grad <- c(d.sig)
    V <- solve(GP$results$hessian, diag(1)) # Only one parameter in Exponential distribution
    return(c(
      exp(-(x - u) / sig) - 1.96 * sqrt(t(grad) %*% V %*% grad),
      exp(-(x - u) / sig),
      exp(-(x - u) / sig) + 1.96 * sqrt(t(grad) %*% V %*% grad)
    ))
  } else {
    stop("Unsupported GP$type. Please specify either 'GP' or 'Exponential'.")
  }
  
}
CI.prolik <- function(x0,GP,alpha=0.95,expand.factor=2,nint=1000,if.plot=FALSE,use.phi=FALSE){
  # x0 is the value of the event at which the probability is of interest
  # GP is the MLE GP model
  # p.guess is the starting value for maximizing the profile likelihood w.r.t prob p
  u <- GP$threshold
  data <- GP$x[GP$x>u]
  p0 <- Tail.prob.GP(x0,model = GP,use.phi = use.phi)
  p.guess <- c(p0*1e-2,p0*10)
  if (p0==0){
    stop("target probability is 0")
  }
  max.iter <- 8
  iter <- 0
  repeat{
    p <- seq(p.guess[1],p.guess[2],length=ifelse(p.guess[2]/p.guess[1]<1e-4,nint,10*nint))
    
    theta <- GP$results$par[2]
    pu <- sum(GP$x>GP$threshold)/GP$n
    
    # profile likelihood of l(gamma)
    l.p.gamma <-function(theta,x0,data,pu,p0,u){
      sig <- theta*(x0-u)/((p0/pu)^(-theta)-1)
      y <- 1 + theta * (data-u)/sig
      if (any(y <= 0)){
        return(1e-6)
      }
      else{
        return(ifelse(use.phi,length(data-u) * (log(pu) - sig) - (theta+1)/theta * sum(log(y)),
                      length(data-u) * (log(pu) - log(sig)) - (theta+1)/theta * sum(log(y))) )
      }
    }
    
    # find argmax_p of l.p.gamma
    l <- sapply(p, function(p0) {
      output <- optim(theta, fn = l.p.gamma, method = "BFGS", # change to BFGS if shit happend
                      p0 = p0, x0 = x0, data = data, u = u, pu = pu,
                      control = list(fnscale = -5))
      
      output$value
    })
    p.mle <- p[which.max(l)]
    l.max <- max(l)
    # critical region of the test statistic
    q <- l.max - 0.5 * qchisq(alpha,1)
    p.lb <- l[l>= q][1]
    p.ub <- l[l <= q & p>p.mle][1]
    ci.lb <- p[match(p.lb,l)]
    ci.ub <- p[match(p.ub,l)]
    if (!any(is.na(c(ci.lb,ci.ub))) || iter >= max.iter) {
      break  # Exit loop if no NA or max iterations reached
    }
    if (is.na(ci.lb) || ci.lb==p.guess[1]){
      p.guess[1] <- p.guess[1] / expand.factor
    }
    if (is.na(ci.ub) || ci.ub==p.guess[2]){
      p.guess[2] <- p.guess[2] * expand.factor
    }
  
    iter <- iter + 1  
  }
  
  
  if (if.plot){
    plot(p,  l, type = "l", xlab = "p", ylab = 
           "Profile Log-likelihood")
    abline(h = q, col = "red")
    abline(v = p.mle, col = 4)
    abline(v = ci.lb, col = "red")
    abline(v = ci.ub, col = "red")
  }
  return(c(ci.lb,p.mle,ci.ub))
}

CI.Bayes <- function(x0,GP,alpha,if.plot=FALSE){
  # GP is a fevd object with method = Bayesian
  pp <- Tail.prob.GP(x0,GP)
  pp.density <- density(pp)
  MAP <- pp.density$x[which.max(pp.density$y)]
  MP <- quantile(pp,0.5)
  BCI <- quantile(pp,probs=c(1-alpha,alpha))
  if (if.plot){
    plot(pp.density,main=sprintf('Posterior distribution of P(X>%f)',x0),xlab='p')
    abline(v = MAP, col = 'red')
    abline(v = MP, col = 'green')
    abline(v = BCI[1], col = "blue")
    abline(v = BCI[2], col = "blue")
  }
  return(c(BCI[1],MP,MAP,BCI[2]))
}
