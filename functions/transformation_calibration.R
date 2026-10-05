# Calibration of transformation parameters


## transformation
SRP_transform <- function(X,alpha,delta){return(1/(X+delta)^alpha)} # shifted reciprocal power transformation
SRP_inverse <- function(Y, alpha, delta) {return(Y^(-1 / alpha) - delta)}

GE_transform <- function(X,theta,kappa){return(kappa*exp( -X^theta/theta) )} # generalized expoential link
GE_inverse <- function(Y, theta, kappa) {return((-theta * log(Y / kappa))^(1/theta))}


trans_quickTest <- function(data, type, obsCI, expo, prange1, prange2,u.prob,use_phi=FALSE) {
  # quickly check if the estimated frequency lies within the observed confidence interval for a grid of transformation parametersf 
  logical_matrix <- matrix(FALSE, nrow = length(prange1), ncol = length(prange2))
  freq_matrix <- matrix(NA, nrow = length(prange1), ncol = length(prange2))
  Qc_matrix <- matrix(NA, nrow = length(prange1), ncol = length(prange2))
  shape_matrix<- matrix(NA, nrow = length(prange1), ncol = length(prange2))
  
  # Qp_matrix <- matrix(NA, nrow = length(prange1), ncol = length(prange2))
  
  for (i in seq_along(prange1)) {
    for (j in seq_along(prange2)) {
      p1 <- prange1[i]
      p2 <- prange2[j]
      
      tdat <- switch(type,SRP = SRP_transform(data, p1, p2),
                     GE = GE_transform(data, p1, p2),
                     stop("Invalid transformation type"))
      x0 <- switch(type,SRP = SRP_transform(0,p1,p2),GE = GE_transform(0,p1,p2))
      u <- quantile(tdat,probs=u.prob)
      M <- fevd(x=tdat,threshold = u,type = 'GP',use.phi = use_phi)
      p0 <- tail_prob_gp(x0,model = M,use.phi=use_phi)
      pred <- p0 * expo * 365
      
      # Check if estimated frequency lies within the observed confidence interval
      in_obs <- pred >= obsCI[1] & pred <= obsCI[2]
      # Store results
      logical_matrix[i, j] <- in_obs
      freq_matrix[i, j] <- pred
      Qc_matrix[i, j] <- Qc(x=x0,data=tdat,u=u)
      shape_matrix[i, j] <- M$result$par[2]
      # Qp_matrix[i, j] <- Qp(tdat,p = u.prob)
    }
  }
  
  # Assign row and column names for clarity
  rownames(logical_matrix) <- prange1
  colnames(logical_matrix) <- prange2
  rownames(freq_matrix) <- prange1
  colnames(freq_matrix) <- prange2
  rownames(Qc_matrix) <- prange1
  colnames(Qc_matrix) <- prange2
  rownames(shape_matrix) <- prange1
  colnames(shape_matrix) <- prange2
  # rownames(Qp_matrix) <- prange1
  # colnames(Qp_matrix) <- prange2
  # Return the results as a list
  return(list(logical_matrix = logical_matrix, freq_matrix = freq_matrix,Qc_matrix=Qc_matrix,shape_matrix=shape_matrix))
}


freq_loss_fun_trans <- function(param,dat,u.prob,obs,type,weight=c(1,5,4),fixed.param=NULL, fixed.which=NULL,
                     expo,expand_factor,if.return=FALSE,use_phi=FALSE){
  # loss function based on the difference between estimated and observed frequencyhe
  # obs = observed CI of frequency
  # expo is exposure = length(data)/duration
  # expand_factor is used in tunning p.guess for CI.proLik
  # u is the threshold of excess prob.
  # fixed.param: fixed value for one param;
  type <- match.arg(type,c("SRP","GE"))
  if (!is.null(fixed.param) && !is.null(fixed.which)) {
    if (fixed.which == 1) {
      p1 <- fixed.param
      p2 <- param
    } 
    else if (fixed.which == 2) {
      p1 <- param
      p2 <- fixed.param
    } 
    else {
      stop("exceeds the length of parameters")
    }
  } 
  else {
    p1 <- param[1]
    p2 <- param[2]
  }
  tdat <- switch(type,SRP = SRP_transform(dat,p1,p2),GE = GE_transform(dat,p1,p2))
  x0 <- switch(type,SRP = SRP_transform(0,p1,p2),GE = GE_transform(0,p1,p2))
  # if (any(is.na(tdat)) || any(is.nan(tdat))) {
  #   print(list(p1 = p1, p2 = p2, tdat = head(tdat, 5)))  # Show the first few elements
  #   stop("NA/NaN detected after transformation.")
  # }
  u <- quantile(tdat,probs=u.prob)
  M <- fevd(x=tdat,threshold = u,type = 'GP',use.phi = use_phi)
  p0 <- tail_prob_gp(x0,model = M,use.phi = use_phi)
  
  # pred <- ci_prolik(x0,GP=M,expand_factor=expand_factor,use_phi = use_phi)
  pred <- ci_sim(x0,GP=M)
  pred <- pred * expo * 365
  total.loss <- weight[1]*(pred[1]-obs[1])+weight[2]*(obs[1]<pred[2] && pred[2]<obs[2])-
    weight[3]*(pred[2]-pred[1])
  
  if (if.return){
    loss.Obj <- list(value=total.loss,GP.model = M,PredCI = pred)
    return(loss.Obj)
  }
  return(total.loss)
}

transParam_select <- function(data,type,prange1,prange2,obsCI,u.prob,expo,expand_factor=10,
                              fixed.param=NULL, fixed.which=NULL,optim.method="L-BFGS-B",
                              user.control,use_phi=FALSE){
  # select trans param that minimizes the estimated Prob CI and the observedCI
  # data: raw data
  # type: "SRP" or "GE"; prange1,2 are the range of trans parameter, order matters, nint is the number of tries within the prange
  # obsCI, u.prob, expo and expand_factor pass down to freq_loss_fun_trans 
  
  if (!is.null(fixed.param) && !is.null(fixed.which)) {
    if (fixed.which == 1) {
      p2.start <- runif(1, min = prange2[1], max = prange2[2])
      start.param <- p2.start
      lower <- prange2[1]
      upper <- prange2[2]
    } 
    else if (fixed.which == "2") {
      p1.start <- runif(1, min = prange1[1], max = prange1[2])
      start.param <- p1.start
      lower <- prange1[1]
      upper <- prange1[2]
    } 
    else {
      stop("exceeds the length of parameters")
    }
  } 
  else {
    # Optimize both parameters if none are fixed
    p1.start <- runif(1, min = prange1[1], max = prange1[2])
    p2.start <- runif(1, min = prange2[1], max = prange2[2])
    start.param <- c(p1.start, p2.start)
    lower <- c(prange1[1], prange2[1])
    upper <- c(prange1[2], prange2[2])
  }
  default_control <- list(maxit = 100, fnscale = -5)
  
  # Merge default and user-provided controls
  final_control <- modifyList(default_control, user.control)
  output <- optim(start.param,
                  fn = freq_loss_fun_trans,
                  dat = data,
                  u.prob = u.prob,
                  obs = obsCI,
                  type = type,
                  expo = expo,
                  expand_factor = expand_factor,
                  fixed.param = fixed.param,
                  fixed.which = fixed.which,
                  use_phi = use_phi,
                  method = "L-BFGS-B",
                  lower = lower,
                  upper = upper,
                  control = final_control)
  
  # minimize the loss function first w.r.t param1, then param2
  # if (optim.method=="L-BFGS-S"){
  #   output <- optim(start.param,fn = freq_loss_fun_trans,
  #         dat = data,u.prob = u.prob,obs = obsCI,type = type,
  #         expo = expo,expand_factor = 10,
  #         fixed.param = fixed.param,fixed.which = fixed.which,
  #         method = optim.method,
  #         lower = lower,
  #         upper = upper,
  #         control = list(maxit = 100,fnscale=-5))
  # }
  # else{
  #   output<-optim(start.param,fn = freq_loss_fun_trans,
  #     dat = data,u.prob = u.prob,obs = obsCI,type = type,
  #     expo = expo,expand_factor = 10,
  #     fixed.param = fixed.param,fixed.which = fixed.which,
  #     method = optim.method,control = list(maxit = 100,fnscale=-5)) }
  return(output)
}

Qc <- function(x,data,u){return((x-max(data))/(max(data)-u))}
Qp <- function(data,p){
  pp <- quantile(data,probs = c(1-p,1-2*p,1-4*p))
  return( (pp[1]-pp[2])/(pp[2]-pp[3]))
}

