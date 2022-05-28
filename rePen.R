library(quantreg)

check <- function(x,tau=.5){
  x*(tau - (x<0))
}

pos_part <- function(x){
  ifelse( x < 0, x, 0 ) # min(x,0) # 
}

lasso <- function(x,lambda=1){
  lambda*abs(x)
}

scad <- function(x, lambda=1, a=3.7){
  absx <- abs(x)
  ifelse( absx < lambda, # Evaluate this
          lambda*absx, # If true (absx < lambda)
          ifelse( absx < a*lambda, # If false, evaluate this
                  ( (a^2-1)*lambda^2 - (absx-a*lambda)^2 ) / ( 2*(a-1) ), # If true (absx < a*lambda)
                  (a+1)*lambda^2 / 2 # If false (absx > a*lambda)
          ) 
  )
}


scad_deriv <- function(x, lambda=1,a=3.7){
  absx <- u <- abs(x)
  u[] <- 0
  index <- absx < a*lambda & absx > 0
  u[ index ] <-
    ifelse( absx[ index ] <= lambda, 
            lambda,
            ( a*lambda - absx[ index ] )/( a-1 ) 
    )
  u[index] <- u[index]*sign( x[index] )
  u[ x == 0 ] <- lambda # because we take derivative as x approaces 0 from above
  
  u
}

#scad_1_deriv <- function(x,lambda=1,a=3.7){
#  sign(x)*lambda
#}
#
#scad_2_deriv <- function(x,lambda=1,a=3.7){
#  sign(x)*lambda*(1-( pos_part(a*lambda-abs(x)) / ( lambda*(a-1))))*(abs(x) > lambda)
#}

mcp <- function(x, lambda=1, a=3){
  absx <- abs(x)
  ifelse( absx < a*lambda, # Evaluate this
          lambda*(absx - absx^2/(2*a*lambda)), # If true (absx < a*lambda)
          a*lambda^2/2 # If false
  )
}


mcp_deriv <- function(x, lambda=1, a=3){
  u <- x
  u[] <- 0
  index <- abs(x) < a*lambda
  u[ index ] <- ifelse( x[index] == 0,
                        lambda,
                        lambda*sign(x[index]) - x[index]/a
  )
  
  u
}

# MCL penalty
mcl <- function(x, lambda =1, gamma = 0.5, a = 3){
  absx <- abs(x)
  ifelse( absx < a*(lambda - gamma),
          lambda*absx - absx^2/(2*a),
          gamma*absx + a*(lambda^2 - gamma^2)/2
  )
}

# MCL derivation
mcl_deriv <- function(x, lambda =1, gamma = .5, a = 3){
  u <- x
  u[] <- lambda
  index <- abs(x) > 0
  u[index] <- ifelse(abs(x[index]) < a*(lambda - gamma),
                     lambda*sign(x[index]) - x[index]/a,
                     gamma*sign(x[index])
  )
  u
}


square <- function(x){
  x^2
}

randomly_assign <- function(n,k){
  #randomly assign n samples into k groups
  small_set <- floor(n/k)
  group_assign <- NULL
  if(n %% k == 0){
    group_assign <-  rep(seq(1,k),n/k)
  } else{
    remainder <- n %% k
    for(i in 1:remainder){
      group_assign <- c(group_assign, rep(i,small_set+1))
    }
    group_assign <- c(group_assign, rep(seq((i+1),k),small_set))
  }
  sample(group_assign)
}

rq.lasso.fit.mult <- function(x,y,tau_seq=c(.1,.3,.5,.7,.9),lambda=NULL,weights=NULL,intercept=TRUE,coef.cutoff=.00000001,...){
  model_list <- list()
  iter <- 1
  for(tau in tau_seq){
    model_list[[iter]] <- rq.lasso.fit(x,y,tau,lambda,weights,intercept,coef.cutoff,...)
    iter <- iter+1
  }
  model_list
}

transform_coefs <- function(coefs,mu_x,sigma_x,intercept=TRUE){
  new_coefs <- coefs
  if(intercept){
    intercept <- coefs[1]
    for(j in 2:length(coefs)){
      new_coefs[j] <- coefs[j]/sigma_x[j-1]
      intercept <- intercept - coefs[j]*mu_x[j-1]/sigma_x[j-1]
    }
    new_coefs[1] <- intercept
  } else{
    for(j in 1:length(coefs)){
      new_coefs[j] <- coefs[j]/sigma_x[j]
    }
  }
  new_coefs
}

get_coef_pen <- function(coefs,lambda,intercept,penVars,penalty="LASSO"){
  if(intercept){
    coefs <- coefs[-1]
  }
  if(is.null(penVars)==FALSE){
    coefs <- coefs[penVars]
  }
  if(penalty=="LASSO"){
    sum(abs(coefs))*lambda
  }
}

rq.lasso.fit <- function(x,y,tau=.5,lambda=NULL,weights=NULL,intercept=TRUE,
                         coef.cutoff=1e-08, method="br",penVars=NULL,scalex=TRUE, ...){
  # x is a n x p matrix without the intercept term
  # y is a n x 1 vector
  # lambda takes values of 1 or p
  # coef.cutoff is a threshold to set to zero. 
  # Choose the method used to estimate the coefficients ("br" or "fn")
  ### According to quantreg manual and my experience, "fn" is much faster for big n
  ### The "n" can grow rapidly using lin. prog. approach  
  # penVars - variables to be penalized, doesn't work if lambda has multiple entries (Ben: I think it does though it is a little bit strange to do)
  if(is.null(dim(x))){
    stop('x needs to be a matrix with more than 1 column')
  }
  p <- dim(x)[2]
  if(p == 1){
    stop('x needs to be a matrix with more than 1 column')
  }
  n <- dim(x)[1]
  if(n != length(y)){
    stop('length of y and rows of x do not match')
  }
  if(is.null(lambda)==TRUE | (length(lambda) != 1 & length(lambda) != dim(x)[2])){
    stop(paste('input of lambda must be of length 1 or', dim(x)[2]))
  }
  if( sum(lambda < 0) > 0){
    stop(paste('lambda must be positive and we have a lambda of ', lambda, sep=""))
  }
  if(scalex){
    original_x <- x
    x <- scale(x)
    mu_x <- attributes(x)$`scaled:center`
    sigma_x <- attributes(x)$`scaled:scale`
  }
  #if(is.null(method)){
  # if(n + 2*p < 500){
  #    method <- "br"
  # } else{
  #	method <- "fn"
  # }
  #}
  
  # ##############################################################################################
  # ### This uses the LASSO.fit or LASSO.fit.nonpen functions to obtain coefficient estimates. ###
  # ### Note that the coefficients might need to be reordered to match x.                      ###
  # ##############################################################################################
  # if( !is.null(weights) & length(weights) != n )
  #    stop("Length of weights does not match length of y")
  
  # ### Find indices of penalized and nonpenalized coefficients
  # nonpenVars <-     ### indices of nonpenalized coefficients (lambdas of 0 or not included in penVars), NULL if no nonpenalized oefficients
  # penVars    <-     ### indices of penalized coefficients
  
  # xnew <- as.matrix( x[,penVars] )
  # if( is.null(nonpenVars) ){
  #   coefs <- LASSO.fit(y, xnew, tau, intercept, coef.cutoff, weights)
  # } else {
  #   znew <- as.matrix( x[,nonpenVars] )
  #   coefs <- LASSO.fit.nonpen(y, xnew, znew, tau, intercept, coef.cutoff, weights)
  #   coefs[ intercept + 1:p ] <- coefs[ intercept + order(c(penVars, nonpenVars)) ]
  # }
  
  # ### coefs are the coefficients in the correct order with the intercept first
  # ##############################################################################################
  # ##############################################################################################
  # ##############################################################################################
  
  if(is.null(penVars) !=TRUE){# & length(lambda) == 1){
    if(length(lambda)==1){
      mult_lambda <- rep(0,p)
      mult_lambda[penVars] <- lambda
      lambda <- mult_lambda
    } else{
      lambda[-penVars] <- 0
    }
  }
  lambda <- lambda*n # need this to account for the fact that rq does not normalize the objective function
  if(length(lambda)==1){
    pen_x <- rbind(diag(rep(lambda,p)),diag(rep(-lambda,p)))
  } else{
    pen_x <- rbind(diag(lambda), diag(-lambda))
    pen_x <- pen_x[rowSums(pen_x==0)!=dim(pen_x)[2],]#drop all zero rows
  }
  aug_n <- dim(pen_x)[1]
  aug_x <- rbind(x,pen_x)
  if(intercept){
    aug_x <- cbind(c(rep(1,n),rep(0,aug_n)), aug_x)
  }
  aug_y <- c(y, rep(0,aug_n))
  if(is.null(weights)){
    model <- rq(aug_y ~ aug_x+0, tau=tau, method=method)
  } else{
    if(length(weights) != n){
      stop("Length of weights does not match length of y")
    }
    orig_weights <- weights
    weights <- c(weights, rep(1,aug_n))
    model <- rq(aug_y ~ aug_x+0, tau=tau, weights=weights, method=method)
  }
  p_star <- p+intercept
  coefs <- coefficients(model)[1:p_star]
  return_val <- NULL
  return_val$coefficients <- coefs
  if(is.null(colnames(x))){
    x_names <- paste("x",1:p,sep="")
  } else{
    x_names <- colnames(x)
  }
  if(intercept){
    x_names <- c("intercept",x_names)
  }
  attributes(return_val$coefficients)$names <- x_names
  return_val$coefficients[abs(return_val$coefficients) < coef.cutoff] <- 0
  if(scalex){
    #need to update for penVars
    return_val$coefficients <- transform_coefs(return_val$coefficients,mu_x,sigma_x,intercept)
    if(intercept){
      fits <- cbind(1,original_x) %*% return_val$coefficients
    } else{
      fits <- original_x %*% return_val$coefficients
    }
    return_val$residuals <- y - fits
    return_val$PenRho <- sum(sapply(return_val$residuals,check,tau))+get_coef_pen(return_val$coefficients,lambda,intercept,penVars)	 
  } else{
    return_val$PenRho <- model$rho
    return_val$residuals <- model$residuals[1:n]   
  }
  if(is.null(weights)){   
    return_val$rho <- sum(sapply(return_val$residuals,check,tau))
  } else{
    return_val$rho <- sum(orig_weights*sapply(return_val$residuals,check,tau))
  }
  return_val$tau <- tau
  return_val$n <- n                  
  return_val$intercept <- intercept
  class(return_val) <- c("rq.pen", "rqLASSO")
  return_val
}

predict.rq.pen <- function(object, newx,...){
  coefs <- object$coefficients
  if(object$intercept){
    newx <- cbind(1,newx)
  }
  newx %*% coefs
}

predict.cv.rq.pen <- function(object, newx, lambda="lambda.min",...){
  if(lambda == "lambda.min"){
    target_pos <- which(object$cv$lambda == object$lambda.min)
  } else{
    target_pos <- which(object$cv$lambda == lambda)
  }
  predict(object$models[[target_pos]],newx,...)
}

getRho <- function(model){
  model$rho
}

coef.cv.rq.group.pen <- function(object, lambda='min',...){
  if(lambda=='min'){
    lambda <- object$lambda.min
  }
  target_model <- which(object$cv[,1] == lambda)
  object$beta[,target_model]
}

group_derivs <- function(deriv_func,groups,coefs,lambda,a=3.7){
  if(length(lambda)==1){
    lambda <- rep(lambda,length(groups))
  }
  derivs <- NULL
  for(g in 1:length(unique(groups))){
    g_index <- which(groups==g)
    current_lambda <- lambda[g]
    coefs_l1 <- sum(abs(coefs[g_index]))
    derivs <- c(derivs, deriv_func(coefs_l1,current_lambda,a))
  }
  derivs
}

rq.group.lin.prog <- function(x,y,groups,tau,lambda,intercept=TRUE,eps=1e-05,penalty="SCAD", a=3.7, coef.cutoff=1e-08,
                              initial_beta=NULL,iterations=10,converge_criteria=.0001,penGroups=NULL,...){
  group_num <- length(unique(groups))
  if(length(lambda) == 1){
    lambda <- rep(lambda,group_num)
  }
  if (length(lambda) != group_num) {
    stop("lambdas do not match with group number")
  }
  if (sum(groups == 0) > 0) {
    stop("0 cannot be used as a group")
  }
  if (dim(x)[2] != length(groups)) {
    stop("length of groups must be equal to number of columns in x")
  }
  if (penalty == "SCAD") {
    deriv_func <- scad_deriv
  }
  if (penalty == "MCP") {
    deriv_func <- mcp_deriv
  } 
  if (penalty == "MCL") {
    deriv_func <- mcl_deriv
  }  
  new_lambda <- NULL
  group_count <- xtabs(~groups)
  for (g in 1:group_num) {
    new_lambda <- c(new_lambda, rep(lambda[g], each = group_count[g]))
  }
  if(is.null(penGroups)==FALSE){
    zero_pen_spots <- which(!groups %in% penGroups)
    new_lambda[zero_pen_spots] <- 0
  }
  if(is.null(initial_beta)){
    initial_beta <- rq.lasso.fit(x,y,tau,new_lambda, intercept=intercept, coef.cutoff=coef.cutoff,...)$coefficients
  }
  
  coef_by_group_deriv <- group_derivs(deriv_func, groups, initial_beta,lambda,a)
  lambda_update <- coef_by_group_deriv[groups]
  old_beta <- initial_beta
  
  iter_complete <- FALSE
  iter_num <- 0
  
  #pen_range <- (1+intercept):(dim(x)[2]+intercept)
  coef_range <- (1+intercept):(dim(x)[2]+intercept)
  
  while(!iter_complete){
    sub_fit <- rq.lasso.fit(x=x,y=y,tau=tau,lambda=lambda_update,intercept=intercept,...)
    coef_by_group_deriv <- group_derivs(deriv_func,groups,sub_fit$coefficients[coef_range],lambda,a)
    lambda_update <- coef_by_group_deriv[groups]
    if(is.null(penGroups)==FALSE){
      lambda_update[zero_pen_spots] <- 0
    }
    iter_num <- 1
    new_beta <- sub_fit$coefficients
    beta_diff <- sum( (old_beta - new_beta)^2)
    if(iter_num == iterations | beta_diff < converge_criteria){
      iter_complete <- TRUE
      if(iter_num == iterations & beta_diff > converge_criteria){
        warning(paste("did not converge after ", iterations, " iterations", sep=""))
      }
    } else{
      old_beta <- new_beta
    }
  }
  sub_fit$penalty <- penalty
  class(sub_fit) <-  c("rq.pen", "rqNC")
  sub_fit
}



groupMultLambda <- function (x, y, groups, tau = 0.5, lambda, intercept = TRUE, penalty="LASSO", 
                             #initial_beta = NULL,
                             alg="QICD_warm",penGroups=NULL, ...) 
{
  if(alg != "QICD_warm"){
    #don't know how to do warm start with linear programming approach 
    return_val <- list()
    pos <- 1
    for (lam in lambda) {
      return_val[[pos]] <- rq.group.fit(x = x, y = y, groups = groups, 
                                        tau = tau, lambda = lam, intercept = intercept, penalty=penalty,alg=alg, penGroups=penGroups,
                                        ...)
      #initial_beta <- return_val[[pos]]$coefficients
      pos <- pos + 1
    }
  } else{
    p <- dim(x)[2]
    pos <- 1
    alg = "QICD"
    
    return_val <- list()
    if(intercept){
      initial_beta <- c(quantile(y,tau), rep(0,p))
    } else{
      initial_beta <- rep(0,p)
    }
    
    for(lam in lambda){
      return_val[[pos]] <- rq.group.fit(x=x, y=y, groups=groups, tau=tau, lambda= lam, intercept=intercept, penalty="LASSO", alg=alg, initial_beta=initial_beta, penGroups=penGroups, ...)
      initial_beta <- coefficients(return_val[[pos]])
      pos <- pos + 1
    }
    
    #if penalty is not lasso then update those initial estimates
    if(penalty != "LASSO"){
      pos <- 1
      for(lam in lambda){
        initial_beta <- coefficients(return_val[[pos]]) #use lasso estimate as initial estimate
        return_val[[pos]] <- rq.group.fit(x=x, y=y, groups=groups, tau=tau, lambda= lam, intercept=intercept, penalty=penalty, alg=alg, initial_beta=initial_beta, penGroups=penGroups, ...)
        pos <- pos + 1
      }
    }
  }
  return_val
}

nonzero <- function (obj) 
{
  UseMethod("nonzero")
}

nonzero.cv.rq.group.pen <- function (obj) 
{
  coefs <- coefficients(obj)
  if (obj$intercept) {
    coefs <- coefs[-1]
  }
  tapply(coefs, obj$groups, sum) != 0
}

print.cv.rq.pen <- function(x,...){
  cat("\nCoefficients:\n")
  print(coefficients(x,...))
  cat("\nCross Validation (or BIC) Results\n")
  print(x$cv)
}

print.rq.pen <- function(x,...){
  cat("\nCoefficients:\n")
  print(coefficients(x,...))
}


### This is a stripped down version of rq that
### only estimates betas without bells and whistles uses "br" method
#### Useful for problems with sample size up to several thousand
shortrq.fit.br <- function (x, y, tau = 0.5)
{
  tol <- .Machine$double.eps^(2/3)
  eps <- tol
  big <- .Machine$double.xmax
  x <- as.matrix(x)
  p <- ncol(x)
  n <- nrow(x)
  ny <- NCOL(y)
  nsol <- 2
  ndsol <- 2
  # # Check for Singularity of X since br fortran isn't very reliable about this
  # if (qr(x)$rank < p)
  #     stop("Singular design matrix")
  lci1 <- FALSE
  qn <- rep(0, p)
  cutoff <- 0
  
  z <- .Fortran("rqbr", as.integer(n), as.integer(p), as.integer(n +
                                                                   5), as.integer(p + 3), as.integer(p + 4), as.double(x),
                as.double(y), as.double(tau), as.double(tol), flag = as.integer(1),
                coef = double(p), resid = double(n), integer(n), double((n +
                                                                           5) * (p + 4)), double(n), as.integer(nsol), as.integer(ndsol),
                sol = double((p + 3) * nsol), dsol = double(n * ndsol),
                lsol = as.integer(0), h = integer(p * nsol), qn = as.double(qn),
                cutoff = as.double(cutoff), ci = double(4 * p), tnmat = double(4 *
                                                                                 p), as.double(big), as.logical(lci1))
  # if (z$flag != 0)
  #     warning(switch(z$flag, "Solution may be nonunique", "Premature end - possible conditioning problem in x"))
  coef <- z$coef
  coef
}

### This is a stripped down version of rq that
### only estimates betas without bells and whistles uses "fn" method
#### Useful for problems with huge sample size
shortrq.fit.fnb <- function (x, y, tau = 0.5, beta = 0.99995, eps = 1e-06)
{
  n <- length(y)
  p <- ncol(x)
  rhs <- (1 - tau) * apply(x, 2, sum)
  d   <- rep(1,n)
  u   <- rep(1,n)
  wn <- rep(0,10*n)
  wn[1:n] <- (1-tau) #initial value of dual solution
  z <- .Fortran("rqfnb", as.integer(n), as.integer(p), a = as.double(t(as.matrix(x))),
                c = as.double(-y), rhs = as.double(rhs), d = as.double(d),as.double(u),
                beta = as.double(beta), eps = as.double(eps),
                wn = as.double(wn), wp = double((p + 3) * p),
                it.count = integer(3), info = integer(1))
  if (z$info != 0)
    stop(paste("Error info = ", z$info, "in stepy: singular design"))
  coefficients <- -z$wp[1:p]
  names(coefficients) <- dimnames(x)[[2]]
  return ( coefficients )
}




QICD <- function(y, x, tau=.5, lambda, gamma, intercept=TRUE, penalty="SCAD", 
                 initial_beta=NULL, maxin=100, maxout=20, eps = 1e-05, coef.cutoff=1e-08,  
                 a=3.7, scalex=TRUE, ...)
  #y: response variable, length n vector
  #x: input nxp matrix, of dimension nobs x nvars; each row is an observation vector. 
  #tau is the quantile value
  #lambda is the tuning parameter (numeric value > 0)
  #intercept is a logical value,should intercept be fitted (default=TRUE) or set to zero (FALSE)
  #maxin: maximum number of iterations for inside coordinate descent,default value is 100
  #maxout: maximum number of iterations for outside MM step, default value is 20
  #initial_beta: initial value for x-covariates, the default value is NULL (lasso estimates will be used)
  #eps is the convergence threshold for coordinate descent and majorization minimization step
  #penalty is the name of the penalty function ("SCAD", "MCP", "LASSO")
  #a is scale parameter, the default value is 3.7 (>2 for SCAD, >1 for MCP, not used in LASSO)
#coef.cutoff is threshold for determining nonzero coefficients
#initial_beta is vector containing initial values for intercept (if included) and coefficients.
### Should be in the form (intercept, coefficients) intercept should be left out if intercept=FALSE
{
  cleanInputs(y, x, lambda, gamma, initial_beta, intercept, penalty, a)
  
  if(scalex){
    x <- scale(x)
    mu_x <- attributes(x)$`scaled:center`
    sigma_x <- attributes(x)$`scaled:scale`
  }
  
  # Set penalty function
  if(penalty == "SCAD"){
    pentype <- as.integer(0)  
  } else if (penalty == "MCP"){
    pentype <- as.integer(1)
  } else if (penalty == "MCL"){
    pentype <- as.integer(3) 
  } else{
    pentype <- as.integer(2)
  }
  
  if( is.null(initial_beta) ){
    # initial_beta <- LASSO.fit(y, x, tau, lambda, intercept, coef.cutoff)
    initial_beta <- coefficients( cv.rq.pen(x, y, tau=tau, intercept=intercept, 
                                            penalty="LASSO", criteria="BIC") )
  }
  
  if( intercept ){
    beta <- initial_beta[-1]
    intval <- initial_beta[1]
  } else{
    beta <- initial_beta
    intval <- 0    
  }
  
  
  n         <- length(y)
  y         <- as.double(y)
  xdoub     <- as.double(x)
  p         <- as.integer( ncol(x) )
  tau       <- as.double(tau)
  int       <- as.integer(intercept)
  nint      <- as.integer( length(y) )
  a         <- as.double(a)
  eps    <- as.double(eps)
  maxin     <- as.integer(maxin)
  
  if( length(lambda) == 1 ){
    lambda <- rep( lambda, p)
  } else if ( length(lambda) != p ){
    stop("lambda must be of length 1 or p")
  }
  lambda    <- as.double(lambda)
  
  if( length(gamma) == 1 ){
    gamma <- rep( gamma, p)
  } else if ( length(gamma) != p ){
    stop("gamma must be of length 1 or p")
  }
  gamma    <- as.double(gamma)
    
  i=0
  distance <- eps+1
  groupl1 <- rep(0, p)
  beta0 <- beta
  intval0 <- intval
  residuals <- as.double(y - x%*%beta - intval)
  
  while( (i < maxout) & (distance >= eps) ){
    
    #lambda <- lambda*n
    out <- .C("penderiv", as.double(beta), p, a, lambda, gamma, pentype)
    #penweight <- as.double(out[[1]])
    penweight <- as.double( n*out[[1]] )
    
    out <- .C("QCD", xdoub, as.double(beta), as.double(intval), penweight, residuals,
              nint, p, int, tau, eps, maxin)
    beta <- out[[2]]
    intval <- out[[3]]
    residuals <- as.double( out[[5]] )
    i <- i+1
    distance <- sqrt( sum((beta - beta0)^2) + (intval0 - intval)^2 )
    beta0 <- beta
    intval0 <- intval
  }
  
  if(i == maxout & distance > eps){
    warning(paste("did not converge after ", maxout, " iterations", sep=""))
  }    
  
  beta[ abs(beta) < coef.cutoff ] <- 0
  coefficients <- beta
  if(intercept){
    coefficients <- c(intval, beta)
  }
  if(scalex){
    coefficients <- transform_coefs(coefficients,mu_x,sigma_x,intercept)
  }
  return( coefficients )
}



QICD.nonpen <- function(y, x, z, tau=.5, lambda, intercept=TRUE, penalty="SCAD", 
                        initial_beta=NULL, maxin=100, maxout=20, eps = 1e-05, coef.cutoff=1e-08,  
                        a=3.7, method="br", scalex=TRUE, ...)
  #y: response variable, length n vector
  #x: input nxp matrix, of dimension nobs x nvars; each row is an observation vector.
  #z is nxq matrix of bases; the coefficients for these columns will be unpenalized
  #tau is the quantile value
  #lambda is the tuning parameter (numeric value > 0)
  #intercept is a logical value,should intercept be fitted (default=TRUE) (intercept should be included when using splines)
  #maxin: maximum number of iterations for inside coordinate descent,default value is 100
  #maxout: maximum number of iterations for outside MM step, default value is 20
  #eps is the convergence threshold for coordinate descent and majorization minimization step
  #penalty is the name of the penalty function ("SCAD", "MCP", "LASSO")
  #a is scale parameter, the default value is 3.7 (>2 for SCAD, >1 for MCP, not used in LASSO)
#coef.cutoff is threshold for determining nonzero coefficients
#initial_beta is vector containing initial values for intercept (if included) and x coefficients.
### Should be in the form (intercept, coefficients) intercept should be left out if intercept=FALSE.
### The intercept should be included to be consistent with other methods, but intercept and z coefficients
### will be initialized to rq( y-x%*%beta ~ z ).
#method for quantile regression: can be "br" or "fn", see top for description.
{
  cleanInputs(y, x, lambda, initial_beta[ 1:(intercept+ncol(x)) ], intercept, penalty, a)
  if(scalex){
    x <- scale(x)
    mu_x <- attributes(x)$`scaled:center`
    sigma_x <- attributes(x)$`scaled:scale`
  }
  
  if(is(z,"matrix") == FALSE ){                                                                                    
    stop('z needs to be a matrix')
  }
  
  if( nrow(z) != length(y) ){
    stop('length of y and rows of z do not match')
  }
  
  if( method == "br"){
    zreg <- shortrq.fit.br
  } else if ( method == "fn" ){
    zreg <- shortrq.fit.fnb 
  } else {
    stop("Incorrect method.  Choose br or fn")
  }
  
  # Set penalty function
  if(penalty == "SCAD"){
    pentype <- as.integer(0)  
  } else if (penalty == "MCP"){
    pentype <- as.integer(1)
  } else if (penalty == "MCL"){
    pentype <- as.integer(3)
  } else{
    pentype <- as.integer(2)
  }
  
  zmat <- z
  if( intercept )
    zmat <- cbind(1,z)
  
  # Check for Singularity of X since br fortran isn't very reliable about this
  if (qr(zmat)$rank < ncol(zmat))
    stop("z is a singular matrix (make sure intercept column is not included)")  
  
  if( is.null(initial_beta) ){
    initial_beta <- LASSO.fit.nonpen(y, x, z, tau, lambda, intercept, coef.cutoff) 
    initial_beta <- initial_beta[ 1:(intercept+ncol(x)) ] ### Only keep the coefficients for x 
  } else {
    initial_beta <- initial_beta[ 1:(intercept+ncol(x)) ] ### Only keep the coefficients for x 
  }
  
  if( intercept ){
    beta <- initial_beta[-1]
    beta <- beta
  } else{
    beta <- initial_beta
  }
  
  
  
  n         <- length(y)
  xdoub     <- as.double(x)
  p         <- as.integer( ncol(x) )
  tau       <- as.double(tau)
  nint      <- as.integer( length(y) )
  a         <- as.double(a)
  eps       <- as.double(eps)
  maxin     <- as.integer(maxin)
  
  if( length(lambda) == 1 ){
    lambda <- rep( lambda, p)
  } else if ( length(lambda) != p ){
    stop("lambda must be of length 1 or p")
  }
  lambda    <- as.double(lambda)
  
  i=0
  distance <- distance.inner <-  eps+1
  beta0 <- beta00 <- beta
  
  xb <- x%*%beta
  xresiduals <- y - xb
  zbeta <- zreg(zmat, xresiduals, tau = tau)
  zbeta0 <- zbeta00 <- zbeta
  zb <- zmat%*%zbeta
  residuals <- as.double(y - xb - zb)
  
  while( (i < maxout) & (distance >= eps) ){
    
    ii=0
    out <- .C("penderiv", as.double(beta), p, a, lambda, pentype)
    penweight <- as.double( n*out[[1]] )
    
    while( (ii < maxin) & distance.inner >= eps ){
      
      
      out <- .C("QCD", xdoub, as.double(beta), as.double(0), penweight, residuals,
                nint, p, as.integer(0), tau, eps, as.integer(1))
      beta <- out[[2]]
      
      xb <- x%*%beta
      xresiduals <- y - xb
      zbeta <- zreg(zmat, xresiduals, tau = tau)
      zb <- zmat%*%zbeta
      residuals <- as.double(y - xb - zb)
      
      ii <- ii+1
      distance.inner <- sqrt( sum((beta - beta00)^2) + sum((zbeta - zbeta00)^2) )
      beta00 <- beta
      zbeta00 <- zbeta
    }
    
    i <- i+1
    distance <- sqrt( sum((beta - beta0)^2) + sum((zbeta - zbeta0)^2) )
    beta0 <- beta
    zbeta0 <- zbeta
  }
  
  if(i == maxout & distance > eps){
    warning(paste("did not converge after ", maxout, " iterations", sep=""))
  }    
  
  beta[ abs(beta) < coef.cutoff ] <- 0
  if(scalex){
    beta <- transform_coefs(beta,mu_x,sigma_x,intercept)
  }
  
  coefficients <- c(beta, zbeta)
  if(intercept){
    coefficients <- c(zbeta[1], beta, zbeta[-1])
  }
  
  return( coefficients )
}







QICD.group <- function(y, x, groups, tau=.5, lambda, intercept=TRUE, penalty="SCAD", 
                       initial_beta=NULL, maxin=100, maxout=20, eps = 1e-05, coef.cutoff=1e-08,  
                       a=3.7, scalex=TRUE, ...)
  #y: response variable, length n vector
  #x: input nxp matrix, of dimension nobs x nvars; each row is an observation vector. 
  #groups: numeric vector of length ncol(x) with the group number of the coefficient (can be unique)
  #tau is the quantile value
  #lambda is the tuning parameter (numeric value > 0)
  #intercept is a logical value,should intercept be fitted (default=TRUE) or set to zero (FALSE)
  #maxin: maximum number of iterations for inside coordinate descent,default value is 100
  #maxout: maximum number of iterations for outside MM step, default value is 20
  #initial_beta: initial value for x-covariates, the default value is NULL (lasso estimates will be used)
  #eps is the convergence threshold for coordinate descent and majorization minimization step
  #penalty is the name of the penalty function ("SCAD", "MCP", "LASSO")
#coef.cutoff is threshold for determining nonzero coefficients
#initial_beta is vector containing initial values for intercept (if included) and coefficients.
#a is scale parameter, the default value is 3.7 (>2 for SCAD, >1 for MCP, not used in LASSO)

### Should be in the form (intercept, coefficients) intercept should be left out if intercept=FALSE
{
  cleanInputs(y, x, lambda, initial_beta, intercept, penalty, a)
  
  if(scalex){
    x <- scale(x)
    mu_x <- attributes(x)$`scaled:center`
    sigma_x <- attributes(x)$`scaled:scale`
  }
  
  # Set penalty function
  if(penalty == "SCAD"){
    pentype <- as.integer(0)  
  } else if (penalty == "MCP"){
    pentype <- as.integer(1)
  } else if (penalty == "MCL"){
    pentype <- as.integer(3)
  } else{
    pentype <- as.integer(2)
  }
  
  if( is.null(initial_beta) ){
    initial_beta <- LASSO.fit(y, x, tau, lambda, intercept, coef.cutoff) 
  }
  
  if( intercept ){
    beta <- initial_beta[-1]
    intval <- initial_beta[1]
  } else{
    beta <- initial_beta
    intval <- 0    
  }
  
  
  n         <- length(y)
  y         <- as.double(y)
  xdoub     <- as.double(x)
  p         <- as.integer( ncol(x) )
  tau       <- as.double(tau)
  int       <- as.integer(intercept)
  nint      <- as.integer( length(y) )
  a         <- as.double(a)
  eps    <- as.double(eps)
  maxin     <- as.integer(maxin)
  
  if( length(lambda) == 1 ){
    lambda <- rep( lambda, p)
  } else if ( length(lambda) != p ){
    stop("lambda must be of length 1 or p")
  }
  lambda    <- as.double(lambda)
  
  i=0
  distance <- eps+1
  groupl1 <- rep(0, p)
  beta0 <- beta
  intval0 <- intval
  residuals <- as.double(y - x%*%beta - intval)
  
  badValues <- FALSE
  
  while( (i < maxout) & (distance >= eps) ){
    
    for(grps in unique(groups)){
      groupl1[groups==grps] <- sum( abs(beta[groups==grps]) )
    }
    #lambda <- n*lambda
    out <- .C("penderiv", as.double(groupl1), p, a, lambda, pentype)
    #penweight <- as.double(out[[1]])
    penweight <- as.double( n*out[[1]] )
    
    out <- .C("QCD", xdoub, as.double(beta), as.double(intval), penweight, residuals,
              nint, p, int, tau, eps, maxin)
    beta <- out[[2]]
    intval <- out[[3]]
    residuals <- as.double( out[[5]] )
    i <- i+1
    distance <- sqrt( sum((beta - beta0)^2) + (intval0 - intval)^2 )
    beta0 <- beta
    intval0 <- intval
    
    if( max(abs(beta)) > 1e10 ){
      badValues <- TRUE
      break
    }
  }
  
  if( badValues ){
    warning("Some coefficients diverged to infinity (bad results)")    
  }
  
  if(i == maxout & distance > eps){
    warning(paste("did not converge after ", maxout, " iterations", sep=""))
  }    
  
  beta[ abs(beta) < coef.cutoff ] <- 0
  coefficients <- beta
  if(intercept){
    coefficients <- c(intval, beta)
  }
  
  if(scalex){
    coefficients <- transform_coefs(coefficients,mu_x,sigma_x,intercept)
  }
  
  return( coefficients )
}


### Can fit models for many different lambdas using warm start.  
## If no intial_beta is provided, QICD.master will find an appropriate starting value and fit the model for each lambda
QICD.master <- function(y, x, z=NULL, groups=NULL, tau=.5, lambda, intercept=TRUE, penalty="SCAD", 
                        initial_beta, maxin=100, maxout=20, 
                        eps = 1e-05, coef.cutoff=1e-08, a=3.7, ...){
  if( !is.null(z) & !is.null(groups) )
    stop("Currently not supporting group penalty and some nonpenalized variables \n  z or groups (or both) must be NULL")
  
  lambda <- sort( unique(lambda) )
  if( lambda[1] <= 0)
    stop("lambda must be positive")
  nlam <- length(lambda)
  p <- ncol(x)
  q <- ifelse( is.null(z), 0, ncol(z) )
  penVars <- intercept + 1:p ### Keep track of penalized coefficients
  
  ### Get names of output correct
  out <- matrix(0, nrow=intercept+p+q, ncol=nlam) ### Each column is coefficients for 1 lambda
  rnms <- paste("x",1:p, sep="")
  if( intercept )
    rnms <- c("(Intercept)", rnms)
  if( !is.null(z) )
    rnms <- c(rnms, paste("z",1:q, sep=""))
  
  rownames(out) <- rnms
  colnames(out) <- 1:nlam
  
  ### Set correct QICD function
  QF <- match.call()
  if( is.null(z) & is.null(groups) ){ ### QICD
    QF[[1]] <- as.name("QICD")
  } else if (is.null(z) ){            ### QICD.group
    QF[[1]] <- as.name("QICD.group")
  } else {                            ### QICD.nonpen
    QF[[1]] <- as.name("QICD.nonpen")
  }
  
  #####################################
  ### Don't do warm start (for now) ###
  #####################################
  ### Do first iteration
  QF$lambda <- lambda[1]
  # if( penalty == "LASSO" & is.null(groups) & !is.null(initial_beta) ){     ### Use initial_beta as coefficients for smallest lambda when LASSO with no group penalty
  #     out[,1] <- initial_beta
  # } else if (penalty == "LASSO" & is.null(groups) & is.null(initial_beta)) {                                        ### SCAD, MCP, group LASSO penalty use initial_beta as starting values
  #     QF$maxout <- -1 
  #     out[,1] <- eval.parent(QF)
  # } else {
  #     out[,1] <- eval.parent(QF)
  # }
  
  # QF$maxout <- maxout
  #####################################
  out[,1] <- eval.parent(QF)
  
  ibeta <- out[,1]
  if( nlam > 1 ){
    for( i in 2:nlam ){
      
      if( max( abs(ibeta[penVars]) ) == 0 ){ ### If all penalized betas are 0, then we are done 
        out[,i:nlam] <- matrix(ibeta, nrow=length(ibeta), ncol=nlam-i+1) ### Assign coefficients for remaining lambdas to current value 
        break ### End the for loop
      }
      
      QF$lambda <- lambda[i]
      # QF$initial_beta <- ibeta # For now, use initial_beta as starting values for all lambda
      out[,i] <- eval.parent(QF)
      ibeta <- out[,i]
    }
  }
  
  return(list( coefficients=out, lambda=lambda ))
}



### LASSO estimates for initial values in QICD and QICD.group functions when initial_beta = NULL
LASSO.fit <- function(y, x, tau, lambda, intercept, coef.cutoff, weights=NULL)
{
  p <- ncol(x)
  n <- nrow(x)
  ynew <- c( y, rep(0, 2*p) )
  xnew <- rbind( x, diag(n*lambda, p), -diag(n*lambda, p) )
  if( intercept )
    xnew <- cbind( c( rep(1,n), rep(rep(0, 2*p)) ), xnew )
  
  if( !is.null(weights) ){
    xnew[1:n,] <- xnew[1:n,] * weights
    ynew[1:n]  <-  ynew[1:n] * weights
  }
  
  #Ben: had some problems with fnb, setting all to use br for now
  #if( n + 2*p < 500 ){ ### If problem is small, use "br" method
  out <- shortrq.fit.br(xnew, ynew, tau)
  #} else {             ### Else use "fn" method
  #  out <- shortrq.fit.fnb(xnew, ynew, tau)
  #}
  
  out[ abs(out) < coef.cutoff ] <- 0
  return(out)
}


### LASSO estimates for initial values in QICD.nonpen function when initial_beta = NULL
LASSO.fit.nonpen <- function(y, x, z, tau, lambda, intercept, coef.cutoff, weights=NULL)
{
  p <- ncol(x)
  pxz <- ncol(z) + p
  n <- nrow(x)
  ynew <- c( y, rep(0, 2*p) )
  xz <- cbind( x, z )
  xz <- rbind( xz, diag(n*lambda, p, pxz), -diag(n*lambda, p, pxz) )
  if( intercept )
    xz <- cbind( c( rep(1,n), rep(rep(0, 2*p)) ), xz )
  
  if( !is.null(weights) ){
    xz[1:n,]   <-   xz[1:n,] * weights
    ynew[1:n]  <-  ynew[1:n] * weights
  }
  
  #Ben: had some problems with fnb, setting all to use br for now
  #if( n + 2*p < 500 ){ ### If problem is small, use "br" method
  out <- shortrq.fit.br(xz, ynew, tau)
  #} else {             ### Else use "fn" method
  #  out <- shortrq.fit.fnb(xz, ynew, tau)
  #}
  
  out <- out
  
  out[ abs(out) < coef.cutoff ] <- 0
  return(out)
}







#################
### Cleaning up the functions
#################
### This function just checks to make sure that the inputs to the QICD functions are appropriate
# Easier to create this function to include in each QICD function
cleanInputs <- function(y, x, lambda, gamma, initial_beta=NULL, intercept=TRUE,
                        penalty, a, ...){
  if( is(x,"matrix") == FALSE){                                                                                    
    stop('x needs to be a matrix')
  }
  
  if( any(lambda <= 0)){
    stop("lambda must be positive")
  }

  if( any(gamma <= 0)){
    stop("gamma must be positive")
  }

  if( any(gamma > lambda)){
    stop("gamma must be smaller than lambda")
  }
  
  # Make sure we use br or fn method ("?rq.fit.br" or "?rq.fit.fnc")
  # if( method != "br" & method != "fn"){
  #   stop("Incorrect method.  Choose br or fn")
  # }
  
  if( nrow(x) != length(y) ){
    stop('length of y and rows of x do not match')
  }
  
  if( !is.null(initial_beta) & (length(initial_beta) < (ncol(x) + intercept)) ){
    stop("initial_beta must contain initial value for intercept (if TRUE) and each coefficient")
  }
  
  if( penalty == "SCAD" ){
    if( a <= 2 )
      stop("a must be > 2 for SCAD penalty")
  } else if ( penalty == "MCP" ){
    if( a <= 1 )
      stop("a must be > 1 for MCP penalty")
  } else if ( penalty == "MCL" ){
    if( a <= 1 )
      stop("a must be > 1 for MCL penalty")
  } else {
    if( penalty != "LASSO" )
      stop("wrong penalty function")
  }
  
  return(NULL)
}

model_eval <- function(model, test_x, test_y, test_w=NULL, func="check",...){
  #func: "check" (Quantile Check), "SqErr" (Squared Error), "AE" (Absolute Value)
  if(model$intercept){
    test_x <- cbind(1,test_x)
  }
  fits <- test_x %*% coefficients(model)
  eval_func <- switch(which(c("check","SqErr","AE")==func), check, square, abs)
  if(is.null(test_w)){
    mean(eval_func(test_y-fits,...)) 
  } else{
    weighted.mean(eval_func(test_y-fits,...), test_w)
  }
}


qbic <- function(model, largeP=FALSE){
  tau <- model$tau
  n <- model$n
  df <- sum(model$coefficients !=0)
  if(largeP){
    bic <- log(model$rho) + df*log(n)*log(length(model$coefficients))/(2*n)
  }else{
    bic <- log(model$rho) + df*log(n)/(2*n)
  }
  bic
}

coef.cv.rq.pen <- function(object, lambda='min',...){
  if(lambda=='min'){
    lambda <- object$lambda.min
  }
  target_model <- which(object$cv[,1] == lambda)
  coefficients(object$models[[target_model]])
}



cv.rq.pen <- function(x,y,tau=.5,lambda=NULL,weights=NULL,penalty="LASSO",intercept=TRUE,criteria="CV",cvFunc="check",nfolds=10,foldid=NULL,nlambda=100,eps=.0001,init.lambda=1,penVars=NULL,alg=ifelse(ncol(x)<50,"LP","QICD"),...){
  # x is a n x p matrix without the intercept term
  # y is a n x 1 vector
  # criteria used to select lambda is cross-validation (CV), BIC, or PBIC (large P)
  # nfolds: number of folds for cross validation
  # foldid: preset id of folds
  # penVar: variables to be penalized, default is all non-intercept terms
  
  # Pre-algorithm setup/ get convenient values
  m.c <- match.call() # This stores all the arguments in the function call as a list
  
  p <- dim(x)[2]
  if(is.null(penVars)){
    penVars <- 1:p
  }
  p_range <- penVars + intercept
  n <- dim(x)[1]
  pen_func <- switch(which(c("LASSO","SCAD","MCP")==penalty), lasso, scad, mcp)
  
  ### QICD ###
  if( alg=="QICD" & penalty!="LASSO" ){
    if(criteria=="CV"){
      stop("CV criteria not implemented for QICD algorithm with nonconvex penalties. Please use BIC or PBIC instead")
    }
    m.c[["alg"]] <- "LP" #maybe this should be moved inside the is.null initial beta if statement. I don't think it matters, but might be cleaner code
    penname <- penalty
    
    if( !all(penVars==1:p) ){ # Some unpenalized coefficients
      z    <- as.matrix(x[,-penVars])
      xpen <- as.matrix(x[,penVars])
      QICD_func <- "QICD.nonpen"
      mapback <- order( c(penVars, (1:p)[-penVars]) ) # reorders the coefficients properly if some (non-intercept) coefficients are not penalized 
      if( intercept )
        mapback <- c(1, 1+mapback)
    } else { # All penalized coefficients
      z <- NULL
      xpen <- x
      QICD_func <- "QICD"
      mapback <- 1:p # no reordering necessary if all (non-intercept) coefficients are penalized
      if( intercept )
        mapback <- c(1, 1+mapback)
    }
    
    # The QICD algorithm needs good starting values, use LASSO solution 
    ## Speed things up using BIC, not k-fold, to select lambda for LASSO
    ## Can skip this part if starting values are provided
    if( is.null(m.c[["initial_beta"]]) ){
      m.c[["penalty"]] <- "LASSO"
      m.c[["criteria"]] <- "BIC"
      if(is.null(m.c[["lambda"]])==FALSE){
        m.c[["lambda"]] <- NULL
      }
      suppressWarnings(
        m.c[["initial_beta"]] <- coefficients( eval.parent(m.c) )
        # QICD.start <- coefficients( cv.rq.pen(x,y,tau=tau,lambda=lambda,penalty="LASSO",intercept=intercept,criteria="BIC",nlambda=nlambda,eps=eps,init.lambda=lambda,penVars=penVars,...) ) # Use the LASSO with BIC
      )
    }
    
    # Start in middle of lambda vector
    ## Increase lambda until intercept only model (or all penlized coefficients are zero)
    ## Decrease lambda until full model (Sparsity assumption => full model is bad)
    m.c[[1]] <- as.name(QICD_func)
    m.c[["penalty"]] <- penalty
    m.c[["x"]] <- xpen
    m.c$z <- z
    
    if( is.null(lambda) ){
      lambdas <- exp( seq(-7, 1, length=100) ) 
    } else {
      lambdas <- lambda
    }
    
    coefs <- matrix(NA, p+intercept, length(lambdas))
    
    ### Loop for increasing lambda
    for( i in floor(length(lambdas)/2):length(lambdas) ){
      m.c[["lambda"]] <- lambdas[i]
      coefs[,i] <- eval.parent( m.c )[mapback]
      # coefs[,i] <- QICD_func( y=y,x=x, z=z, tau=tau, lambda=lambdas[i], intercept=intercept, 
      #                       penalty=penalty, initial_beta=QICD.start, ... )[mapback]
      if( all(coefs[p_range,i] == 0) )
        break()
    }
    
    ### Loop for decreasing lambda
    for( i in floor(length(lambdas)/2):2 -1 ){
      m.c[["lambda"]] <- lambdas[i]
      coefs[,i] <- eval.parent( m.c )[mapback]
      # coefs[,i] <- QICD_func( y=y,x=x, z=z, tau=tau, lambda=lambdas[i], intercept=intercept, 
      #                       penalty=penalty, initial_beta=QICD.start, ... )[mapback]
      if( all(coefs[p_range,i] != 0) )
        break()
    }
    
    #### Remove the NA columns from coefs and corresponding lambdas
    lambdas.keep <- which( !is.na(coefs[1,]) )
    lambdas <- lambdas[lambdas.keep]
    coefs <- coefs[,lambdas.keep]
    rownames(coefs) <- names( m.c[["initial_beta"]] )
    XB <- x%*%coefs[p_range,]
    if( intercept )
      XB <- XB + matrix(coefs[1,], n, ncol(coefs), byrow=TRUE)
    
    residuals <- y - XB
    rho <- colSums( check(residuals, tau=tau) )
    if( is.null(m.c[["a"]]) )
      a <- 3.7
    PenRho <- rho + colSums(apply( rbind(lambdas, coefs), 2, 
                                   function(xx) pen_func(xx[1+p_range], lambda=xx[1], a=a) ))
    
    cv <- data.frame(lambda=lambdas, cve=NA)
    
    if( criteria=="BIC" ){
      cv$cve <- log(rho) + colSums(coefs!=0)*log(n)/(2*n)
      names(cv)[2] <- "BIC"
    } else { # PBIC
      cv$cve <- log(rho) + colSums(coefs!=0)*log(n)*log(nrow(coefs))/(2*n)
      names(cv)[2] <- "PBIC"
    }
    
    lambda.min <- lambdas[which.min(cv[,2])]
    
    # Final cleanup for QICD
    ## First get models for each lambda
    models <- vector( "list", length(lambdas) )
    for( j in 1:length(models) ){
      models[[j]]$coefficients <- coefs[,j]
      models[[j]]$PenRho <- PenRho[j]
      models[[j]]$residuals <- residuals[,j]
      models[[j]]$rho <- rho[j]
      models[[j]]$tau <- tau
      models[[j]]$n <- n
      models[[j]]$intercept <- intercept
      models[[j]]$penalty <- penalty
      class(models[[j]]) <- c("rq.pen", "rqNC")
    }
    
    return_val <- list( models=models, cv=cv, lambda.min=lambda.min, penalty=penalty )
    class(return_val) <- "cv.rq.pen"
  } else{
    ############
    
    
    
    # If no lambdas provided, find reasonable lambdas to use
    if(is.null(lambda)){
      # find a lambda that sets all coefficients to zero. 
      # Strategy is to get \lambda \sum p_\lambda(|\beta_j}) >= \sum \rho_\tau(y-quantile(y,tau)
      # do this by fitting model with lambda = init.lambda and then set new lambda such that 
      # \lambda* = \sum \rho_\tau(y-quantile(y,tau)) / \sum p_\lambda(|beta_j|) repeat as needed
      sample_q <- quantile(y,tau)
      inter_only_rho <- sum(check(y-sample_q,tau))
      #lambda_star <- rep(0,p)
      #lambda_star[penVars] <- init.lambda
      lambda_star <- init.lambda
      searching <- TRUE
      while(searching){
        if(penalty=="LASSO"){
          init_fit <- rq.lasso.fit(x,y,tau,lambda=lambda_star,weights,intercept,penVars=penVars,...)
        } else{
          init_fit <- rq.nc.fit(x,y,tau,lambda=lambda_star,weights,intercept,penVars=penVars,...)
        }
        if(sum(init_fit$coefficients[p_range])==0){
          searching <- FALSE     
        } else{
          lambda_star <- inter_only_rho / sum(sapply(init_fit$coefficients[p_range],pen_func,1)) 
          #1 used here because solving for lambda
        }
      }
      lambda_min <- eps*lambda_star
      lambda <- exp(seq(log(max(lambda_min)),log(max(lambda_star)),length.out=nlambda))#max is included to handle cases where
      # some variables are unpenalized and thus lambda is a multivalued vector with some zeros
    }
    # lambda is the vector of reasonable choices of lambda to use in the penalty
    
    models <- list()
    fit_models <- TRUE
    lam_pos <- 1
    if(penalty=="LASSO"){
      while(fit_models){
        if(fit_models){
          models[[lam_pos]] <- rq.lasso.fit(x,y,tau,lambda[lam_pos],weights,intercept,penVars=penVars,...)
          
        }
        if(sum(abs(coefficients(models[[lam_pos]])[p_range]))==0 || lam_pos==length(lambda)){
          #if we got a fully sparse model, no need to fit more sparse models
          fit_models <- FALSE
          lambda <- lambda[1:lam_pos]
        }
        lam_pos <- lam_pos + 1
      }
    } else{
      while(fit_models){
        if(fit_models){
          models[[lam_pos]] <- rq.nc.fit(x,y,tau,lambda[lam_pos],weights,intercept,penalty=penalty,penVars=penVars,...)
        }
        if(sum(abs(coefficients(models[[lam_pos]])[p_range]))==0 || lam_pos==length(lambda)){
          fit_models <- FALSE
          lambda <- lambda[1:lam_pos]
        }
        lam_pos <- lam_pos + 1
      }
      #models <- lapply(lambda,rq.nc.fit, x=x,y=y,tau=tau,weights=weights,intercept=intercept,penalty=penalty,
      #                                   penVars=penVars,...)
    }
    cv_results <- NULL
    if(criteria=="CV"){
      if(is.null(foldid)){
        foldid <- randomly_assign(n,nfolds)
      }
      for(i in 1:nfolds){
        train_x <- x[foldid!=i,]
        train_y <- y[foldid!=i]
        test_x <- x[foldid==i,,drop=FALSE]
        test_y <- y[foldid==i]
        train_weights <- weights[foldid!=i] #not sure this line is needed
        if(is.null(weights)){
          train_weights <- test_weights <- NULL
        } else{
          train_weights <- weights[foldid!=i]
          test_weights <- weights[foldid==i]
        }
        if(penalty=="LASSO"){
          cv_models <- lapply(lambda,rq.lasso.fit, x=train_x,y=train_y,tau=tau,weights=train_weights,intercept=intercept,penVars=penVars,...)
        } else{
          cv_models <- lapply(lambda,rq.nc.fit, x=train_x,y=train_y,tau=tau,weights=train_weights,intercept=intercept,penalty=penalty,penVars=penVars,...)
        }
        if(cvFunc=="check"){
          cv_results <- cbind(cv_results, sapply(cv_models,model_eval, test_x, test_y, test_weights, tau=tau))
        } else{
          cv_results <- cbind(cv_results, sapply(cv_models,model_eval, test_x, test_y, test_weights, func=cvFunc))
        } 
      }
      cv_results <- apply(cv_results,1,mean)
    }
    if(criteria=="BIC"){
      cv_results <- sapply(models,qbic)
    }
    if(criteria=="PBIC"){
      cv_results <- sapply(models,qbic,largeP=TRUE)
    }
    lambda.min <- lambda[which.min(cv_results)]
    return_val <- NULL
    return_val$models <- models
    return_val$cv <- data.frame(lambda=lambda, cve=cv_results)
    colnames(return_val$cv)[2] <- criteria
    return_val$lambda.min <- lambda.min
    return_val$penalty <- penalty
    class(return_val) <- "cv.rq.pen"
  }
  
  return_val
}


re_order_nonpen_coefs <- function(nonpen_coefs, penVars, intercept=TRUE){
  p <- length(nonpen_coefs)
  new_coefs <- rep(NA,p)
  if(intercept){
    penVars <- penVars+1
    pen_output <- 2:(length(penVars)+1)
  } else{
    pen_output <- 1:length(penVars)
  }
  new_coefs[penVars] <- nonpen_coefs[pen_output]
  new_coefs[-penVars] <- nonpen_coefs[-pen_output]
  new_coefs
}


rq.nc.fit <- function(x,y,tau=.5,lambda=NULL,weights=NULL,intercept=TRUE,
                      penalty="SCAD",a=3.7,iterations=10,converge_criteria=1e-06,
                      alg=ifelse(p<50,"LP","QICD"),penVars=NULL,...){
  # x is a n x p matrix without the intercept term
  # y is a n x 1 vector
  # lambda takes values of 1 or p
  # penalty SCAD or MCP
  # penVars - variables to be penalized, doesn't work if lambda has multiple entries
  p <- ncol(x)
  n <- nrow(x)
  
  if( alg=="QICD" ){
    ### QICD Algorithm ###
    coefnames <- paste("x",1:p, sep="") ### Coefficient names
    if( length(lambda) != 1 )
      stop( "QICD Algorithm only allows 1 lambda value")	
    ### Check if we are using QICD or QICD.nonpen
    if( is.null(penVars) | length(penVars) == p){ ### No unpenalized coefficients
      
      coefs <- QICD(y, x, tau, lambda, intercept, penalty, eps=converge_criteria, a=a, ...)
      penbeta <- intercept + 1:p ### Use later to calculate objective function
    } else { ### Some unpenalized coefficients
      z    <- as.matrix(x[,-penVars])
      xpen <- as.matrix(x[,penVars])
      #coefnames <- paste("x",1:ncol(xpen), sep="") ### Coefficient names
      #coefnames <- c( coefnames, paste("z",1:ncol(z), sep="") )
      penbeta <- intercept + penVars ### Use later to calculate objective function
      coefs <- QICD.nonpen(y, xpen, z, tau, lambda, intercept, penalty, eps=converge_criteria, a=a, ...)
      coefs <- re_order_nonpen_coefs(coefs, penVars, intercept)
    }
    
    ### Add extra information to QICD output
    
    if( intercept ){ ### Residuals
      residuals <- c( y - x%*%(coefs[-1]) - coefs[1] )
      names(coefs) <- c("intercept",coefnames)
    } else {
      residuals <- c( y - x%*%coefs )
      names(coefs) <- coefnames
    }
    rho <- sum( check(residuals) )
    #1/n*sum( check(residuals) ) ### rho
    if( penalty == "LASSO" ){ ### PenRho for LASSO
      PenRho <- sum( abs( coefs[penbeta] )*lambda )
    } else if( penalty == "SCAD" ){ ### PenRho for SCAD
      PenRho <- sum( scad( coefs[penbeta], lambda, a ))
    } else { ### PenRho for MCP
      PenRho <- sum( mcp( coefs[penbeta], lambda, a ))
    }
    PenRho <- rho + PenRho
    
    sub_fit <- list( coefficients=coefs,  PenRho=PenRho, residuals=residuals,
                     rho=rho, tau=tau, n=n , intercept=intercept)
    ######################
  } else {  
    ###  LP Algorithm  ###
    if(penalty=="SCAD"){
      deriv_func <- scad_deriv
    }
    if(penalty=="MCP"){
      deriv_func <- mcp_deriv
    }
    if(is.null(dim(x))){                                                                                    
      stop('x needs to be a matrix with more than 1 column')
    }
    if(n != length(y)){
      stop('length of y and rows of x do not match')
    }
    if(is.null(lambda)==TRUE | (length(lambda) != 1 & length(lambda) != dim(x)[2])){
      stop(paste('input of lambda must be of length 1 or', dim(x)[2]))
    }
    if( sum(lambda < 0) > 0){
      stop('lambda must be positive')
    }
    if(is.null(penVars) !=TRUE & length(lambda) == 1){
      mult_lambda <- rep(0,p)
      mult_lambda[penVars] <- lambda
      lambda <- mult_lambda
    }
    if(length(lambda) != 1){
      penVars <- (1:p)[lambda != 0]
      pen_range <- intercept + penVars
    } else{
      pen_range <- intercept + 1:p
      if(is.null(penVars)){
        penVars <- 1:p
      }
    }
    #lambda_update <- n*lambda
    lambda_update <- lambda
    iter_complete <- FALSE
    iter_num <- 0
    old_beta <- rep(0, p+intercept)
    while(!iter_complete){
      sub_fit <- rq.lasso.fit(x=x,y=y,tau=tau,lambda=lambda_update,weights=weights,intercept=intercept,...)
      if(length(lambda)==1){
        lambda_update <- sapply(abs(sub_fit$coefficients[pen_range]),deriv_func, lambda=lambda, a=a)
      } else{
        lambda_update[-penVars] <- 0
        lambda_update[penVars] <- mapply(deriv_func, x=abs(sub_fit$coefficients[pen_range]),
                                         lambda=lambda[penVars],
                                         MoreArgs=list(a=a))
      }
      #lambda_update <- n*lambda_update
      iter_num <- iter_num + 1
      new_beta <- sub_fit$coefficients
      beta_diff <- sum( (old_beta - new_beta)^2)
      if(iter_num == iterations | beta_diff < converge_criteria){
        iter_complete <- TRUE
        if(iter_num == iterations & beta_diff > converge_criteria){
          warning(paste("did not converge after ", iterations, " iterations", sep=""))
        }
      } else{
        old_beta <- new_beta
      }
    }
    ######################
  }
  sub_fit$penalty <- penalty
  class(sub_fit) <-  c("rq.pen", "rqNC")
  return_val <- sub_fit
  return(return_val)
}

beta_plots <- function(model,voi=NULL,logLambda=TRUE,loi=NULL,...){
  #voi - index variables of interest
  #logLambda - lambdas plotted on log scale
  #loi - index of target lambdas
  if( "cv.rq.group.pen" %in% class(model)){
    betas <- t( model$beta)
    if(model$intercept){
      betas <- betas[,-1]
    }
  }
  else{
    betas <- t(sapply(model$models, coefficients))
    if(is.null(voi)==FALSE){
      betas <- betas[,voi]
    }
    if(colnames(betas)[1]=="intercept"){
      betas <- betas[,-1]
    }
  }
  if(logLambda){
    lambdas <- log(model$cv$lambda)
  } else{
    lambdas <- model$cv$lambda
  }
  
  if(is.null(loi)==FALSE){
    lambdas <- lambdas[loi]
  }                                    
  plot(lambdas, betas[,1], type="n",ylim=c(min(betas),max(betas)),ylab="Coefficient Value",xlab="Log Lambda",...)
  for(i in 1:dim(betas)[2]){
    lines(lambdas, betas[,i],col=i)
  }  
}

cv_plots <- function(model,logLambda=TRUE,loi=NULL,...){
  #logLambda - lambdas plotted on log scale
  #loi - index of target lambdas
  cv_data <- model$cv
  if(logLambda){
    cv_data$lambda <- log(cv_data$lambda)
    colnames(cv_data)[1] <- "logLambda"
  }
  if(is.null(loi)==FALSE){
    cv_data <- cv_data[loi,]                                        
  }                                    
  plot(cv_data[,1], cv_data[,2], ylab=colnames(cv_data)[2], xlab=colnames(cv_data)[1],...)
}

getRho <- function(model){
  model$rho
}

cv.rq.group.pen <- function (x, y, groups, tau = 0.5, lambda = NULL, penalty = "SCAD", 
                             intercept = TRUE, criteria = "CV", cvFunc = "check", nfolds = 10, 
                             foldid = NULL, nlambda = 100, eps = 1e-04, init.lambda = 1,alg="QICD",penGroups=NULL,
                             ...) 
{
  if(penalty=="LASSO"){
    warning("Group penalties use the L1 norm and the Lasso group penalty is the same as the standard Lasso penalty and therefore does not account for group structure. The group lasso method is only implemented because it is needed for the SCAD and MCP algorithms. Otherwise it should be avoided. ")
  }
  if(is.null(penGroups)){
    p_range <- 1:dim(x)[2] + intercept
  } else{
    p_range <- which(groups %in% penGroups) + intercept
  }
  n <- dim(x)[1]
  pen_func <- switch(which(c("LASSO", "SCAD", "MCP") == penalty), 
                     lasso, scad, mcp)
  if(alg=="QICD" & penalty != "LASSO"){
    if(criteria== "CV"){
      stop("QICD algorithm wtih non-convex penalties setup only to use BIC or PBIC as the criteria")
    }
    #start with lasso fit
    lasso_fit <- cv.rq.group.pen(x,y,groups,tau,lambda,penalty="LASSO",intercept,criteria="BIC",cvFunc,nfolds,foldid,
                                 nlambda, eps, init.lambda,alg="LP",penGroups,...)
    #then iterate through lasso models to get new models
    model_coefs <- NULL
    lambda_vals <- lasso_fit$cv$lambda
    lasso_beta <- lasso_fit$beta
    for(model_num in 1:dim(lasso_beta)[2]){
      model_coefs <- cbind(model_coefs, QICD.group(y, x, groups, tau, lambda_vals[model_num], intercept, penalty, 
                                                   initial_beta=lasso_beta[,model_num], eps = eps,...))
    }
    return_val <- NULL
    return_val$beta <- model_coefs
    if(intercept){
      fits <- cbind(1,x) %*% model_coefs
    } else{
      fits <- x %*% model_coefs
    }
    return_val$residuals <- y - fits
    return_val$rho <- apply(check(return_val$residuals,tau),2,sum)
    non_zero_count <- apply(model_coefs!=0,2,sum)
    if( criteria=="BIC" ){
      cve <- log(return_val$rho) + non_zero_count*log(n)/(2*n)
    } else { # PBIC
      cve <- log(return_val$rho) + non_zero_count*log(n)*log(nrow(model_coefs))/(2*n)
    }
    
    return_val$cv <- data.frame(lambda = lambda_vals, cve = cve)
    colnames(return_val$cv)[2] <- criteria
    return_val$lambda.min <- lambda_vals[which.min(cve)]
    return_val$penalty <- penalty
    return_val$intercept <- intercept
    return_val$groups <- groups
    class(return_val) <- c("cv.rq.group.pen", "cv.rq.pen")    
  }      
  else{
    if (is.null(lambda)) {
      sample_q <- quantile(y, tau)
      inter_only_rho <- sum(check(y - sample_q, tau))
      lambda_star <- init.lambda
      searching <- TRUE
      search_num <- 1
      while(searching){
        if (search_num == 1) {
          init_fit <- rq.group.fit(x, y, groups, tau, lambda_star, 
                                   intercept, penalty,alg,penGroups=penGroups, ...)
          search_num <- 2
        }
        else {
          init_fit <- rq.group.fit(x, y, groups, tau, lambda_star, 
                                   intercept, penalty,alg, initial_beta = beta_update, penGroups=penGroups,
                                   ...)
        }
        beta_update <- init_fit$coefficients
        if (sum(init_fit$coefficients[p_range]) == 0) {
          searching <- FALSE
        }
        else {
          option_1 <- (inter_only_rho-init_fit$rho)/ sum(sapply(init_fit$coefficients[p_range],pen_func, 1))
          option_2 <- lambda_star*search_num
          if(option_1 >= option_2){
            lambda_star <- option_1
          } else{
            lambda_star <- option_2
            search_num <- search_num + 1
          }
          # lambda_star = max((inter_only_rho-init_fit$rho)/ sum(sapply(init_fit$coefficients[p_range],pen_func, 1)),
          #                    lambda_star+1)
          #this is sort of a hack, need to think of a better way to pick lambda_star
          #lambda_star = inter_only_rho/sum(sapply(init_fit$coefficients[p_range], 
          #  pen_func, 1))
          #weird behavior if we don't set penalty lambda to 1
          #in general this idea needs to be improved upon
        }
      }
      
      lambda_min <- eps * lambda_star
      lambda <- exp(seq(log(lambda_star), log(lambda_min),  length.out = nlambda))
      fine_tune <- TRUE
      fine_tune_pos <- length(lambda)-1
      while(fine_tune){
        test_fit <- rq.group.fit(x,y,groups,tau,lambda[fine_tune_pos],intercept,penalty,alg,penGroups=penGroups,...)
        if (sum(test_fit$coefficients[p_range]) != 0) {
          #want only one intercept only model
          fine_tune <- FALSE
        } else{
          fine_tune_pos <- fine_tune_pos - 1
        }
      }
      if(fine_tune_pos != (length(lambda)-1)){
        lambda_star <- lambda[fine_tune_pos+1]
        lambda_min <- eps*lambda_star
        lambda <- exp(seq(log(lambda_star), log(lambda_min),  length.out = nlambda))
      }
    }
    
    models <- groupMultLambda(x = x, y = y, groups = groups, 
                              tau = tau, lambda = lambda, intercept = intercept, penalty=penalty,alg=alg,penGroups=penGroups, ...)
    model_coefs  <- sapply(models,coefficients)
    model_resids <- sapply(models,residuals)
    model_rhos <- sapply(models, getRho)
    
    cv_results <- NULL
    if (criteria == "CV") {
      if (is.null(foldid)) {
        foldid <- randomly_assign(n, nfolds)
      }
      for (i in 1:nfolds) {
        train_x <- x[foldid != i, ]
        train_y <- y[foldid != i]
        #test_x <- x[foldid == i, ]
        test_x <- x[foldid==i,,drop=FALSE]
        test_y <- y[foldid == i]
        cv_models <- groupMultLambda(x = train_x, y = train_y, 
                                     groups = groups, tau = tau, lambda = lambda, 
                                     intercept = intercept,penalty=penalty,alg=alg,penGroups=penGroups, ...)
        if (cvFunc == "check") {
          cv_results <- cbind(cv_results, sapply(cv_models, 
                                                 model_eval, test_x, test_y, tau = tau))
        }
        else {
          cv_results <- cbind(cv_results, sapply(cv_models, 
                                                 model_eval, test_x, test_y, func = cvFunc))
        }
      }
      cv_results <- apply(cv_results, 1, mean)
    }
    if (criteria == "BIC") {
      cv_results <- sapply(models, qbic)
    }
    if (criteria == "PBIC") {
      cv_results <- sapply(models, qbic, largeP = TRUE)
    }
    lambda.min <- lambda[which.min(cv_results)]
    return_val <- NULL
    #return_val$models <- models
    return_val$beta <- model_coefs
    return_val$residuals <- model_resids
    return_val$rho <- model_rhos   
    return_val$cv <- data.frame(lambda = lambda, cve = cv_results)
    colnames(return_val$cv)[2] <- criteria
    return_val$lambda.min <- lambda.min
    return_val$penalty <- penalty
    return_val$intercept <- intercept
    return_val$groups <- groups
    class(return_val) <- c("cv.rq.group.pen", "cv.rq.pen")
  }
  return_val
}

rq.group.fit <- function (x, y, groups, tau = 0.5, lambda, intercept = TRUE, 
                          penalty = "SCAD", alg="QICD", a=3.7,penGroups=NULL, ...) 
{
  ### Some cleaning/checking before getting to the algorithms
  p <- ncol(x)
  n <- nrow(x)
  #if(is.null(penGroups) & max(penGroups) > max(groups)){ stop("penalize groups not coefficients")}  
  if (!penalty %in% c("LASSO", "SCAD", "MCP")) {
    stop("Penalty must be LASSO, SCAD or MCP")
  }
  if(penalty=="LASSO"){
    warning("Group penalties use the L1 norm and the Lasso group penalty is the same as the standard Lasso penalty and therefore does not account for group structure. The group lasso method is only implemented because it is needed for the SCAD and MCP algorithms. Otherwise it should be avoided. ")
  }
  if(is.null(dim(x))){ stop("x must be matrix with at least 1 column") }
  if(length(groups)!=ncol(x)){
    stop("length(groups) must be equal to ncol(x)")
  }
  if( lambda <= 0 ){ stop("lambda must be positive")}
  
  if(penalty=="LASSO"){
    pen_func <- lasso
  }
  if(penalty=="SCAD"){
    pen_func <- scad
  }
  if(penalty=="MCP"){
    pen_func <- mcp
  }
  
  if (alg == "QICD") {
    ### QICD Algorithm ###
    if( length(lambda) != 1 )
      stop( "QICD Algorithm only allows 1 lambda value")
    
    coefs <- QICD.group(y, x, groups, tau, lambda, intercept, penalty,a=a, ...)
    
    ### Add extra information to QICD output
    coefnames <- paste("x",1:p, sep="") ### Coefficient names
    if(intercept)
      coefnames <- c("(Intercept)", coefnames)
    names(coefs) <- coefnames
    if( intercept ){ ### Residuals
      residuals <- c( y - x%*%(coefs[-1]) - coefs[1] )
      pen_vars <- coefs[-1]
    } else {
      residuals <- c( y - x%*%coefs )
      pen_vars <- coefs
    }
    if(penalty=="LASSO"){
      pen_val <- sum(pen_func(tapply(abs(pen_vars),groups,sum),lambda=lambda))
    } else{
      pen_val <- sum(pen_func(tapply(abs(pen_vars),groups,sum),lambda=lambda,a=a))
    }
    
    rho <- sum( check(residuals) ) # rho (similiar to quantreg)
    #1/n*sum( check(residuals) ) ### rho
    PenRho <- rho+pen_val
    
    return_val <- list( coefficients=coefs, PenRho=PenRho, residuals=residuals, rho=rho, tau=tau, n=n, intercept=intercept, penalty=penalty)
    class(return_val) <- c("rq.group.pen", "rq.pen")
    ######################
  } else {
    group_num <- length(unique(groups))
    if (length(lambda) == 1) {
      lambda <- rep(lambda, group_num)
    }
    if (length(lambda) != group_num) {
      stop("lambdas do not match with group number")
    }
    if (sum(groups == 0) > 0) {
      stop("0 cannot be used as a group")
    }
    if (dim(x)[2] != length(groups)) {
      stop("length of groups must be equal to number of columns in x")
    }
    
    if (penalty == "LASSO") {			
      new_lambda <- NULL
      group_count <- xtabs(~groups)
      for (g in 1:group_num) {
        new_lambda <- c(new_lambda, rep(lambda[g], each = group_count[g]))
      }
      if(is.null(penGroups)==FALSE){
        new_lambda[-which(groups %in% penGroups)] <- 0
      }
      return_val <- rq.lasso.fit(x, y, tau, new_lambda, 
                                 intercept = intercept, ...)
      class(return_val) <- c("rq.group.pen", "rq.pen", 
                             "rqLASSO")
    }
    else {
      return_val <- rq.group.lin.prog(x,y,groups,tau,lambda,intercept=intercept,penalty=penalty,penGroups=penGroups,a=a,...)
      class(return_val) <- c("rq.group.pen", "rq.pen")
    }
  }
  return_val
}

plot.cv.rq.group.pen <- function (x,y=NULL,...) 
{
  plot(x$cv[, 1], x$cv[, 2])
}

qaSIS <- function(x,y,tau=.5,linear=FALSE,...){#n.cores=1,...){
  if(linear){
    eval_function<- function(x,y,tau){ 
      q1 <- rq(y ~ x, tau)
      sum((fitted(q1)-quantile(y,tau))^2)
    }
    eval_results <- apply(x,2,eval_function,y,tau,...)
  } else{
    eval_function2 <- function(x,y,tau,...){ 
      b <- bs(x,...)
      q1 <- rq(y ~ b, tau)
      sum((fitted(q1)-quantile(y,tau))^2)
    }
    eval_results <- apply(x,2,eval_function2,y,tau,...)
  }
  #if(n.cores==1){
  
  #} else{
  #	p <- dim(x)[2]
  #	mc_func <- function(idx,...){ eval_function(x[,idx],y,...)}
  #	mc_results <- mclapply(1:p, mc_func, mc.cores=n.cores, ...)
  #	eval_results <- do.call(c,mc_results)
  #}
  order( eval_results, decreasing=TRUE)
}


####################################################################3
library(Rcpp)

setwd("C:/Users/sypar/sy_work/R_work")
getwd()

Rcpp::sourceCpp("C:/Users/sypar/sy_work/R_work/src/QCD.cpp")

# 예제 적용
library(mlbench)

data("BostonHousing")
head(BostonHousing, 5)

X = data.matrix(BostonHousing[, -14])
y = BostonHousing[, 14]

head(x)

QICD(y, X, tau = 0.75, 2, 1, intercept = TRUE, penalty = "MCL", initial_beta = NULL, a = 3.7, scalex = TRUE, )

