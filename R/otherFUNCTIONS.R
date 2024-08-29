## The function performs penalized G-estimation for a given data, working
## correlation structure and a single value of the tuning parameter

#' Function to perform penalized G-estimation for a single value of tuning parameter
#' @description This function performs penalized G-estimation for a given longitudinal data,
#' a specific working correlation structure and a single value of the tuning parameter.
#' @param data A data frame containing the variables in longitudinal format. The name of the
#' variables ID, exposure and outcome must be "id", "a" and "y", respectively, in the data.
#' @param wc.str A character string specifying the working correlation structure. The 
#' following are currently allowed: "independence", "exchangeable", "ar1", and "unstructured".
#' @param tf.model A single formula object specifying the covariates of a (linear)
#' treatment-free model.
#' @param treat.model A single formula object specifying the covariates of the logistic
#' regression model for the binary treatment.
#' @param lambda A single value of the tuning parameter.
#' @param maxitr Maximum number of iterations allowed.
#' @param penalty The penalty type to be used, available options include "SCAD" and "MCP". 
#' "SCAD" refers to the Smoothly Clipped Absolute Deviation penalty and "MCP" refers to the Minimax Concave Penalty.

#'
#' @return A list containing the estimates, asymptotic variance etc.
#' @export

penG <- function(data, wc.str, tf.model, treat.model, lambda, maxitr, penalty){
  # n: number of patients
  # p: number of variables in the treatment-free model including the intercept
  # ni, l, a, y, l.mat.split, e: all are list of length n
  # estimate.current: a vector of length 2*p

  n <- length(split(data, data$id))
  dat <- split(data, data$id)
  l.mat <- data.frame(id=data$id,model.matrix(tf.model, data)) # cov+treat history
  l.mat.split <- split(l.mat, l.mat$id)

  ni <- lapply(1:n, function(i) length(dat[[i]]$id))
  l <- lapply(1:n, function(i) as.matrix(l.mat.split[[i]][,-1]))
  a <- lapply(1:n, function(i) as.matrix(dat[[i]]$a))
  y <- lapply(1:n, function(i) as.matrix(dat[[i]]$y))
  E.a <- lapply(1:n, function(i) as.matrix(dat[[i]]$E.a))
  p <- dim(l.mat[,-1])[2] # dimension of psi/delta including intercept

  sum1 <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%l[[i]])) # p*p
  sum2 <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%(c(a[[i]])*l[[i]]))) # p*p
  sum3 <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%y[[i]])) # p*1
  sum4 <- Reduce("+", lapply(1:n, function(i) t(c(a[[i]]-E.a[[i]])*l[[i]])%*%(c(a[[i]])*l[[i]])))# p*p
  sum5 <- Reduce("+", lapply(1:n, function(i) t(c(a[[i]]-E.a[[i]])*l[[i]])%*%y[[i]])) #p*1
  sum6 <- Reduce("+", lapply(1:n, function(i) t(c(a[[i]]-E.a[[i]])*l[[i]])%*%l[[i]])) #p*p


  ## initial estimate, iteration = 0
  estimate.current <- solve(cbind(rbind(sum4, sum2), rbind(sum6, sum1)))%*%rbind(sum5, sum3)
  e <- lapply(1:n, function(i)
    y[[i]]-(c(a[[i]])*l[[i]])%*%estimate.current[1:p]-l[[i]]%*%estimate.current[(p+1):(2*p)])
  sigma2.hat <- phi.hat.fun(n, ni, e)
  alpha.hat <- switch(wc.str,
                      "exchangeable" = corr.exch.fun(n, ni, e, phi = sigma2.hat),
                      "ar1" = corr.ar1.fun(n, ni, e, phi = sigma2.hat),
                      "independence" = 0,
                      "unstructured" = corr.un.fun(n, ni, e, phi = sigma2.hat))
  itr <- 0

  ### iteration = 1

  ## corr.mat, V, V.inv >> list of length n
  V <- lapply(1:n, function(i) sigma2.hat*corMatF(alpha = alpha.hat, ni = ni[[i]], corstr = wc.str))

  rcond <- lapply(V, function(x) rcond(x))
  if (sum(unlist(rcond) <= 1e-15) == 0){
    V.inv <- lapply(1:n, function(i) solve(V[[i]]))
    sum1.n <- Reduce("+", lapply(1:n, function(i) t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%e[[i]]))
    sum2.n <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%V.inv[[i]]%*%e[[i]]))
    sum3.n <- Reduce("+", lapply(1:n, function(i)
      t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%(c(a[[i]])*l[[i]])))
    sum4.n <- Reduce("+", lapply(1:n, function(i) t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%l[[i]]))
    sum6.n <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%V.inv[[i]]%*%l[[i]]))
    sum5.n <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%V.inv[[i]]%*%(c(a[[i]])*l[[i]])))
    penalty.der <- switch(penalty,
                          "SCAD"=q.scad(lambda, abs(estimate.current[2:p])),
                          "MCP"=q.mcp(lambda, abs(estimate.current[2:p])))
    En.mat <- diag(c(0, penalty.der/(10^(-6)+abs(estimate.current[2:p])),
                     rep(0,p)))
    estimate.new <- estimate.current + solve(rbind(cbind(sum3.n, sum4.n), cbind(sum5.n, sum6.n))+n*En.mat)%*%(rbind(sum1.n, sum2.n)-n*En.mat%*%estimate.current)
    e <- lapply(1:n, function(i)
      y[[i]]-(c(a[[i]])*l[[i]])%*%estimate.new[1:p]-l[[i]]%*%estimate.new[(p+1):(2*p)])
    sigma2.hat <- phi.hat.fun(n, ni, e)
    alpha.hat <- switch(wc.str,
                        "exchangeable" = corr.exch.fun(n, ni, e, phi = sigma2.hat),
                        "ar1" = corr.ar1.fun(n, ni, e, phi = sigma2.hat),
                        "independence" = 0,
                        "unstructured" = corr.un.fun(n, ni, e, phi = sigma2.hat))
    itr <- 0+1

    ### iterations: 2, 3, ... , till convergence/maxitr
    while(max(abs(estimate.new-estimate.current)) > 0.000001 & itr<= maxitr &
          sum(unlist(rcond) <= 1e-15) == 0){

      estimate.current <- estimate.new
      penalty.der <- switch(penalty,
                            "SCAD"=q.scad(lambda, abs(estimate.current[2:p])),
                            "MCP"=q.mcp(lambda, abs(estimate.current[2:p])))
      En.mat <- diag(c(0, penalty.der/(10^(-6)+abs(estimate.current[2:p])),
                       rep(0,p)))
      ## corr.mat, V, V.inv >> list of length n
      V <- lapply(1:n, function(i) sigma2.hat*corMatF(alpha = alpha.hat, ni = ni[[i]], corstr = wc.str))

      rcond <- lapply(V, function(x) rcond(x))
      if (sum(unlist(rcond) <= 2.220446e-16) == 0){
        V.inv <- lapply(1:n, function(i) solve(V[[i]]))

        sum1.n <- Reduce("+", lapply(1:n, function(i) t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%e[[i]]))
        sum2.n <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%V.inv[[i]]%*%e[[i]]))
        sum3.n <- Reduce("+", lapply(1:n, function(i)
          t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%(c(a[[i]])*l[[i]])))
        sum4.n <- Reduce("+", lapply(1:n, function(i) t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%l[[i]]))
        sum6.n <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%V.inv[[i]]%*%l[[i]]))
        sum5.n <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%V.inv[[i]]%*%(c(a[[i]])*l[[i]])))

        estimate.new <- estimate.current + solve(rbind(cbind(sum3.n, sum4.n), cbind(sum5.n, sum6.n))+n*En.mat)%*%(rbind(sum1.n, sum2.n)-n*En.mat%*%estimate.current)
        e <- lapply(1:n, function(i)
          y[[i]]-(c(a[[i]])*l[[i]])%*%estimate.new[1:p]-l[[i]]%*%estimate.new[(p+1):(2*p)])
        sigma2.hat <- phi.hat.fun(n, ni, e)
        alpha.hat <- switch(wc.str,
                            "exchangeable" = corr.exch.fun(n, ni, e, phi = sigma2.hat),
                            "ar1" = corr.ar1.fun(n, ni, e, phi = sigma2.hat),
                            "independence" = 0,
                            "unstructured" = corr.un.fun(n, ni, e, phi = sigma2.hat))
        itr <- itr+1
      }
    }
  }

  ##Asymptotic variance calculation
  if (itr >0 & itr <= maxitr & sum(unlist(rcond) <= 1e-15) == 0){
    l.mat.treat <- data.frame(id=data$id, model.matrix(treat.model,data)) # cov+treat history
    l.mat.treat.split <- split(l.mat.treat, l.mat.treat$id)
    l.treat <- lapply(1:n, function(i) as.matrix(l.mat.treat.split[[i]][,-1]))
    V <- lapply(1:n, function(i) sigma2.hat*corMatF(alpha = alpha.hat, ni = ni[[i]], corstr = wc.str))
    V.inv <- lapply(1:n, function(i) solve(V[[i]]))
    sum3.n <- Reduce("+", lapply(1:n, function(i)
      t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%(c(a[[i]])*l[[i]])))
    ss1 <- Reduce("+", lapply(1:n, function(i)
      t(l.treat[[i]])%*%(a[[i]]-E.a[[i]])%*%t(a[[i]]-E.a[[i]])%*%l.treat[[i]]))
    ss2 <- Reduce("+", lapply(1:n, function(i)
      t(l.treat[[i]])%*%(a[[i]]-E.a[[i]])%*%t(e[[i]])%*%V.inv[[i]]%*%(c(a[[i]]-E.a[[i]])*l[[i]])))
    ss4 <- Reduce("+", lapply(1:n, function(i)
      t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%e[[i]]%*%t(e[[i]])%*%
        V.inv[[i]]%*%(c(a[[i]]-E.a[[i]])*l[[i]])))

    I.psi.hat <- ss4 - t(ss2)%*%solve(ss1)%*%ss2
    B.hat.inv <- solve(sum3.n+n*En.mat[1:p,1:p])
    var.psi.hat <- B.hat.inv%*%I.psi.hat%*%B.hat.inv
    return(list(estimate = estimate.new, phi = sigma2.hat, alpha = alpha.hat,
                asymp.var.psi = var.psi.hat, nitr=itr, error = 0, En.mat = En.mat))
  } else {
    #print("Estimation did not converge")
    return(list(error = 1, nitr=itr))
  }
}

### Other functions to be called

## Method of moments estimate for the variance parameter sigma^2
#' Function to perform method of moments estimation for the variance parameter
#'
#' @param n A scalar value indicating the number of subjects.
#' @param ni A list of length n containing the number of measurements per subject.
#' @param e A list of length n containing the residuals.
#'
#' @return Method of moments estimate for the variance parameter
#' @export

phi.hat.fun <- function(n, ni, e){
  sum.phi.i <- lapply(1:n, function(i) (1/ni[[i]])*sum((e[[i]])^2))
  sigma2.hat <- sum(unlist(sum.phi.i))/n
  return(sigma2.hat)
}

## Method of moments estimate for the correlation parameter alpha under exchangeable structure
#' Function to perform method of moments estimation for the correlation parameter alpha under exchangeable structure
#'
#' @param n A scalar value indicating the number of subjects.
#' @param ni A list of length n containing the number of measurements per subject.
#' @param e A list of length n containing the residuals.
#' @param phi The estimate of variance parameter.
#'
#' @return Method of moments estimate for the correlation parameter alpha under exchangeable structure
#' @export

corr.exch.fun <- function(n, ni, e, phi){
  sum.corr.i <- lapply(1:n, function(i) (1/(ni[[i]]*(ni[[i]]-1)))*(sum(e[[i]]%*%t(e[[i]]))
                                                                   -sum(diag(e[[i]]%*%t(e[[i]])))))
  corr <- sum(unlist(sum.corr.i))/(phi*n)
  return(corr) ####
}

## Method of moments estimate for the correlation parameter alpha under AR(1) structure
#' Function to perform method of moments estimation for the correlation parameter alpha under AR(1) structure
#'
#' @param n A scalar value indicating the number of subjects.
#' @param ni A list of length n containing the number of measurements per subject.
#' @param e A list of length n containing the residuals.
#' @param phi The estimate of variance parameter.
#'
#' @return Method of moments estimate for the correlation parameter alpha under AR(1) structure
#' @export
corr.ar1.fun <- function(n, ni, e, phi){
  sum.corr.i <- lapply(1:n, function(i) (1/(ni[[i]]-1))*(t(e[[i]][1:(ni[[i]]-1)])%*%(e[[i]][2:ni[[i]]])))
  corr <- sum(unlist(sum.corr.i))/(phi*n)
  return(corr) ####
}

## Method of moments estimate for the correlation parameters under unstructured correlation structure
#' Function to perform method of moments estimation for the correlation parameter alpha under unstructured correlation structure
#'
#' @param n A scalar value indicating the number of subjects.
#' @param ni A list of length n containing the number of measurements per subject.
#' @param e A list of length n containing the residuals.
#' @param phi The estimate of variance parameter.
#'
#' @return Method of moments estimate for the correlation parameter alpha under unstructured correlation structure
#' @export
corr.un.fun <- function(n, ni, e, phi){
  corr.i <- lapply(1:n, function(i) (1/phi)*e[[i]]%*%t(e[[i]]))
  corr.i.upper <- lapply(corr.i, function(x) cbind(x[upper.tri(x)], 1:length(x[upper.tri(x)])))
  allMat <- as.data.frame(do.call("rbind", corr.i.upper))
  names(allMat) <- c("corr", "pos")
  corr.pars <- tapply(allMat$corr, allMat$pos, mean)
  return(corr.pars)
}

## Creates correlation matrix from correlation parameters under a working correlation structure
#' Function to construct the correlation matrix for a specific subject
#'
#' @param alpha The value of correlation parameter (a scalar value for "exchangeable" or "AR(1)" 
#' structure, or a vector with the elements of the upper triangular part for "unstructured").
#' @param ni A scalar value showing the number of measurements for a subject.
#' @param corstr The working correlation structure.
#'
#' @return The correlation matrix of dimension ni*ni.
#' @export

corMatF <- function(alpha, ni, corstr){
  switch(corstr,
         "exchangeable" = toeplitz(c(1, rep(alpha, ni-1))),
         "ar1" = toeplitz(alpha^(0:(ni-1))),
         "unstructured" = vec2uMat(alpha=alpha[1:(ni*(ni-1)/2)], ni),
         "independence" = diag(ni),
  )
}

## Creates correlation matrix from correlation parameters under unstructured correlation structure
#' Function to construct correlation matrix from correlation parameters for unstructured correlation structure
#'
#' @param alpha A vector with the elements of the upper triangular part for "unstructured".
#' @param ni A scalar value showing the number of measurements for a subject.
#'
#' @return An unstructured correlation matrix
#' @export

vec2uMat <- function(alpha, ni){
  x <- matrix(1, ni, ni)
  x[upper.tri(x)] <- alpha
  x[lower.tri(x)] <- t(x)[lower.tri(x)]
  return(x)
}

### Derivative function of SCAD penalty
#' Function to construct the derivative function of SCAD penalty
#'
#' @param lambda A scalar value of the tuning parameter.
#' @param par A scalar value for a single regression parameter estimate
#'
#' @return The derivative function of SCAD penalty
#' @export

q.scad <- function(lambda, par){
  b <- 3.7
  term1 <- 1
  term2 <- ((b*lambda-par)*as.numeric((b*lambda-par)>0))/((b-1)*lambda)
  out <- lambda*(term1*as.numeric(par <= lambda)+term2*as.numeric(par > lambda))
  return(out)
}

### Derivative function of MCP penalty
#' Function to construct the derivative function of MCP penalty
#'
#' @param lambda A scalar value of the tuning parameter.
#' @param par A scalar value for a single regression parameter estimate
#'
#' @return The derivative function of MCP penalty
#' @export
q.mcp <- function(lambda, par){
  b <- 3.7
  out <- ifelse(abs(par) <= b*lambda, (lambda - abs(par)/b)*sign(par), 0)
  return(out)
}

## Calculate DRIC for a single value of lambda (tuning parameter)
#' Function to calculate the value of Double Robust Information Criterion (DRIC)
#'
#' @param out.penalized.G A list corresponding to the output of penalized G-estimation for a 
#' single value of the tuning parameter. This is the output of the function penG().
#' @param n A scalar value indicating the number of subjects.
#' @param ni A list of length n containing the number of measurements per subject.
#' @param p A scalar value showing the total number of main effects to estimate (including the intercept).
#' @param y A list of length n, where each item in the list is a vector of dimension ni*1 containing the responses from a subject.
#' @param a A list of length n, where each item in the list is a vector of dimension ni*1 containing the exposure status of a subject.
#' @param E.a A list of length n, where each item in the list is a vector of dimension ni*1 containing the propensity scores of a subject.
#' @param l A list of length n, where each item in the list is a matrix of dimension ni*p containing the covariate values for a subject.
#' @param wc.str A working correlation structure.
#'
#' @return The value of DRIC under a working correlation structure
#' @export

calcDRIC <- function(out.penalized.G, n, ni, p, y, a, E.a, l, wc.str){
  estimate <- out.penalized.G$estimate
  sigma2.hat <- out.penalized.G$phi
  alpha.hat <- out.penalized.G$alpha
  En.mat <- out.penalized.G$En.mat
  V <- lapply(1:n, function(i) sigma2.hat*corMatF(alpha = alpha.hat, ni = ni[[i]], corstr = wc.str))
  V.inv <- lapply(1:n, function(i) solve(V[[i]]))
  e <- lapply(1:n, function(i)
    y[[i]]-(c(a[[i]])*l[[i]])%*%estimate[1:p]-l[[i]]%*%estimate[(p+1):(2*p)])

  sum3.n <- Reduce("+", lapply(1:n, function(i)
    t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%(c(a[[i]])*l[[i]])))
  sum4.n <- Reduce("+", lapply(1:n, function(i) t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%l[[i]]))
  sum6.n <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%V.inv[[i]]%*%l[[i]]))
  sum5.n <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%V.inv[[i]]%*%(c(a[[i]])*l[[i]])))

  df.lambda <- sum(diag(solve(rbind(cbind(sum3.n, sum4.n), cbind(sum5.n, sum6.n)) + n*En.mat)%*%
                          rbind(cbind(sum3.n, sum4.n), cbind(sum5.n, sum6.n))))

  sum.DRIC <- Reduce("+", lapply(1:n, function(i) sum((c(abs(a[[i]]-E.a[[i]]))*e[[i]]^2))))
  DRIC <- log(sum.DRIC/sum(unlist(ni))) + log(log(n))*log(2*p)*df.lambda/n ##for high-dimension
  # DRIC <- log(sum.DRIC/sum(unlist(ni))) + log(n)*df.lambda/n ##works for low-dimension
  return(DRIC)
}
