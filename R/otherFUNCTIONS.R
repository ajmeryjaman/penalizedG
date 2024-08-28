## The function performing penalized G-estimation for a given data, working
## correlation structure and fixed tuning parameter

penG <- function(data, wc.str, tf.model, treat.model, lambda, maxitr = maxitr){
  # n: number of patients
  # p: number of variables in the treatment-free model including the intercept
  # ni, l, a, y, l.mat.split, e: all are list of length n
  # estimate.current: a vector of length 2*p

  n <- length(split(data, data$id))
  dat <- split(data, data$id)
  l.mat <- data.frame(id=data$id,model.matrix(tf.model,data)) # cov+treat history
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


  ## initial estimate at iteration = 0
  estimate.current <- solve(cbind(rbind(sum4,sum2),rbind(sum6,sum1)))%*%rbind(sum5,sum3)
  e <- lapply(1:n, function(i)
    y[[i]]-(c(a[[i]])*l[[i]])%*%estimate.current[1:p]-l[[i]]%*%estimate.current[(p+1):(2*p)])
  sigma2.hat <- phi.hat.fun(n, ni, e)

  if(wc.str == "unstructured"){
    sum.cov <- Reduce("+", lapply(e, function(res.i) res.i%*%t(res.i)))
    cov.un <- sum.cov/n
  }
  cov.hat <- switch(wc.str,
                    "exchangeable" = cov.exch.fun(n, ni, e),
                    "ar1" = cov.ar1.fun(n, ni, e),
                    "independence" = 0,
                    "unstructured" = cov.un)
  # alpha.hat <- cov.hat/sigma2.hat
  # cov.hat/diag(cov.hat)
  itr <- 0

  ### iteration = 1

  ## corr.mat, V, V.inv >> list of length n
  V <- lapply(1:n, function(i)
    switch(wc.str,
           "exchangeable" = toeplitz(c(sigma2.hat, rep(cov.hat, ni[[i]]-1))),
           "ar1" = toeplitz(sigma2.hat*c(1,(cov.hat/sigma2.hat)^(1:(ni[[i]]-1)))),
           "independence" = diag(ni[[i]]),
           "unstructured" = cov.hat)
  )
  rcond <- lapply(V, function(x) rcond(x))
  if (sum(unlist(rcond) <= 1e-15) == 0){
    V.inv <- lapply(1:n, function(i) solve(V[[i]]))
    ## replacing a[[i]] with a[[i]]-E.a[[i]] in the last part of H(theta), which is (H_i A_i.H_i),
    ## in the algorithm helps achieve faster convergence, while the results are same
    sum1.n <- Reduce("+", lapply(1:n, function(i) t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%e[[i]]))
    sum2.n <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%V.inv[[i]]%*%e[[i]]))
    sum3.n <- Reduce("+", lapply(1:n, function(i)
      t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%(c(a[[i]]-E.a[[i]])*l[[i]])))
    sum4.n <- Reduce("+", lapply(1:n, function(i) t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%l[[i]]))
    sum6.n <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%V.inv[[i]]%*%l[[i]]))
    sum5.n <- t(sum4.n)

    En.mat <- diag(c(0,q.scad(lambda,abs(estimate.current[2:p]))/(10^(-6)+abs(estimate.current[2:p])),
                     rep(0,p)))

    estimate.new <- estimate.current + solve(rbind(cbind(sum3.n,sum4.n),cbind(sum5.n,sum6.n))+n*En.mat)%*%(rbind(sum1.n,sum2.n)-n*En.mat%*%estimate.current)
    e <- lapply(1:n, function(i)
      y[[i]]-(c(a[[i]])*l[[i]])%*%estimate.new[1:p]-l[[i]]%*%estimate.new[(p+1):(2*p)])
    sigma2.hat <- phi.hat.fun(n, ni, e)
    if(wc.str == "unstructured"){
      sum.cov <- Reduce("+", lapply(e, function(res.i) res.i%*%t(res.i)))
      cov.un <- sum.cov/n
    }
    cov.hat <- switch(wc.str,
                      "exchangeable" = cov.exch.fun(n, ni, e),
                      "ar1" = cov.ar1.fun(n, ni, e),
                      "independence" = 0,
                      "unstructured" = cov.un)
    # alpha.hat <- cov.hat/sigma2.hat
    itr <- 0+1

    ## next iterations 2, 3, ...
    while(max(abs(estimate.new-estimate.current)) > 0.0001 & itr<= maxitr &
          sum(unlist(rcond) <= 1e-15) == 0){ ### itr -> [2, maxitr]

      estimate.current <- estimate.new

      En.mat <- diag(c(0,q.scad(lambda,abs(estimate.current[2:p]))/(10^(-6)+abs(estimate.current[2:p])),
                       rep(0,p)))
      ## corr.mat, V, V.inv >> list of length n
      V <- lapply(1:n, function(i)
        switch(wc.str,
               "exchangeable" = toeplitz(c(sigma2.hat, rep(cov.hat, ni[[i]]-1))),
               "ar1" = toeplitz(sigma2.hat*c(1,(cov.hat/sigma2.hat)^(1:(ni[[i]]-1)))),
               "independence" = diag(ni[[i]]),
               "unstructured" = cov.hat)
      )
      rcond <- lapply(V, function(x) rcond(x))
      if (sum(unlist(rcond) <= 2.220446e-16) == 0){
        V.inv <- lapply(1:n, function(i) solve(V[[i]]))

        sum1.n <- Reduce("+", lapply(1:n, function(i) t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%e[[i]]))
        sum2.n <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%V.inv[[i]]%*%e[[i]]))
        sum3.n <- Reduce("+", lapply(1:n, function(i)
          t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%(c(a[[i]]-E.a[[i]])*l[[i]])))
        sum4.n <- Reduce("+", lapply(1:n, function(i) t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%l[[i]]))
        sum6.n <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%V.inv[[i]]%*%l[[i]]))
        sum5.n <- t(sum4.n)

        estimate.new <- estimate.current + solve(rbind(cbind(sum3.n,sum4.n),cbind(sum5.n,sum6.n))+n*En.mat)%*%(rbind(sum1.n,sum2.n)-n*En.mat%*%estimate.current)
        e <- lapply(1:n, function(i)
          y[[i]]-(c(a[[i]])*l[[i]])%*%estimate.new[1:p]-l[[i]]%*%estimate.new[(p+1):(2*p)])
        sigma2.hat <- phi.hat.fun(n, ni, e)
        if(wc.str == "unstructured"){
          sum.cov <- Reduce("+", lapply(e, function(res.i) res.i%*%t(res.i)))
          cov.un <- sum.cov/n
        }
        cov.hat <- switch(wc.str,
                          "exchangeable" = cov.exch.fun(n, ni, e),
                          "ar1" = cov.ar1.fun(n, ni, e),
                          "independence" = 0,
                          "unstructured" = cov.un)
        # alpha.hat <- cov.hat/sigma2.hat
        itr <- itr+1
      }
    }
  }

  ##Asymptotic (sandwich) variance calculation
  if (itr > 0 & itr <= maxitr & sum(unlist(rcond) <= 1e-15) == 0){
    l.mat.treat <- data.frame(id=data$id, model.matrix(treat.model,data)) # cov+treat history
    l.mat.treat.split <- split(l.mat.treat, l.mat.treat$id)
    l.treat <- lapply(1:n, function(i) as.matrix(l.mat.treat.split[[i]][,-1]))
    V <- lapply(1:n, function(i)
      switch(wc.str,
             "exchangeable" = toeplitz(c(sigma2.hat, rep(cov.hat, ni[[i]]-1))),
             "ar1" = toeplitz(sigma2.hat*c(1,(cov.hat/sigma2.hat)^(1:(ni[[i]]-1)))),
             "independence" = diag(ni[[i]]),
             "unstructured" = cov.hat)
    )
    V.inv <- lapply(1:n, function(i) solve(V[[i]]))
    sum3.n <- Reduce("+", lapply(1:n, function(i)
      t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%(c(a[[i]]-E.a[[i]])*l[[i]])))
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
    return(list(estimate = estimate.new, phi = sigma2.hat, cov = cov.hat,
                asymp.var.psi = var.psi.hat, data = data, nitr = itr, error = 0))
  } else {
    return(list(error = 1, nitr = itr)) #Estimation did not converge
  }
}

### Other functions to be used

## Method of moments estimate for the variance parameter sigma^2
phi.hat.fun <- function(n, ni, e){
  sum.phi.i <- lapply(1:n, function(i) (1/ni[[i]])*sum((e[[i]])^2))
  sigma2.hat <- sum(unlist(sum.phi.i))/n
  return(sigma2.hat)
}

## Method of moments estimate of the covariance under exchangeable structure
cov.exch.fun <- function(n, ni, e){
  sum.corr.i <- lapply(1:n, function(i) (1/(ni[[i]]*(ni[[i]]-1)))*(sum(e[[i]]%*%t(e[[i]]))
                                                                   -sum(diag(e[[i]]%*%t(e[[i]])))))
  cov <- sum(unlist(sum.corr.i))/n
  return(cov) ####
}

## Method of moments estimate of the covariance under AR(1) structure
cov.ar1.fun <- function(n, ni, e){
  sum.corr.i <- lapply(1:n, function(i) (1/(ni[[i]]-1))*(t(e[[i]][1:(ni[[i]]-1)])%*%(e[[i]][2:ni[[i]]])))
  cov <- sum(unlist(sum.corr.i))/n
  return(cov) ####
}

### Derivative function of the SCAD penalty
q.scad <- function(lambda,par){
  b <- 3.7
  term1 <- 1
  term2 <- ((b*lambda-par)*as.numeric((b*lambda-par)>0))/((b-1)*lambda)
  out <- lambda*(term1*as.numeric(par <= lambda)+term2*as.numeric(par > lambda))
  return(out)
}

## calculate BIC for a single value of the tuning parameter lambda
calcBIC <- function(out.penalized.G, n, ni, p, y, a, E.a, l, wc.str, lambda.seq){
  estimate <- out.penalized.G$estimate
  sigma2.hat <- out.penalized.G$phi
  cov.hat <- out.penalized.G$cov
  V <- lapply(1:n, function(i)
    switch(wc.str,
           "exchangeable" = toeplitz(c(sigma2.hat, rep(cov.hat, ni[[i]]-1))),
           "ar1" = toeplitz(sigma2.hat*c(1,(cov.hat/sigma2.hat)^(1:(ni[[i]]-1)))),
           "independence" = diag(ni[[i]]),
           "unstructured" = cov.hat))
  V.inv <- lapply(1:n, function(i) solve(V[[i]]))
  e <- lapply(1:n, function(i)
    y[[i]]-(c(a[[i]])*l[[i]])%*%estimate[1:p]-l[[i]]%*%estimate[(p+1):(2*p)])

  En.mat <- diag(c(0,q.scad(lambda.seq,abs(estimate[2:p]))/(10^(-6)+abs(estimate[2:p])), rep(0,p)))
  sum3.n <- Reduce("+", lapply(1:n, function(i)
    t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%(c(a[[i]]-E.a[[i]])*l[[i]])))
  sum4.n <- Reduce("+", lapply(1:n, function(i) t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%l[[i]]))
  sum6.n <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%V.inv[[i]]%*%l[[i]]))
  sum5.n <- t(sum4.n)
  df.lambda <- sum(diag(solve(rbind(cbind(sum3.n,sum4.n),cbind(sum5.n,sum6.n))+n*En.mat)%*%
                          rbind(cbind(sum3.n,sum4.n),cbind(sum5.n,sum6.n))))

  sum.bic <- Reduce("+", lapply(1:n, function(i) t(e[[i]])%*%e[[i]]))
  bic <- log(sum.bic/sum(unlist(ni))) + df.lambda*log(sum(unlist(ni)))/sum(unlist(ni))
  return(bic)
}

