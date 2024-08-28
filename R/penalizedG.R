#' Penalized G-estimation
#' @description This function performs penalized G-estimation for a given longitudinal data,
#' a specific working correlation structure and a sequence of tuning parameters, and returns result
#' under the optimal value (selected by a BIC criterion) of the tuning parameter in the given range.
#' Currently, the function allows only a continuous outcome and a binary treatment/exposure.
#' @param data A data frame containing the variables in longitudinal format.
#' @param wc.str A character string specifying the working correlation structure. The
#' following are currently allowed: "independence", "exchangeable", "ar1", and "unstructured".
#' @param id.var The column name in data that corresponds to the variable id (unique identifier).
#' @param response.var The column name in data that corresponds to the response variable.
#' @param treat.var The column name in data that corresponds to the treatment/exposure variable.
#' @param tf.model A single formula object specifying the covariates of a (linear)
#' treatment-free model.
#' @param treat.model A single formula object specifying the covariates of the logistic
#' regression model for the binary treatment.
#' @param lambda.seq A sequence  values for the tuning parameter.
#' @param maxitr Maximum number of iterations allowed. The default value is 50.
#'
#' @return A list containing the following:
#' \item{estimate}{The vector of parameter estimates. First half corresponds to the blip coefficients and the second half corresponds to the coefficients of the treatment-free model}
#' \item{Selected.EMs}{A vector showing which variables are selected as effect modifiers (EMs)}
#' \item{sigma2.hat}{The estimated variance parameter sigma^2.}
#' \item{alpha.hat}{ The estimated correlation parameter alpha if the provided structure is either "exchangeable" or "ar1.}
#' \item{cov.hat.un}{The estimated covariance matrix if the provided structure is "unstructured".}
#' \item{sandwich.var.psi}{The sandwich variance-covariance matrix of the blip coefficients.}
#' \item{nitr}{The number of iterations at which the estimation converged with the optimal tuning parameter.}
#' \item{lambda.optimal}{The optimal tuning parameter from the given sequence (lambda.seq).}
#' \item{data}{The original dataframe same as provided argument}
#'
#' @export penalizedG
#'
#' @examples
#' library(mvtnorm)
#' expit <- function(x) exp(x)/(1+exp(x))
#'
#' ## data.gen is a function that generates a longitudinal data set for a specific correlation
#' ## structure. Available structures in this function are: independence, exchangeable and ar1.
#'
#' # Arguments(data.gen):
#' #     n = Number of subjects
#' #     ni = A vector containing number of time points for each subject
#' #     sigma2.e = Error variance
#' #     alpha = Correlation parameter
#' #     corstr = The correlation structure among the repeated outcomes
#' #     autocorr.coef = The autocorrelation coefficient for inducing correlation among the
#' #     continuous confounders and the noise covariates
#'
#' data.gen <- function(n, ni, sigma2.e, alpha, corstr, autocorr.coef){
#'   # treatment model: a ~ 1+l1+...+l6 [logit link]
#'   # treatment-free model:  ~ 1+l1+...+l6+exp(l5)
#'   # blip model: a * (1+l1+...+l5)
#'   # noise covariates: x1,...,x10
#'   beta <- c(0, 1, -1.1, 1.2, 0.75, -0.9, 1.2) # treatment model parameters
#'   delta <- c(1, 1, 1.5, 1.2, -0.9, 0.8, -1, 1) # treatment-free model parameters
#'   psi <- c(1, 1, -1, -0.9, 0.8, 1, 0) # blip parameters
#'
#'   # generating two continuous baseline covariates
#'   l1 <- rnorm(n, 0, 1)
#'   l2 <- rnorm(n, 0, 1)
#'
#'   # V is the covariance matrix of the time-varying confounders (l3,..,l6) and
#'   # noise covariates (x1,...,x10)
#'   V <- toeplitz(autocorr.coef^(0:(14-1)))
#'
#'   lx <- a <- y <- vector(mode="list", length=n)
#'   lx.mat <- NULL
#'   for(i in 1:n){
#'     a[[i]] <- y[[i]] <- rep(NA, ni[i])
#'     lx[[i]] <- matrix(NA, ni[i], 16)
#'     lx[[i]][,1] <- rep(l1[i], ni[i])
#'     lx[[i]][,2] <- rep(l2[i], ni[i])
#'
#'     corr.mat <- switch(corstr,
#'                        "exchangeable" = toeplitz(c(1, rep(alpha, ni[i]-1))),
#'                        "ar1" = toeplitz(alpha^(0:(ni[i]-1))),
#'                        "independence" = diag(ni[i])
#'     )
#'     cov.mat <- diag(sqrt(sigma2.e), ni[i]) %*% corr.mat %*% diag(sqrt(sigma2.e), ni[i])
#'     e <- rmvnorm(1, sigma = cov.mat)
#'
#'     ### j=1
#'     mu.l <- rep(0, 4) #rep(0.3*lx[[i]][1, 1]+0.3*lx[[i]][1, 2],4)
#'     mu.x <- rep(0,10)
#'     lx[[i]][1,3:16] <- mvtnorm::rmvnorm(1, mean=c(mu.l,mu.x), sigma = V)
#'     a[[i]][1] <- rbinom(1, 1, expit(sum(c(1,lx[[i]][1,1:6])*beta)))  ## no correlation
#'     tf.mean <- sum(c(1, lx[[i]][1,1:6], exp(lx[[i]][1,5])) * delta)
#'     blip <- sum(c(a[[i]][1], a[[i]][1]*c(lx[[i]][1,1:6])) * psi)
#'     y[[i]][1] <- (tf.mean + blip) + e[1]
#'
#'     ## j=2:ni
#'     for(j in 2:ni[i]){
#'       mu.l <- 0.3*lx[[i]][j-1, 3:6] + 0.3*a[[i]][j-1]
#'       mu.x <- 0.5*lx[[i]][j-1,7:16]
#'       lx[[i]][j,3:16] <- mvtnorm::rmvnorm(1, mean=c(mu.l,mu.x), sigma = V)
#'       a[[i]][j] <- rbinom(1, 1, expit(sum(c(1,lx[[i]][j,1:6])*beta)))
#'       tf.mean <- sum(c(1, lx[[i]][j,1:6], exp(lx[[i]][j,5])) * delta)
#'       blip <- sum(c(a[[i]][j], a[[i]][j]*c(lx[[i]][j,1:6])) * psi)
#'       y[[i]][j] <- (tf.mean + blip) + e[j]
#'     }
#'     lx.mat <- rbind(lx.mat, lx[[i]])
#'   }
#'
#'   colnames(lx.mat) <- c(paste("l", 1:6, sep=""), paste("x", 1:10, sep=""))
#'   data <- data.frame(id=rep(1:n, times=ni), a=unlist(a), lx.mat, y=round(unlist(y),3))
#'   return(data)
#' }
#'
#' data.s <- data.gen(n = 200, ni = rep(6, 200), sigma2.e = 1, alpha = 0.8,
#'                    corstr = "exchangeable", autocorr.coef = 0)
#'
#' #treatment-free model is misspecified, because we did not include exp(l5) in it
#' tf.model <- as.formula(paste("~",paste(c(paste("l",1:6,sep=""), paste("x",1:10,sep="")),
#'                                                    collapse = "+"), collapse=""))
#' treat.model <- ~l1+l2+l3+l4+l5+l6
#' p <- length(all.vars(tf.model))+1
#' maxitr <- 50
#'
#' lam_max <- 1
#' lam_min <- 0.01*lam_max
#' lambda.seq <- sort(seq(lam_min,lam_max,(lam_max-lam_min)/99), decreasing=TRUE)
#'
#' out <- penalized.G(data = data.s, wc.str = "exchangeable", id.var="id", response.var="y",
#'                    treat.var="a", tf.model=tf.model, treat.model = treat.model,
#'                    lambda.seq = lambda.seq, maxitr = 50)
#'
#' names(out)
#' out$lambda.optimal
#' out$nitr
#' out$selected.EMs
#'
#' effect.modification <- out$estimate[2:p]
#' names(effect.modification) <- row.names(out$estimate)[2:p]
#' effect.modification <- effect.modification[abs(effect.modification) > 0.001]
#'
#' main.effect <- out$estimate[1]
#'
#' round(main.effect, 2)
#' round(effect.modification, 2)


## The main function performing penalized G-estimation for a given data, working
## correlation structure and a sequence of tuning parameter values

penalized.G <- function(data, wc.str, id.var, response.var, treat.var, tf.model, treat.model,
         lambda.seq, maxitr = 50){
  if(!wc.str %in% c("independence", "exchangeable", "ar1", "unstructured")){
    stop("Working correlation structure must be one among independence,
         exchangeable, ar1, and unstructured\n")
  }
  data.org <- data
  if(any(c(is.na(data.org)))){
    stop("There can not be any missing values in the data\n")
  }

  names(data)[names(data)==id.var] <- "id"
  names(data)[names(data)==treat.var] <- "a"
  names(data)[names(data)==response.var] <- "y"

  if(length(levels(as.factor(data$a))) != 2){
    stop("The treatment/exposure variable must be binary\n")
  }

  ## Next we calculate the propensity scores
  treat.mod.formula <- as.formula(paste("a",paste(treat.model,collapse = ""),collapse=""))
  data$E.a <- predict.glm(glm(treat.mod.formula,family=binomial, data=data), type = "response")

  ##Next we perform penalized G-estimation for a sequence of tuning parameters and
  ##we record if there is any error (i.e., the estimation did not converge)
  out.penG <- lapply(lambda.seq, function(k)
    penG(data=data, wc.str = wc.str, tf.model = tf.model, treat.model = treat.model,
         lambda = k, maxitr = maxitr))
  errors <- unlist(lapply(out.penG, function(x) x$error))

  ##Next we split the data and construct required quantities as a list of length n
  ## for computing the BIC values
  dat <- split(data, data$id)
  l.mat <- data.frame(id=data$id, model.matrix(tf.model, data)) # cov+treat history
  l.mat.split <- split(l.mat, l.mat$id)
  n <- length(unique(data$id))
  ni <- lapply(1:n, function(i) length(dat[[i]]$id))
  l <- lapply(1:n, function(i) as.matrix(l.mat.split[[i]][,-1]))
  a <- lapply(1:n, function(i) as.matrix(dat[[i]]$a))
  y <- lapply(1:n, function(i) as.matrix(dat[[i]]$y))
  E.a <- lapply(1:n, function(i) as.matrix(dat[[i]]$E.a))
  p <- dim(l.mat[,-1])[2] # dimension of psi/delta including intercept

  ###BIC calculation for the tuning parameters for which the estimation converged (i.e., error = 0)
  bic <- lapply(which(errors == 0), function(k)
    calcBIC(out.penG[[k]], n, ni, p, y, a, E.a, l, wc.str = wc.str, lambda.seq[k])
  )
  lambda.no.error <- lambda.seq[which(errors == 0)]
  lambda.selected <- lambda.no.error[which.min(unlist(bic))]
  res <- out.penG[[which(lambda.seq==lambda.selected)]]
  estimate <- res$estimate
  sigma2.hat <- res$phi
  alpha.hat <- ifelse(wc.str %in% c("exchangeable", "ar1"), res$cov/sigma2.hat, NA)
  cov.hat.un <- ifelse(wc.str == "unstructured", res$cov, NA)
  asymp.var.psi <- res$asymp.var.psi
  nitr <- res$nitr

  row.names(estimate) <- c(treat.var, paste(treat.var, "*", all.vars(tf.model), sep=""),
                           "Intercept", all.vars(tf.model))
  row.names(asymp.var.psi) <- c(treat.var, paste(treat.var, "*", all.vars(tf.model), sep=""))
  colnames(asymp.var.psi) <- c(treat.var, paste(treat.var, "*", all.vars(tf.model), sep=""))

  as.formula(paste("a",paste(treat.model,collapse = ""),collapse=""))
  c(paste("l", 1:6, sep=""), paste("x", 1:10, sep=""))
  selected.EMs <- all.vars(tf.model)[abs(res$estimate[2:p]) > 0.001]

  return(list(estimate = estimate, selected.EMs = selected.EMs, sigma2.hat = sigma2.hat,
              alpha.hat = alpha.hat, cov.hat.un = cov.hat.un, sandwich.var.psi = asymp.var.psi,
              nitr=nitr, lambda.optimal = lambda.selected, data = data.org))
}
