#' Penalized G-estimation
#' @description This function performs penalized G-estimation for a given longitudinal data,
#' a specific working correlation structure and a sequence of tuning parameters, and returns result
#' under the optimal value (selected by a double-robust information  criterion) of the tuning parameter in the given range.
#' Currently, the function allows a continuous outcome and a binary treatment/exposure. The outcome, the exposure and the potential
#' confounders, all can be time-varying, but the potential confounders should be continuous or binary.
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
#' @param lambda.seq A sequence of values (in decreasing order) for the tuning parameter. Typically,
#' the first element in this sequence is a lowest positive number such that all the effect modifiers (EMs) are eliminated by the method, while
#' the last element is a value close to zero for which none of the EMs are eliminated.
#' @param maxitr Maximum number of iterations allowed. The default value is 100.
#' @param penalty The penalty type to be used, available options include "SCAD" and "MCP". The default choice is "SCAD".
#' "SCAD" refers to the Smoothly Clipped Absolute Deviation penalty and "MCP" refers to the Minimax Concave Penalty.
#'
#' @return A list containing the following:
#' \item{estimate}{The vector of parameter estimates. First half corresponds to the blip coefficients
#' and the second half corresponds to the coefficients of the treatment-free model}
#' \item{Selected.EMs}{A vector showing which variables are selected as effect modifiers (EMs)}
#' \item{sigma2.hat}{The estimated variance parameter sigma^2.}
#' \item{alpha.hat}{The estimated correlation parameter(s) alpha(s) if the provided structure is either
#' "exchangeable", "ar1, or "unstructured". For unstructured, the elements of alpha.hat correspond
#' to the upper triangular portion of working correlation matrix having dimension equal to
#' the largest cluster size.}
#' \item{asymp.var.psi}{The sandwich variance-covariance matrix of the blip coefficients. Although, the
#' function will provided
#' sandwich variance, a data analyst should note that tests or confidence intervals based on this
#' sandwich estimates are expected to exhibit
#' inflated type I errors in finite samples.}
#' \item{nitr}{The number of iterations at which the estimation converged with the optimal tuning parameter.}
#' \item{lambda.optimal}{The optimal tuning parameter from the given sequence (lambda.seq).}
#' \item{data}{The original dataframe with an additional column that contains the estimated
#' propensity scores.}
#'
#' @export
#' 
#' @importFrom stats model.matrix as.formula binomial glm predict.glm toeplitz
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
#'   delta <- c(1, 1, 1.2, 1.2, -0.9, 0.8, -1, 1) # treatment-free model parameters
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
#' # Generate the data
#' data.s <- data.gen(n = 500, ni = rep(6, 500), sigma2.e = 1, alpha = 0.8,
#'                    corstr = "exchangeable", autocorr.coef = 0)
#'
#' # Treatment-free model is misspecified, because we did not include exp(l5) in it
#' tf.model <- as.formula(paste("~",paste(c(paste("l",1:6,sep=""), paste("x",1:10,sep="")),
#'                                                    collapse = "+"), collapse=""))
#' # Treatment model is correctly specified                                                 
#' treat.model <- ~l1+l2+l3+l4+l5+l6
#' p <- length(model.matrix(tf.model, data=data.s)[1,])
#' 
#' # Creating a sequence of tuning parameter values in decreasing order
#' lam_max <- 1
#' lam_min <- 0.01*lam_max
#' lambda.seq <- sort(seq(lam_min,lam_max,(lam_max-lam_min)/99), decreasing=TRUE)
#'
#' # Perfom the penalized G-estimation for all lambda values and get the fit corresponding
#' # to the optimal lambda chosen by a Doubly Robust Information Criterion (DRIC)
#' out <- penalizedG(data = data.s, wc.str = "exchangeable", id.var="id", response.var="y",
#'                    treat.var="a", tf.model=tf.model, treat.model = treat.model,
#'                    lambda.seq = lambda.seq, maxitr = 50, penalty = "SCAD")
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
#' 
## The main function performing penalized G-estimation for a given data, working
## correlation structure and a sequence of tuning parameter values
penalizedG <- function(data, wc.str, id.var, response.var, treat.var, tf.model, treat.model,
                        lambda.seq, maxitr = 100, penalty = "SCAD"){
  if(!wc.str %in% c("independence", "exchangeable", "ar1", "unstructured")){
    stop("Working correlation structure must be one among independence,
         exchangeable, ar1, and unstructured\n")
  }
  data.org <- data
  if(any(c(is.na(data.org)))){
    stop("Data can not contain any missing values\n")
  }

  names(data)[names(data)==id.var] <- "id"
  names(data)[names(data)==treat.var] <- "a"
  names(data)[names(data)==response.var] <- "y"

  if(length(levels(as.factor(data$a))) != 2){
    stop("The treatment/exposure variable must be binary\n")
  }

  ## Next we calculate the propensity scores from pooled data
  treat.mod.formula <- as.formula(paste("a",paste(treat.model, collapse = ""), collapse=""))
  data$E.a <- predict.glm(glm(treat.mod.formula,family=binomial, data=data), type = "response")
  data.org$E.a <- data$E.a

  ##Next we perform penalized G-estimation for a sequence of tuning parameters and
  ##we record if there is any error (i.e., the estimation did not converge)
  out.penG <- lapply(lambda.seq, function(k)
    penG(data=data, wc.str = wc.str, tf.model = tf.model, treat.model = treat.model,
         lambda = k, maxitr = maxitr, penalty = penalty))
  errors <- unlist(lapply(out.penG, function(x) x$error))

  ##Next we split the data and construct required quantities as a list of length n
  ## for computing the values of double-robust information criterion (DRIC)
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

  ###DRIC calculation for the tuning parameters for which the estimation converged (i.e., error = 0)
  DRIC <- lapply(which(errors == 0), function(k)
    calcDRIC(out.penG[[k]], n, ni, p, y, a, E.a, l, wc.str = wc.str)
  )
  lambda.no.error <- lambda.seq[which(errors == 0)]
  lambda.selected <- lambda.no.error[which.min(unlist(DRIC))]
  res <- out.penG[[which(lambda.seq==lambda.selected)]]

  estimate <- res$estimate
  sigma2.hat <- res$phi
  alpha.hat <- ifelse(wc.str %in% c("exchangeable", "ar1", "unstructured"), res$alpha, NA)
  asymp.var.psi <- res$asymp.var.psi
  nitr <- res$nitr

  row.names(estimate) <- c(treat.var, paste(treat.var, "*", colnames(l.mat)[-c(1,2)], sep=""),
                           "Intercept", colnames(l.mat)[-c(1,2)])
  row.names(asymp.var.psi) <- c(treat.var, paste(treat.var, "*", colnames(l.mat)[-c(1,2)], sep=""))
  colnames(asymp.var.psi) <- c(treat.var, paste(treat.var, "*", colnames(l.mat)[-c(1,2)], sep=""))
  all_vars <- names(model.matrix(tf.model, data=data)[1, -1])
  selected.EMs <- all_vars[abs(res$estimate[2:p]) > 0.001]

  return(list(estimate = estimate, selected.EMs = selected.EMs, sigma2.hat = sigma2.hat,
              alpha.hat = alpha.hat, asymp.var.psi = asymp.var.psi,
              nitr=nitr, lambda.optimal = lambda.selected, data = data.org))
}
