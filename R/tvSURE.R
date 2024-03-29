#' Time-Varying Seemingly Unrelated Regression Equations Model
#'
#' Fits a set of balanced linear structural equations using Time-varying Ordinary Least 
#' Squares (tvOLS), Time-varying Seemingly Unrelated Regression (tvGLS), when the error 
#' variance-covariance matrix is known, or Time-varying Feasible Seemingly Unrelated 
#' Regression (tvFGLS), when the error variance-covariance matrix is unknown.
#'
#'
#' This function wraps up the kernel smoothing "tvOLS" and "tvGLS" estimators. The former is used when
#' equations are considered independent while the later assumes that the error term is correlated
#' amongst equations. This relation is given in matrix "Sigma" which is used in the estimation. When
#' "Sigma" is known, the estimates are calculated via the "tvGLS", and via the "tvFGLS" when "Sigma"
#' is unknown and must be estimated.
#'
#' Bandwidth selection is of great importance in kernel smoothing methodologies and it is done
#' automatically by cross-validation. One important aspect in the current packages is that the
#' bandwidth is selected independently for each equation and then the average is taken to use the
#' same bandwidth for each equation. It has been shown in Casas et al. (2017) that using
#' different bandwidths for each equation is in general a bad practice, even for uncorrelated equations.
#' Even though, the user may be able to use different bandwidths calling functions \code{\link{bw}} and
#' \code{\link{tvGLS}} separatedly.
#'
#' A system consists of "neq" number of equations with "obs" number of observations each and a number of
#' variables not necessarily equal for all equations. The matrix notation is:
#' \deqn{Y_{t} = X_t \beta_{t}+u_{t}}
#' where \eqn{Y_t  = (y_{1t}, y_{2t}, \ldots, y_{neq t})'}, \eqn{X_t = diag (x_{1t}, x_{2t}, \ldots, x_{neq t})}
#' and \eqn{\beta_{t} = \left(\beta _{1t}', \ldots, \beta _{neq t}'\right)'} is a vector of order the
#' total number of variables in the system. The error vector \eqn{u_{t} = (u_{1t}, u_{2t}, \ldots, u_{neq t})'}
#' has zero mean and  covariance matrix \eqn{E(u_t u_t') = \Sigma_t}.
#'
#' @references
#' Casas, I., Ferreira, E., and Orbe, S. (2017) Time-Varying Coefficient Estimation 
#' in SURE Models: Application to Portfolio Management. Available at SSRN: 
#' https://ssrn.com/abstract=3043137
#' 
#' Chen, X. B., Gao, J., Li, D., and Silvapulle, P (2017) Nonparametric Estimation and 
#' Forecasting for Time-Varying Coefficient Realized Volatility Models.
#' \emph{Journal of Business and Economic Statistics}, pp.1-13
#'
#' Granger, C. W (2008) Non-Linear Models: Where Do We Go Next - Time Varying
#' Parameter Models? \emph{Studies in Nonlinear Dynamics and Econometrics}, 12, pp. 1-11.
#'
#' Kristensen, D (2012) Non-parametric detection and estimation of structural change.
#' \emph{Econometrics Journal}, 15, pp. 420-461.
#'
#' Orbe, S., Ferreira, E., and Rodriguez-Poo, J (2004) On the estimation and testing of
#' time varying constraints in econometric models, \emph{Statistica Sinica}.
#'
#'
#' @aliases tvsure-class tvsure
#' @rdname tvSURE
#' @param formula A list of formulas, one for each equation.
#' @param z A vector containing the smoothing variable.
#' @param ez (optional) A scalar or vector with the smoothing values. If 
#' values are not included then the vector \code{z} is used instead.
#' @param bw An optional scalar or vector of length the number of equations. It represents the bandwidth in
#' the estimation of trend coefficients. If NULL, it is selected by cross validation.
#' @param cv.block A positive scalar with the size of the block in leave one block out cross-validation.
#' By default 'cv.block = 0' meaning leave one out cross-validation.
#' @param data A matrix or data frame containing variables in the formula.
#' @param method A character, a matrix of dimensions neq x neq or an array of dimensions obs x neq x neq, where
#' \code{obs} is the number of observations and \code{neq} is the number of equations.
#' If method = \code{identity} or \code{tvOLS} (default) then the method used is a time-varying OLS.
#' If method is a matrix (constant over time) or an array, then the \code{tvGLS} is called.
#' If method = \code{tvFGLS}, then the covariance matrix is estimated nonparametrically and the
#' estimation of the system is done as a whole.
#' @param Sigma A matrix of dimensions neq x neq or an array of dimensions neq x neq x obs
#' (neq = number of equations, obs = number of observations). It represents
#' the covariance matrix of the error term. Only necessary for method \code{tvGLS}.
#' @param bw.cov An optional scalar. It represents the bandwidth in the nonparametric estimation of the
#' varying covariance matrix. If NULL, it is selected by cross validation.
#' @param est The nonparametric estimation method, one of "lc" (default) for linear constant or "ll" for local linear.
#' @param tkernel A character, either "Triweight" (default), "Epa" or "Gaussian" kernel function.
#' @param singular.ok	Logical. If FALSE, a singular model is an error.
#' @param R An optional nrest x nvar x neq (nrest =  number of restrictions, nvar = number of variables in each equation,
#' neq = number of equations).
#' @param r An optional vector of length the number of restrictions. By default it contains zeros.
#' @param control list of control parameters.  The default is constructed by
#' the function \code{\link{tvreg.control}}.  See the documentation of
#' \code{\link{tvreg.control}} for details.
#' @param ... Other parameters passed to specific methods.
#' @return \code{tvSURE} returns a list of the class \code{tvsure} containing the results of the whole system, results of the estimation
#' and confidence instervals if chosen.
#' The object of class \code{tvsure} have the following components:
#' \item{coefficients}{An array of dimension obs x nvar x neq (obs = number of observations, nvar = number of variables
#' in each equation, neq = number of equations in the system) with the time-varying coefficients estimates.}
#' \item{Lower}{If \code{level} non equal zero, an array of dimension obs x nvar x neq containing the confidence 
#' interval lower band.}
#' \item{Upper}{If \code{level} non equal zero, an array of dimension obs x nvar x neq containing the confidence 
#' interval upper band.}
#' \item{Sigma}{An array of dimension obs x neq x neq with the estimates of the errors covariance matrix.}
#' \item{fitted}{The fitted values.}
#' \item{residuals}{Estimation residuals.}
#' \item{x}{A list with the regressors data.}
#' \item{y}{A matrix with the dependent variable data.}
#' \item{z}{A vector with the smoothing variable.}
#' \item{ez}{A vector with the smoothing estimation values.}
#' \item{bw}{Bandwidth of mean estimation.}
#' \item{obs}{Integer specifying the number of observations in each equation (balanced sample).}
#' \item{neq}{Integer specifying the number of equations.}
#' \item{nvar}{Vector of integers specifying the number of variables in each equation.}
#' \item{method}{Estimation method.}
#' \item{est}{Nonparemtric estimation methodology.}
#' \item{tkernel}{Kernel type.}
#' \item{bw.cov}{Bandwidht of Sigma estimation.}
#' \item{level}{Confidence interval range.}
#' \item{runs}{Number of bootstrap replications.}
#' \item{tboot}{Type of bootstrap.}
#' \item{BOOT}{List with all bootstrap replications of \code{coefficients}, if done.}
#' \item{R}{Restrictions matrix.}
#' \item{r}{Restrictions vector.}
#' \item{formula}{Initial formula.}
#' 
#' @seealso \code{\link{bw}}, \code{\link{tvCov}},  \code{\link{tvVAR}}, 
#' \code{\link{confint}}, \code{\link{plot}}, \code{\link{print}} and \code{\link{summary}}
#' 
#' @examples
#' \dontrun{
#' data("Kmenta", package = "systemfit")
#' eqDemand <- consump ~ price + income
#' eqSupply <- consump ~ price + farmPrice + trend
#' system <- list(demand = eqDemand, supply = eqSupply)
#' eqSupply2 <- consump ~  price + farmPrice 
#' system2 <- list(demand = eqDemand, supply = eqSupply2)
#' 
#' ##OLS estimation of a system
#' OLS <- systemfit::systemfit(system, method = "OLS", data = Kmenta)
#' ##tvOLS estimation of a system with the local linear estimator
#' ##removing trend because it is included in the intercept changing over time
#' TVOLS <- tvSURE(system2, data = Kmenta,  est = "ll")
#' 
#' ##SUR/FGLS estimation
#' FGLS <- systemfit::systemfit(system, data = Kmenta, method = "SUR")
#' ##tvSURE estimation
#' TVFGLS <- tvSURE(system, data = Kmenta, method = "tvFGLS")
#' }
#'
#'@export
tvSURE <- function (formula, z = NULL, ez = NULL, bw = NULL, cv.block = 0, data,  
                    method = c("tvOLS", "tvFGLS", "tvGLS"), Sigma = NULL, 
                    est = c("lc", "ll"), tkernel = c("Triweight", "Epa", "Gaussian"),
                    bw.cov = NULL, singular.ok = TRUE, R = NULL, r = NULL,
                    control = tvreg.control(...), ...)
{
  is.data <- inherits(data, c("data.frame", "matrix"))
  if(!is.data)
    stop("\nArgument 'data' should be entered and it should be a 'matrix' or a 'data.frame'.\n")
  if(!inherits(data, c("data.frame")))
    data <- as.data.frame(data)
  if(!inherits(formula,  "list"))
    stop("\nArgument 'formula' must be a list of formulas. \n")
  if(!all(lapply(formula, class) == "formula"))
    stop("\nArgument 'formula' must contain only objects of class 'formula'")
  neq <- length(formula)
  if(neq < 2)
    stop("\nThe list 'formula' should contain at least two equations for multivariate analysis.\n")
  if(!is.null(Sigma))
    if(any(is.na(Sigma)))
      stop("\nNAs in Sigma.\n")
  method <- match.arg(method)
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  nvar <- numeric(neq)
  if(is.null(names(formula)))
  {
    eq.names <- paste0("eq", c(1:neq))
  }
  else
  {
    eq.names <- names(formula)
    if(sum(regexpr(" |_", eq.names) != -1) > 0)
      stop("\nEquation labels may not contain blanks (' ') or underscores ('_')")
  }
  results <- list()
  callNoDots <- match.call(expand.dots = FALSE)
  mf <- callNoDots[c(1, match(c("data"), names(callNoDots), 0L))]
  mf$na.action <- as.name("na.pass")
  mf[[1]] <- as.name("model.frame")
  y  <- NULL
  x  <- list()
  y.names <- NULL
  for(i in 1:neq)
  {
    mf.eq <- mf
    mf.eq$formula <- formula[[i]]
    eval.mf <-  eval(mf.eq, parent.frame())
    terms <- attr(eval.mf, "terms")
    y <- cbind(y, stats::model.extract(eval.mf, "response"))
    y.names <- c(y.names, formula[[i]][[2]])
    x[[i]] <- stats::model.matrix(terms, eval.mf)
    nvar[i] <- NCOL(x[[i]])
    if(is.null(colnames(x[[i]])))
      colnames(x[[i]]) <- paste0("X", i, 1:nvar[i])
  }
  names(x) <- eq.names
  colnames(y) <- y.names
  obs <- NROW(y)
  if(!is.null(R))
  {
    R <- as.matrix(R)
    if(NCOL(R) != sum(nvar))
      stop("\nWrong dimension of R, it should have as many columns as variables 
           in the whole system. \n")
    if (is.null(r))
      r <- rep(0, NROW(R))
    else if (length(r) == 1)
      r <- rep(r, NROW(R))
    else if (length(r) != NROW(R) & length(r) != 1)
      stop("\nWrong dimension of r, it should be as long as the number of 
           rows in R. \n")
  }
  if (method == "identity" | method == "tvOLS")
  {
    if (is.null(bw))
    {
      cat("Calculating regression bandwidth...\n")
      bw <- bw(x = x, y = y, z = z, cv.block = cv.block, est = est, tkernel = tkernel, 
               singular.ok = singular.ok)
      cat("bw = ", bw, "\n")
    }
    else
    {
      if (any(bw < 5/obs))
        stop("\nAt least one of your bw bandwidths is smaller than 5/obs, 
             please increase! \n")
      else if (any(is.na(bw)))
        stop("\nThe bandwidth cannot be a no number.\n")
    }
    result <- tvGLS(x = x, y = y, z = z, ez = ez, bw = bw, R = R, r = r, 
                    est = est, tkernel = tkernel)
    Sigma <- array(rep(crossprod(result$residuals)/ (obs - neq), obs), dim = c(neq, neq, obs))
  }
  else if(method == "tvFGLS")
  {
    if (is.null(bw))
    {
      cat("Calculating regression bandwidth...\n")
      bw <- bw(x = x, y = y, z = z, cv.block = cv.block, est = est, tkernel = tkernel, 
               singular.ok = singular.ok)
      cat("bw = ", bw, "\n")
    }
    else
    {
      if (any(bw<5/obs))
        stop("\nAt least one of your bw bandwidths is smaller than 5/obs, please increase! \n")
      else if (any (is.na(bw)))
        stop("\nThe bandwidth cannot be a no number.\n")
    }
    result <- tvGLS(x = x, y = y, z = z, ez = ez, bw = bw, R = R, r = r, est = est, tkernel = tkernel)
    if(is.null(bw.cov))
    {
      cat("Calculating variance-covariance estimation bandwidth...\n")
      bw.cov <- bwCov(x = result$residuals, z = z, cv.block = cv.block, tkernel = tkernel)
      cat("bw = ", bw.cov, "\n")
    }
    cat("\n Sigma estimation...")
    Sigma <- tvCov(x = result$residuals, z = z,  bw = bw.cov, tkernel = tkernel)
    cat("\n FGLS estimation until convergence...")
    cat("\n Step 1, control maxiter ", control$maxiter, " control tolerance, ", control$tol)
    result <- tvGLS(x = x, y = y, z = z, ez = ez, bw = bw, Sigma = Sigma, R = R, r = r,
                    est = est, tkernel = tkernel)
    
    itertemp <- 1
    tol <- control$tol
    maxiter <- control$maxiter
    tolrel <- 1
    while((tolrel>tol) && (itertemp < maxiter))
    {
      cat("\n Step ", itertemp + 1, " tol: ", abs(tolrel))
      Sigma <- try(tvCov(x = result$residuals, z = z, bw = bw.cov, tkernel = tkernel))
      if (!inherits(Sigma, "array"))
      {
        bw.cov <- bwCov(x = result$residuals, z = z, cv.block = cv.block, tkernel = tkernel)
        Sigma <- try(tvCov(x = result$residuals, z = z, bw = bw.cov, tkernel = tkernel))
      }
      temp <- tvGLS(x = x, y = y, z = z, ez = ez, bw = bw, Sigma = Sigma, R = R, r = r,
                    est = est, tkernel = tkernel)
      tolrel <- mean(abs((result$coefficients - temp$coefficients)/result$coefficients))
      result <- temp
      itertemp <- itertemp + 1
    }
  }
  else if(method == "tvGLS")
  {
    if(is.matrix(Sigma))
    {
      if(NCOL(Sigma) != neq | NROW(Sigma) != neq)
        stop("\nWrong dimensions of Sigma. \n.")
      Sigma2 <- array(0, dim = c(neq, neq, obs))
      for (t in 1:obs)
        Sigma2[, , t] <- Sigma
      Sigma <- Sigma2
    }
    else if (is.array(Sigma))
    {
      dimensions <- dim(Sigma)
      if(dimensions[3] != obs | dimensions[2] != neq | dimensions[1] != neq)
        stop("\nWrong dimensions of Sigma. \n.")
    }
    else
      stop("\nSigma must be a matrix of dimensions neq x neq or an array of dimensions
           neq x neq x obs. \n")
    if (is.null(bw))
    {
      cat("Calculating regression bandwidth...\n")
      bw <- bw(x = x, y = y, z = z, Sigma = Sigma, est = est, tkernel = tkernel)
      cat("bw = ", bw, "\n")
      
    }
    else
    {
      if (any(bw < 5/obs))
        stop("\nAt least one of your bw bandwidths is smaller than 5/obs,
             please increase! \n")
      else if (any (is.na(bw)))
        stop("\nThe bandwidth cannot be a no number.\n")
    }
    result <- tvGLS(x = x, y = y, z = z, ez = ez, bw = mean(bw), Sigma = Sigma, R = R, r = r,
                    est = est, tkernel = tkernel)
  }
  coefficients <- result$coefficients
  resid <- result$residuals
  fitted <- result$fitted
  if(length(bw) == 1)
    names(bw) <- "bw.mean"
  else
    names(bw) <- paste("bw.", eq.names, sep = "")
  colnames(resid) <- eq.names
  colnames(fitted) <- eq.names
  var.names <- NULL
  for(i in 1:neq)
    var.names <- c(var.names, paste(colnames(x[[i]]), ".", eq.names[i], sep = ""))
  colnames(coefficients) <- var.names
  result <- list(coefficients =  coefficients, Lower = NULL, Upper = NULL, Sigma = Sigma,
                 fitted = fitted, residuals = resid, x = x, y = y, z = z, ez = ez,
                 bw = bw, cv.block = cv.block, obs = obs, neq = neq, nvar = nvar, 
                 method = method, est =  est, tkernel = tkernel, bw.cov = bw.cov, 
                 level = 0, runs = 0, tboot = NULL, BOOT = NULL,
                 R = R, r = r, control = control, formula = formula, call = match.call())
  class(result) <- "tvsure"
  return(result)
}


