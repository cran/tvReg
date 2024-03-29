#' Time-varying Vector Autoregressive Models
#'
#' Fits a time-varying coefficients vector autorregressive model with p lags.
#' @aliases tvvar-class tvvar
#' @rdname tvVAR
#' @import bvarsv
#' @param y A matrix with dimention obs x neq (obs = number of observations and
#' neq = number of equations)
#' @param p A scalar indicating the number of lags in the model
#' @param type A character 'const' if the model contains an intercept and 'none' otherwise.
#' @param exogen A matrix or data.frame with the exogenous variables (optional)
#' @inheritParams tvSURE
#' @param singular.ok	Logical. If FALSE, a singular model is an error.
#' @return An object of class 'tvvar'
#' The object of class \code{tvvar} have the following components:
#' \item{coefficients}{An array of dimension obs x neq (obs = number of observations,
#' neq = number of equations in the system) with the time-varying coefficients estimates.}
#' \item{fitted}{The fitted values.}
#' \item{residuals}{Estimation residuals.}
#' \item{x}{A list with the regressors data and the dependent variable.}
#' \item{y}{A matrix with the dependent variable data.}
#' \item{z}{A vector with the smoothing variable.}
#' \item{ez}{A vector with the smoothing estimation values.}
#' \item{bw}{Bandwidth of mean estimation.}
#' \item{type}{Whether the model has a constant or not.}
#' \item{exogen}{A matrix or data.frame with other exogenous variables.}
#' \item{p}{Number of lags}
#' \item{neq}{Number of equations}
#' \item{obs}{Number of observations in estimation.}
#' \item{totobs}{Number of observations in the original set.}
#' \item{call}{Matched call.}
#' 
#' @seealso \code{\link{bw}}, \code{\link{tvIRF}}, \code{\link{plot}}, 
#' \code{\link{print}} and \code{\link{summary}}
#' 
#' @examples
#' ##Inflation rate, unemployment rate and treasury bill interest rate for 
#' ##the US, as used in Primiceri (2005).
#' data(usmacro, package = "bvarsv")
#' VAR.fit <- vars::VAR(usmacro, p = 6, type = "const")
#' tvVAR.fit <- tvVAR(usmacro, p = 6, type = "const", bw = c(1.8, 20, 20))
#' plot(tvVAR.fit)
#' 
#' @references 
#' Casas, I., Ferreira, E., and Orbe, S. (2017) Time-Varying Coefficient Estimation 
#' in SURE Models: Application to Portfolio Management. Available at SSRN: 
#' https://ssrn.com/abstract=3043137
#' 
#' Primiceri, G.E. (2005) Time varying structural vector autoregressions 
#' and monetary policy. \emph{Review of Economic Studies}, 72, 821-852.
#' 
#' @export

tvVAR <- function (y, p = 1, z = NULL, ez = NULL, bw = NULL, cv.block = 0, type = c("const", "none"), 
                   exogen = NULL, est = c("lc", "ll"), tkernel = c("Triweight", "Epa", "Gaussian"), 
                   singular.ok = TRUE)
{
  y <- as.matrix(y)
  if (any(is.na(y)))
    stop("\nNAs in y.\n")
  toobs <- dim(y)[1]
  neq <- dim(y)[2]
  obs <- toobs - p
  if(!is.null(z))
  {
   if(!inherits(z, c("numeric", "vector")))
      stop("\nArgument 'z' should be 'numeric' or a 'vector'.\n")
    z <- as.numeric(z)[-c(1:p)]
  }
  if (NCOL(y) < 2)
    stop("\nMatrix 'y' should contain at least two variables. For univariate
          analysis consider the 'tvAR' function.\n")
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  type <- match.arg(type)
  if (is.null(colnames(y)))
  {
    colnames(y) <- paste0("y", 1:NCOL(y))
    warning(paste("No column names supplied in y, using:",
                  paste(colnames(y), collapse = ", "), ", instead.\n"))
  }
  var.names <- make.names(colnames(y))
  colnames(y) <- var.names
  y.orig <- y
  type <- match.arg(type)
  
  rhs <- stats::embed(y, dimension = p + 1)[, -(1:neq)]
  colnames(rhs) <- paste0(colnames(y), paste0(".l", rep(1:p, each=length(var.names))))
  yend <- y[-c(1:p), ]
  if (type == "const")
  {
    tmp <- colnames(rhs)
    rhs <- cbind(rhs, 1L)
    colnames(rhs) <- c(tmp, "(Intercept)")
  }
  if (!(is.null(exogen))) 
  {
    if(NCOL(exogen)==1)
      exogen <- matrix(exogen)
    exogen <- as.matrix(exogen[-c(1:p),])
    names.exo <- colnames(exogen)
    if (!identical(NROW(exogen), NROW(yend)))
      stop("\nDifferent row size of 'y' and exogen.\n")
    if (is.null(names.exo)) 
      names.exo <- c(paste0("exo", 1:NCOL(exogen)))
    names.exo <- c(colnames(rhs), names.exo)
    rhs <- cbind(rhs, exogen)
    colnames(rhs) <- names.exo
  }
  equation <- list()
  if(is.null(bw))
  {
    cat("Calculating regression bandwidths... ")
    bw <- bw(y = yend, x = rhs, z = z, cv.block = cv.block, tkernel = tkernel, 
             est = est, singular = singular.ok)
    cat("bandwidth(s) ", bw, "\n")
  }
  resid = fitted <- matrix(0, nrow = obs, ncol = neq)
  for (i in 1:neq)
  {
      y <- yend[, i]
      results <- tvOLS(x = rhs, y = y, z = z, ez = ez, bw = bw[i], est = est, tkernel = tkernel,
                       singular.ok = singular.ok)
      equation[[colnames(yend)[i]]] <- results$coefficients
      colnames(equation[[colnames(yend)[i]]]) <- colnames(rhs)
      resid[,i] <- results$residuals
      fitted[, i] <- results$fitted
  }
  colnames(resid) <- var.names
  colnames(fitted) <- var.names
  colnames(yend) <- var.names
  if(length(bw) == 1)
    names(bw) <- "bw.mean"
  else
    names(bw) <- paste0("bw.", names(equation))
  result <- list(coefficients = equation, Lower = NULL, Upper = NULL,  fitted = fitted,
                 residuals = resid, y = yend, x = rhs, z = z, ez = ez, y.orig = y.orig,
                 bw = bw, cv.block = cv.block, exogen = exogen, p = p, type = type, obs = obs, 
                 totobs = toobs, neq = neq, est = est, tkernel = tkernel, 
                 singular.ok = singular.ok, call = match.call())
  class(result) <- "tvvar"
  return(result)
}
