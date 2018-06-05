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
#' \item{tvcoef}{An array of dimension obs x neq (obs = number of observations,
#' neq = number of equations in the system) with the time-varying coefficients estimates.}
#' \item{fitted}{The fitted values.}
#' \item{residuals}{Estimation residuals.}
#' \item{x}{A list with the regressors data and the dependent variable.}
#' \item{y}{A matrix with the dependent variable data.}
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
#' model.VAR <- vars::VAR(usmacro, p = 6, type = "const")
#' model.tvVAR <- tvVAR(usmacro, p = 6, type = "const", bw = c(1.8, 20, 20))
#' plot(model.tvVAR)
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

tvVAR <- function (y, p = 1, z = NULL, bw = NULL, type = c("const", "none"), exogen = NULL,
                   est = c("lc", "ll"), tkernel = c("Epa", "Gaussian"), singular.ok = TRUE)
{
  y <- as.matrix(y)
  if (any(is.na(y)))
    stop("\nNAs in y.\n")
  if (ncol(y) < 2)
    stop("\nMatrix 'y' should contain at least two variables. For univariate
          analysis consider the 'tvAR' function.\n")
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  type <- match.arg(type)
  if(!(tkernel %in% c("Epa","Gaussian")))
    tkernel <- "Epa"
  if(!(est %in% c("lc", "ll")))
    est <- "lc"
  if(!(type %in% c("const", "none")))
    type <- "const"
  if (is.null(colnames(y)))
  {
    colnames(y) <- paste("y", 1:ncol(y), sep = "")
    warning(paste("No column names supplied in y, using:",
                  paste(colnames(y), collapse = ", "), ", instead.\n"))
  }
  colnames(y) <- make.names(colnames(y))
  y.orig <- y
  type <- match.arg(type)
  obs <- dim(y)[1]
  neq <- dim(y)[2]
  sample <- obs - p
  ylags <- stats::embed(y, dimension = p + 1)[, -(1:neq)]
  temp1 <- NULL
  for (i in 1:p)
  {
    temp <- paste(colnames(y), ".l", i, sep = "")
    temp1 <- c(temp1, temp)
  }
  colnames(ylags) <- temp1
  yend <- y[-c(1:p), ]
  if (type == "const")
  {
    rhs <- cbind( ylags, rep(1, sample))
    colnames(rhs) <- c(colnames(ylags), "Intercept")
  }
  else if (type == "none")
  {
    rhs <- ylags
    colnames(rhs) <- colnames(ylags)
  }
  if (!(is.null(exogen)))
  {
    exogen <- as.matrix(exogen)
    if (!identical(nrow(exogen), nrow(y))) {
      stop("\nDifferent row size of 'y' and exogen.\n")
    }
    if (is.null(colnames(exogen))) {
      colnames(exogen) <- paste("exo", 1:ncol(exogen),
                                sep = "")
    }
    colnames(exogen) <- make.names(colnames(exogen))
    tmp <- colnames(rhs)
    rhs <- cbind(rhs, exogen[-c(1:p), ])
    colnames(rhs) <- c(tmp, colnames(exogen))
  }
  datamat <- as.data.frame(rhs)
  colnames(datamat) <- colnames(rhs)
  equation <- list()
  if(is.null(bw))
   bw <- bw(y = yend, x = datamat, tkernel = tkernel, est = est, singular = singular.ok)
  resid = fitted <- matrix(0, nrow = sample, ncol = neq)
  for (i in 1:neq)
  {
      y <- yend[, i]
      results <- tvOLS(x = datamat, y = y, bw = bw[i], est = est, tkernel = tkernel,
                       singular.ok = singular.ok)
      equation[[colnames(yend)[i]]] <- results$tvcoef
      colnames(equation[[colnames(yend)[i]]]) <- colnames(rhs)
      resid[,i] <- results$residuals
      fitted[, i] <- results$fitted
  }
  colnames(resid) <- names(equation)
  colnames(fitted) <- colnames(resid)
  if(length(bw) == 1)
    names(bw) <- "bw.mean"
  else
    names(bw) <- paste("bw.", names(equation), sep = "")
  result <- list(tvcoef = equation, Lower = NULL, Upper = NULL,  fitted = fitted,
                 residuals = resid, datamat = data.frame(cbind(yend,rhs)), y = y.orig,
                 exogen = exogen, p = p, type = type, obs = sample, totobs = sample + p,
                 neq = neq, est = est, tkernel = tkernel, bw = bw, call = match.call())
  class(result) <- "tvvar"
  return(result)
}