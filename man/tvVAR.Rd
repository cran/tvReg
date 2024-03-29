% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tvVAR.R
\name{tvVAR}
\alias{tvVAR}
\alias{tvvar-class}
\alias{tvvar}
\title{Time-varying Vector Autoregressive Models}
\usage{
tvVAR(
  y,
  p = 1,
  z = NULL,
  ez = NULL,
  bw = NULL,
  cv.block = 0,
  type = c("const", "none"),
  exogen = NULL,
  est = c("lc", "ll"),
  tkernel = c("Triweight", "Epa", "Gaussian"),
  singular.ok = TRUE
)
}
\arguments{
\item{y}{A matrix with dimention obs x neq (obs = number of observations and
neq = number of equations)}

\item{p}{A scalar indicating the number of lags in the model}

\item{z}{A vector containing the smoothing variable.}

\item{ez}{(optional) A scalar or vector with the smoothing values. If 
values are not included then the vector \code{z} is used instead.}

\item{bw}{An optional scalar or vector of length the number of equations. It represents the bandwidth in
the estimation of trend coefficients. If NULL, it is selected by cross validation.}

\item{cv.block}{A positive scalar with the size of the block in leave one block out cross-validation.
By default 'cv.block = 0' meaning leave one out cross-validation.}

\item{type}{A character 'const' if the model contains an intercept and 'none' otherwise.}

\item{exogen}{A matrix or data.frame with the exogenous variables (optional)}

\item{est}{The nonparametric estimation method, one of "lc" (default) for linear constant or "ll" for local linear.}

\item{tkernel}{A character, either "Triweight" (default), "Epa" or "Gaussian" kernel function.}

\item{singular.ok}{Logical. If FALSE, a singular model is an error.}
}
\value{
An object of class 'tvvar'
The object of class \code{tvvar} have the following components:
\item{coefficients}{An array of dimension obs x neq (obs = number of observations,
neq = number of equations in the system) with the time-varying coefficients estimates.}
\item{fitted}{The fitted values.}
\item{residuals}{Estimation residuals.}
\item{x}{A list with the regressors data and the dependent variable.}
\item{y}{A matrix with the dependent variable data.}
\item{z}{A vector with the smoothing variable.}
\item{ez}{A vector with the smoothing estimation values.}
\item{bw}{Bandwidth of mean estimation.}
\item{type}{Whether the model has a constant or not.}
\item{exogen}{A matrix or data.frame with other exogenous variables.}
\item{p}{Number of lags}
\item{neq}{Number of equations}
\item{obs}{Number of observations in estimation.}
\item{totobs}{Number of observations in the original set.}
\item{call}{Matched call.}
}
\description{
Fits a time-varying coefficients vector autorregressive model with p lags.
}
\examples{
##Inflation rate, unemployment rate and treasury bill interest rate for 
##the US, as used in Primiceri (2005).
data(usmacro, package = "bvarsv")
VAR.fit <- vars::VAR(usmacro, p = 6, type = "const")
tvVAR.fit <- tvVAR(usmacro, p = 6, type = "const", bw = c(1.8, 20, 20))
plot(tvVAR.fit)

}
\references{
Casas, I., Ferreira, E., and Orbe, S. (2017) Time-Varying Coefficient Estimation 
in SURE Models: Application to Portfolio Management. Available at SSRN: 
https://ssrn.com/abstract=3043137

Primiceri, G.E. (2005) Time varying structural vector autoregressions 
and monetary policy. \emph{Review of Economic Studies}, 72, 821-852.
}
\seealso{
\code{\link{bw}}, \code{\link{tvIRF}}, \code{\link{plot}}, 
\code{\link{print}} and \code{\link{summary}}
}
