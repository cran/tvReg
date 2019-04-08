
# tvReg

This R package covers a large range of semiparametric regression methods with time-varying coefficients.


## Installation

You can install the released version of tvReg from [CRAN](https://CRAN.R-project.org) with:

```
install.packages("tvReg")
```

or the development version from GitHub with:
```
devtools::install_github("icasas/tvReg")
```

## Main functions

The five basic functions in this package are `tvLM()`, `tvAR()`, `tvSURE()`, `tvVAR()` and `tvIRF()`. Moreover, this package provides the `confint()`, `fitted()`, `forecast()`, `plot()`, `predict()`, `print()`, `resid()` and `summary()` methods adapted to the class attributes of the `tvReg`. In addition, it includes bandwidth selection methods, time-varying variance-covariance estimators and two estimation procedures: the time-varying ordinary least squares, which are implemented in the `tvOLS()` methods and the time-varying generalised least squares, which are implemented in the `tvGLS()` methods.

## Further information

Details on the theory and applications to finance and macroeconomics can be found in [Casas, Isabel and Fernandez-Casal, 2019](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3363526),
and in the package vignette <https://icasas.github.io/tvReg/articles/tvReg.html>.


## References

Casas, Isabel and Fernandez-Casal, Ruben, *tvReg: Time-varying Coefficient Linear Regression for Single and Multi-Equations in R* (April 1, 2019). Available at SSRN:<https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3363526>.



