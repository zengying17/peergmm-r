# PeerCRA

`PeerCRA` estimates endogenous peer effects under conditional random assignment
to peer groups within urns, or selection pools. It implements a two-step GMM
estimator with optional own-effect and peer-effect covariates and a
heteroskedasticity-robust, bias-corrected covariance estimator.

The package is written in base R and contains no compiled code.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("zengying17/peergmm-r")
```


## Data structure

The input must be a data frame with one row per individual and columns for:

- the outcome;
- the urn or selection-pool identifier;
- the peer-group identifier, nested within the urn;
- any covariates included as own effects or peer effects.

For a variable named in `peer_vars`, `PeerCRA` constructs its leave-one-out
peer-group average. Missing estimator variables are removed before the peer
structure is checked. The fitted object records the retained observations and
any covariates dropped because they are absorbed or numerically dependent.

## Example

```r
library(PeerCRA)

fit <- peer_cra(
  data      = mydata,
  y_var     = "Y",
  urn_var   = "urn_id",
  group_var = "group_id",
  own_vars  = c("x1", "x2"),
  peer_vars = c("x1")
)

print(fit)
summary(fit)
coef(fit)
vcov(fit)
confint(fit)
head(residuals(fit))
head(fitted(fit))
```

The coefficient named `lambda` is the endogenous peer-effect parameter.
Coefficients prefixed by `ave_` correspond to leave-one-out peer averages.

## Returned object

`peer_cra()` returns an object of class `peer_cra`. Its main components include:

- point estimates and their covariance matrix;
- standard errors and Wald inference when covariance inference is available;
- estimated urn fixed effects;
- fitted values and structural residuals;
- the estimation-sample indicator;
- stage-one and stage-two criteria and optimizer diagnostics;
- covariance and numerical-rank diagnostics.

Run `?peer_cra` after installation for the complete argument and return-value
documentation.

## Repository contents

- `DESCRIPTION`: package metadata and dependencies.
- `NAMESPACE`: exported function and S3 method registrations.
- `R/`: package source.
- `man/`: generated R help files.

`NAMESPACE` and the files in `man/` are generated from the roxygen comments in
`R/`, but they are committed so the package can be installed directly from
GitHub without regenerating documentation.


## Companion Stata command

The companion Stata implementation is available at
[zengying17/peergmm-stata](https://github.com/zengying17/peergmm-stata).

## Citation

Ying Zeng, *Estimation and Inference for Peer Effects under Conditional Random
Assignment*. Citation details will be updated when the paper is released.

## License

GPL (>= 3).
