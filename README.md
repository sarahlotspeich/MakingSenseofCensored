Making sense of censored covariates: Statistical methods for studies of
Huntington’s disease
================
### Lotspeich, Ashner, Vazquez, Richardson, Grosser, Bodek and Garcia (2023+) 
## Software Supplement

### LOAD PACKAGES

``` r
## Load the censCov package (Qian et al 2018)
### contains thlm function for thresholding methods, complete case, and reverse survival regression
library(censCov)

# Load the survival package 
### used to fit models for the weights for IPW
library(survival) 
```

### DATA GENERATION

``` r
# Set seed for reproducibility
set.seed(128)
```

Simulate a dataset for a sample of *n* = 1000 observations under heavy
(\~40%) **random right censoring**. Specifically,

-   *Z* ∼ Bernoulli(0.5),
-   *X* ∼ Weibull(0.75, 0.25),
-   *Y* = 1 + 0.5*X* + 0.25*Z* + *ϵ* (for *ϵ* ∼ N(0, 1)),
-   *C* ∼ Expo(2.9),
-   *W* = min (*X*, *C*), and
-   *Δ* = I(*X* ≤ *C*).

``` r
# Create simulated dataset 
z = rbinom(n = 1000, size = 1, prob = 0.5) ## Uncensored covariate
x = rweibull(n = 1000, shape = 0.75, scale = 0.25)  ## To-be-censored covariate
e = rnorm(n = 1000, mean = 0, sd = 1) ## Random errors
y = 1 + 0.5 * x + 0.25 * z + e ## Continuous outcome
c = rexp(n = 1000, rate = 2.9) ## Random censoring mechanism
w = pmin(x, c) ## Observed covariate value
d = as.numeric(x <= c) ## "Event" indicator
random_right_dat = data.frame(z, w, y, d) ## Construct data set
```

### SECTION 3: REVERSE SURVIVAL REGRESSION

``` r
# Reverse survival regression "by hand"
## Fit survival model where x is treated as the outcome rather than a covariate
fit_reverse_survival = coxph(formula = Surv(time = w, event = d) ~ y + z, data = random_right_dat)
## Report results for a test of association between x and y, while controlling for z
## Note: parameter interpretations change when we use survival regression, so we
##       only report the hypothesis test result here
summary(fit_reverse_survival)$coefficients[1,5]
```

    ## [1] 1.795412e-08

``` r
# Reverse survival regression with thlm
## Fit survival model where x is treated as the outcome rather than a covariate
fit_reverse_survival_thlm = thlm(y ~ w + z, cens = d, data = random_right_dat,
                            method = "reverse", control = list(t0.plot = FALSE))
## Report results for a test of association between x and y, while controlling for z
## Note: parameter interpretations change when we use survival regression, so we
##       only report the hypothesis test result here
fit_reverse_survival_thlm
```

    ## 
    ##  Call: thlm(formula = y ~ w + z, data = random_right_dat, cens = d, 
    ##     method = "reverse", control = list(t0.plot = FALSE))
    ## 
    ##  Hypothesis test of association, H0: a1 = 0
    ## p-value = 0.0000

### SECTION 4: AD HOC METHODS THAT PRODUCE BIAS

#### 4.1 Naive Analysis

``` r
# Naive analysis
## Fit the model using observed w in place of censored x
fit_naive = lm(y ~ w + z, data = random_right_dat)
coeff_naive <- fit_naive$coefficients
se_naive <- sqrt(diag(vcov(fit_naive)))
## Inspect results
data.frame(coeff = coeff_naive, se = se_naive)
```

    ##                 coeff         se
    ## (Intercept) 0.9640558 0.05657472
    ## w           0.9586777 0.19591964
    ## z           0.2408399 0.06669586

#### 4.2 Substitution

### SECTION 5: COMPLETE CASE ANALYSIS

``` r
# Complete case analysis "by hand"
## Remove observations with censored x
random_right_dat_complete = random_right_dat[random_right_dat$d == 1, ]
## Fit the model to "complete" cases (with uncensored x)
fit_complete = lm(y ~ w + z, data = random_right_dat_complete)
coeff_complete <- fit_complete$coefficients
se_complete <- sqrt(diag(vcov(fit_complete)))
## Inspect results
data.frame(coeff = coeff_complete, se = se_complete)
```

    ##                 coeff         se
    ## (Intercept) 0.9530811 0.07529996
    ## w           0.8380126 0.30051911
    ## z           0.1377386 0.08724157

``` r
# Complete case analysis with thlm
fit_complete_thlm = thlm(y ~ w + z, cens = d, data = random_right_dat,
                            method = "cc", control = list(t0.plot = FALSE))
## estimating std dev of coefficient for x
coeff_complete_thlm <- c(fit_complete_thlm$a1, fit_complete_thlm$a2)
se_complete_thlm <- c(fit_complete_thlm$a1.sd, fit_complete_thlm$a2.sd)
## Inspect results
data.frame(coeff = coeff_complete_thlm, se = se_complete_thlm, row.names = c("x", "z"))
```

    ##       coeff         se
    ## x 0.8380126 0.30051911
    ## z 0.1377386 0.08724157

### SECTION 6: WEIGHTING METHODS

#### 6.1 Inverse probability weighting (IPW)

``` r
# IPW analysis
## Estimate probabilities of being uncensored
### Fit a model for the hazard function of being censored given z 
fit_weights = coxph(formula = Surv(time = w, event = (1-d)) ~ z, 
                    data = random_right_dat)
### Use the model to estimate the survival function of c
### (i.e., probability of being uncensored) for each complete observation
prob_uncensor = survfit(formula = fit_weights, 
                        newdata = random_right_dat_complete)
### Use the summary() function to get p(d = 1|z) for each complete case



weights = data.frame(x = summary(object = prob_uncensor, 
                                 times = extend)$time, 
                     ipw = 1 / diag(summary(object = prob_uncensor, extend = TRUE)$surv))

weights = data.frame()
weights_df <- data.frame(X = summary(G, times = data_sim$X, extend = TRUE)$time,
                        weights = summary(G, times = data_sim$X, extend = TRUE)$surv %>% diag())
      data_sim <- data_sim %>% left_join(weights_df, by = "X")

## Remove observations with censored x
random_right_dat_complete = random_right_dat[random_right_dat$d == 1, ]
## Fit the model to "complete" cases (with uncensored x)
fit_complete = lm(y ~ w + z, data = random_right_dat_complete)
coeff_complete <- fit_complete$coefficients
se_complete <- sqrt(diag(vcov(fit_complete)))
## Inspect results
data.frame(coeff = coeff_complete, se = se_complete)
```

#### 6.2 Augmented inverse probability weighting (AIPW)

There is no simple implementation of AIPW for censored covariates that
we are aware of. However, Ahn et al 2018 share their code with their
manuscript (this is for a model with an interval censored covariate and
a time-to-event outcome). The code can be found in the Supporting
Information section [here](https://doi.org/10.1002/bimj.201700090).

### SECTION 7: IMPUTATION METHODS

#### Single imputation

#### Multiple imputation

### SECTION 8: MAXIMUM LIKELIHOOD METHODS

#### Parametric MLE

```{r}
# Parametric MLE analysis
## Define the log-likelihood function based on Equation (3) 
## Per the data generation, assume a normal for Y|X,Z and a Weibull for X|Z
joint_dens = function(params) {

}

loglik = function(params) {
    # Separate the parameter vector into the two models 
    ## Model parameters for Y|X,Z
    theta = params[1:4]
    
    ## Model parameters for X|Z
    eta = params[-c(1:4)] 
    
    # Compute log-likelihood contributions for uncensored data 
    ll = d * log(dnorm(y, mean = (theta[1] + theta[2] * w + theta[3] * z), sd = theta[4]) * dweib(w, shape = eta[1], scale = (eta[2] + eta[3] * z)))
    
    # Add log-likelihood contributions for censored data
    ## Write internal function to integrate censored x out of the joint for each observation
    integrate_joint = function(data_row) {
      Wi <- as.numeric(data_row[W])
      Yi <- data_row[Y]
      Zi <- data_row[Z]
      return(
        tryCatch(expr = integrate(f = function(x) dnorm(y, mean = (theta[1] + theta[2] * x + theta[3] * z), sd = theta[4]) * dweib(x, shape = eta[1], scale = (eta[2] + eta[3] * z)),
                                  lower = Wi, 
                                  upper = Inf)$value,
                 error = function(err) {0})
      )
    }
    
    # Create simulated dataset 
    z = rbinom(n = 1000, size = 1, prob = 0.5) ## Uncensored covariate
    x = rweibull(n = 1000, shape = 0.75, scale = 0.25)  ## To-be-censored covariate
    e = rnorm(n = 1000, mean = 0, sd = 1) ## Random errors
    y = 1 + 0.5 * x + 0.25 * z + e ## Continuous outcome
    c = rexp(n = 1000, rate = 2.9) ## Random censoring mechanism
    w = pmin(x, c) ## Observed covariate value
    d = as.numeric(x <= c) ## "Event" indicator
    random_right_dat = data.frame(z, w, y, d) ## Construct data set
}
```

```{r}
loglik <- function(params, Y, X, W, D, Z = NULL, data, subdivisions = 100, distY = "normal", distX = "normal", cens = "right") {
  ####################################################
  # Pre-processing ###################################
  ####################################################
  # < number of uncensored subjects > ----------------
  n1 <- sum(data[, D]) # -----------------------------
  # ---------------- < number of uncensored subjects >
  # Reordered data to be uncensored first ------------
  data <- data[order(data[, D], decreasing = TRUE), ]
  # ------------ Reordered data to be uncensored first
  # Create subset of uncensored subjects' data -------
  uncens_data <- data[1:n1, ]
  # ------- Create subset of uncensored subjects' data
  # Create subset of censored subjects' data -------
  cens_data <- data[-c(1:n1), ]
  # ------- Create subset of censored subjects' data

  ####################################################
  # Joint density P(Y,X,Z) ###########################
  ####################################################
  pYXandZ_uncens <- calc_pYXandZ(x = uncens_data[, X],
                                 y = uncens_data[, Y],
                                 z = uncens_data[, Z],
                                 lengthZ = length(Z),
                                 distY = distY,
                                 distX = distX,
                                 params = params)

  ####################################################
  # Likelihood (Uncensored) ##########################
  ####################################################
  if (any(is.na(pYXandZ_uncens))) {
    # If params are out of domain, calc_pYXandZ returns NA
    ## And the log-likelihood needs to be arbitrarily "huge"
    return(1E8)
  } else {
    # Replace P(Y,X,Z) = 0 with P(Y,X,Z) = 1 so that
    ## log P(Y,X,Z) = 0.
    pYXandZ_uncens[pYXandZ_uncens == 0] = 1
    ll <- sum(log(pYXandZ_uncens))
  }

  ####################################################
  # Likelihood (Censored) ############################
  ####################################################
  if (nrow(cens_data) > 0) {
    integrate_pYXandZ <- function(data_row) {
      Wi <- as.numeric(data_row[W])
      Yi <- data_row[Y]
      Zi <- data_row[Z]
      return(
        tryCatch(expr = integrate(f = calc_pYXandZ,
                                  lower = ifelse(test = cens == "right", Wi, -Inf),
                                  upper = ifelse(test = cens == "right", Inf, Wi),
                                  subdivisions = subdivisions,
                                  y = Yi,
                                  z = Zi,
                                  lengthZ = length(Z),
                                  distY = distY,
                                  distX = distX,
                                  params = params)$value,
                 error = function(err) {0})
      )
    }
    int_pYXandZ_cens <- apply(X = cens_data,
                              MARGIN = 1,
                              FUN = integrate_pYXandZ)
    log_int_pYXandZ_cens <- log(int_pYXandZ_cens)
    log_int_pYXandZ_cens[log_int_pYXandZ_cens == -Inf] <- 0
    ll <- ll + sum(log_int_pYXandZ_cens)
  }

  # Return (-1) x log-likelihood for use with nlm() --
  return(- ll)
  # -- Return (-1) x log-likelihood for use with nlm()
}
```

#### Semiparametric MLE

There is no simple implementation of semiparametric MLE for censored covariates that
we are aware of. However, Kong and Nan give an overview of the necessary steps to coding on in Section 3 of their paper, which can be found [here](https://academic.oup.com/biomet/article/103/1/161/2389883#113398541).

### SECTION 9: BAYESIAN METHODS

There is no simple implementation of Bayesian Methods for censored
covariates that we are aware of. However, Wu et al 2012 share their
WinBUGS code in the Appendix of their manuscript, which can be found
[here](https://doi.org/10.1080/02664763.2012.681362).

### SECTION 10: THRESHOLD METHODS (DICHOTOMIZATION)

#### Deletion Thresholding

``` r
# Deletion Thresholding analysis
## Fit model using dichotomized covariate and with indeterminate observations deleted
fit_deletion_threshold = thlm(y ~ w + z, cens = d, data = random_right_dat,
                              method = "dt", control = list(t0.plot = FALSE),
                              B = 100) # 100 bootstrap replicates for estimating std dev of coefficient for x
coeff_deletion_threshold <- c(fit_deletion_threshold$a1, fit_deletion_threshold$a2)
se_deletion_threshold <- c(fit_deletion_threshold$a1.sd, fit_deletion_threshold$a2.sd)
## Inspect results
data.frame(coeff = coeff_deletion_threshold, se = se_deletion_threshold, row.names = c("x", "z"))
```

    ##       coeff         se
    ## x 0.7118696 0.16830116
    ## z 0.1944697 0.07969871

#### Complete Thresholding

``` r
# Complete Thresholding analysis
## Fit model using dichotomized covariate derived from w
fit_complete_threshold = thlm(y ~ w + z, cens = d, data = random_right_dat,
                              method = "ct", control = list(t0.plot = FALSE),
                              B = 100) # 100 bootstrap replicates for estimating std dev of coefficient for x
coeff_complete_threshold <- c(fit_complete_threshold$a1, fit_complete_threshold$a2)
se_complete_threshold <- c(fit_complete_threshold$a1.sd, fit_complete_threshold$a2.sd)
## Inspect results
data.frame(coeff = coeff_complete_threshold, se = se_complete_threshold, row.names = c("x", "z"))
```

    ##       coeff        se
    ## x 0.7682436 0.2359501
    ## z 0.2391911 0.0670703
