Making sense of censored covariates: Statistical methods for studies of
Huntington’s disease
================

## Software Supplement

## Lotspeich et al. 

``` r
# Set seed for reproducibility
set.seed(128)
```

``` r
# packages needed 

## Load the censCov package (Qian et al 2018)
### contains thlm function for thresholding methods, complete case, and reverse survival regression
library(censCov)

# Load the survival package 
### used to fit models for the weights for IPW
library(survival) 
```

### DATA GENERATION

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

### SECTION 3: USING REVERSE SURVIVAL REGRESSION TO HANDLE A CENSORED COVARIATE

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

### SECTION 4: METHODS THAT PRODUCE BIAS

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

``` r
# Substitution
## Define function which will be used to replace censored W
f.w <- function(w){
  return(w / 2) # for illustrative purposes, we replace w with w / 2
}
## Create copy of dataset to modify
random_right_dat_copy <- random_right_dat
## replace censored w with w / 2
random_right_dat_copy$fw <- ifelse(d==1, f.w(w),w)
## fit linear model on modified dataset
fit_substitution = lm(y ~ fw + z, data = random_right_dat_copy)
coeff_substitution <- fit_substitution$coefficients
se_substitution <- sqrt(diag(vcov(fit_substitution)))
## Inspect results
data.frame(coeff = coeff_substitution, se = se_substitution)
```

    ##                 coeff         se
    ## (Intercept) 0.9779523 0.05406085
    ## fw          1.1668621 0.22103866
    ## z           0.2369362 0.06657209

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
## Estimate probabilities of being uncensored. Can use a logisitc regression, 
## Cox-proportional hazards model, Kapla-Meier curve

### In this example we fit a logisitc regression model using the logit link 
### function given z on the probability of being observed. 
mylogisitc = glm(d ~ z, data=random_right_dat, family = "binomial")

# Estimate weights for only observed cases: 1/Pr(Delta=1)
random_right_dat_complete$weights_logisitc = 
  1/predict(mylogisitc, random_right_dat_complete)

## Fit the model to "complete" weighted cases (with uncensored x)
fit_ipw = lm(y ~ w + z, data = random_right_dat_complete, weights = weights_logisitc)
coeff_ipw <- fit_ipw$coefficients
se_ipw <- sqrt(diag(vcov(fit_ipw)))
## Inspect results
data.frame(coeff = coeff_ipw, se = se_ipw)
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

### SECTION 8: MAXIMUM LIKELIHOOD ESTIMATION (MLE)

#### Parametric MLE

#### Semiparametric MLE

There is no simple implementation of semiparametric MLE for censored covariates that we are aware of. However, Kong and Nan give an overview of the necessary steps to coding on in Section 3 of their paper, which can be found [here](https://academic.oup.com/biomet/article/103/1/161/2389883#113398541).

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
