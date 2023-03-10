Making sense of censored covariates: Statistical methods for studies of
Huntington’s disease
================

## Software Supplement

## Lotspeich et al. 

``` r
# Set seed for reproducibility
set.seed(128)
```

### DATA GENERATION

Simulate a dataset for a sample of $n = 1000$ observations under heavy
(\~40%) **random right censoring**. Specifically,

- $Z \sim \textrm{Bernoulli}(0.5)$,
- $X \sim \textrm{Weibull}(0.75, 0.25)$,
- $Y = 1 + 0.5X + 0.25Z + \epsilon$ (for
  $\epsilon \sim \textrm{N}(0, 1)$),
- $C \sim \textrm{Expo}(2.9)$,
- $W = \min(X, C)$, and
- $\Delta = \textrm{I}(X \leq C)$.

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

### NAIVE ANALYSIS (Section 4)

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

### COMPLETE CASE ANALYSIS (Section 5)

``` r
# Complete case analysis
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

### WEIGHTED METHODS (Section 6)

``` r
# Load the survival package 
library(survival) ## used to fit models for the weights 
```

#### Inverse probability weighting (IPW)

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

#### Augmented inverse probability weighting (AIPW)

### IMPUTATION METHODS (Section 7)

#### Single imputation

#### Multiple imputation

### MAXIMUM LIKELIHOOD ESTIMATION (MLE) (Section 8)

#### Parametric MLE

#### Semiparametric MLE

### BAYESIAN METHODS (Section 9)

### OTHER METHODS (Section 10)

#### Thresholding (dichotomization)

#### Subtitution

#### Reverse survival regression
