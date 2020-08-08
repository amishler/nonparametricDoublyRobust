packages <- c("boot", "dplyr", "mvtnorm", "readr", "sandwich", "SuperLearner",
              "tibble", "tidyr", "tmle")
for (package in packages) {
  library(package, character.only=T)
}


#' Expit function, aka inverse logit.
#'
#' @param x Real number
expit <- function(x){ exp(x)/(1+exp(x)) }

#' Generate indices for cross validates folds.
#'
#' @param n Number of rows.
#' @param cv_folds Number of folds.
#' @return Vector of indices indicating fold membership.
generate_folds <- function(n, cv_folds) {
  inds <- sample(rep(1:cv_folds, ceiling(n/cv_folds))[1:n])
  return(inds)
}

#' Generate covariates X and transform it to induce model misspecification.
#'
#' Data is 4-dimensional multivariate normal with mean 0 and diagonal
#' covariance.
#' @param size Number of samples.
#' @return A list containing two n-by-4 matrices: X, the baseline covariates, and
#' Z, the transformed version of X.
generate_covariates <- function(size) {
  X <- rmvnorm(size, mean = rep(0, 4), sigma = diag(4))
  colnames(X) <- paste0("X", 1:4)
  Z <- X
  Z[,1] <- exp(X[,1]/2)
  Z[,2] <- X[,2]/(1 + exp(X[,1])) + 10
  Z[,3] <- (X[,1]*X[,3]/25 + 0.6)^3
  Z[,4] <- (X[,2] + X[,4] + 20)^2
  colnames(Z) <- paste0("Z", 1:4)

  return(list(X = X, Z = Z))
}

#' Generate treatments and outcomes from covariates.
#'
#' Propensity scores pi are a logistic function of the covariates X. The
#' treatment is binary, 0 or 1. The outcomes are linear in the covariates, plus
#' independent Gaussian noise, plus a constant treatment effect for those who
#' received treatment. Given this structure, the treatment and the potential
#' outcomes are independent conditional on X, but not marginally because X
#' confounds the treatment and the outcome.
#' @param X n-by-p matrix of covariate data.
#' @param beta_a Vector of logistic coefficients of length p, relating the
#' propensity scores to X.
#' @param beta_y Vector of linear coefficients of length p, relating the outcome
#' Y to X.
#' @param effect_size The treatment effect, which is constant (not
#' heterogeneous).
#' @return A list containing treatments A and outcomes Y.
generate_outcomes <- function(X, beta_a, beta_y, effect_size=6) {
  n <- nrow(X)
  X <- cbind(1, X)                  # add an intercept column
  pi <- expit(X %*% beta_a)         # propensity scores
  A <- 1 - rbinom(n, 1, pi)         # treatment levels
  mu0 <- X %*% beta_y               # E[Y(A = 0)|X]
  Y <- A*effect_size + mu0 + rnorm(n, 0, 6)   # Y

  return(list(A = A, Y = c(Y)))
}


#' Estimate propensity scores parametrically, via logistic regression.
#'
#' @param A A vector of treatment values, in {0, 1}.
#' @param X A matrix or data frame of covariate values.
#' @return A vector of estimated treatment propensities P(A = 1|X) in [0, 1].
estimate_pi_par <- function(A, X) {
  dat <- as.data.frame(cbind(A, X))
  form_pi <- as.formula(paste("A ~", paste(colnames(X), collapse = " + ")))
  pihat <- glm(form_pi, family = binomial(link = "logit"), dat)$fitted.values

  return(pihat)
}


#' Estimate outcome regression function parametrically, via linear regression.
#'
#' @param Y A vector of continuous outcomes of length n.
#' @param A A vector of treatment values of length n, in {0, 1}.
#' @param X An n-by-p matrix or data frame of covariate values.
#' @return An n-by-2 data frame of estimated regression values E[Y|A=a, X], for
#' a = 0 (column 1) and a = 1 (column 2).
estimate_mu_par <- function(Y, A, X) {
  dat <- as.data.frame(cbind(Y, A, X))
  form_mu <- as.formula(paste("Y ~ A +", paste(colnames(X), collapse = " + ")))
  mumod <- glm(form_mu, family = gaussian, dat)
  dat$A <- 0; mu0 <- predict(mumod, newdata = dat)
  dat$A <- 1; mu1 <- predict(mumod, newdata = dat)
  mu <- A*mu1 + (1 - A)*mu0
  muhat <- as.data.frame(cbind(mu, mu0, mu1))

  return(muhat)
}


#' Estimate propensity scores nonparametrically using SuperLearner.
#'
#' @param A A vector of treatment values, in {0, 1}.
#' @param X A matrix or data frame of covariate values.
#' @param sl_lib A character vector of learner functions created via
#' SuperLearner, as in sl_lib <- SuperLearner::create.Learner(...)$names, where
#' ... are the learner type and learner parameters. These are passed to
#' SuperLearner() to estimate the propensity score function.
#' @param nsplits The number of splits to use for sample splitting, if fold
#' indices aren't provided to the split_inds argument.
#' @param split_inds Fold indices for sample splitting, if nsplits is NULL. One
#' of split_inds or nsplits must be provided. split_inds takes precedence if
#' it's not NULL.
#' @param ... Additional arguments passed to SuperLearner.
#' @return A vector of estimated treatment propensities P(A = 1|X) in [0, 1].
estimate_pi_nonpar <- function(A, X, sl_lib, nsplits=NULL, split_inds=NULL, ...) {
  if (is.null(nsplits) & is.null(split_inds)) {
    stop("One of nsplits or split_inds must be provided.")
  }
  if (is.null(split_inds)) {
    split_inds <- generate_folds(length(A), nsplits)
  }
  else {
    nsplits = length(unique(split_inds))
  }
  X <- as.data.frame(X)
  n <- length(A)
  pihat <- rep(NA, n)
  coefs <- rep(0, length(sl_lib))

  for (vfold in 1:nsplits) {
    train <- split_inds != vfold
    test <- split_inds == vfold
    if (nsplits == 1) {
      train <- test
    }
    model <- SuperLearner(A[train], X[train, ], newX = X[test, ],
                          family=binomial(link = "logit"),
                          SL.library = sl_lib, ...)
    pihat[test] <- model$SL.predict
    coefs <- coefs + model$coef
  }
  coefs <- coefs/nsplits

  return(list(preds = pihat, coefs = coefs))
}


#' Estimate outcome regression function nonparametrically using SuperLearner.
#'
#' Outcomes are treated as continuous, meaning that the error distribution in
#' SuperLearner() is set to "gaussian."
#' @param Y A vector of continuous outcomes of length n.
#' @param A A vector of treatment values, in {0, 1}.
#' @param X A matrix or data frame of covariate values.
#' @param sl_lib A character vector of learner functions created via
#' SuperLearner, as in sl_lib <- SuperLearner::create.Learner(...)$names, where
#' ... are the learning type and learner parameters. These are passed to
#' SuperLearner() to estimate the propensity score function.
#' @param nsplits The number of splits to use for sample splitting, if fold
#' indices aren't provided to the split_inds argument.
#' @param split_inds Fold indices for sample splitting, if nsplits is NULL. One
#' of split_inds or nsplits must be provided. split_inds takes precedence if
#' it's not NULL.
#' @param ... Additional arguments passed to SuperLearner.
#' @return An n-by-3 data frame of estimated regression values, where column
#' 1 = E[Y|A, X], column 2 = E[Y|A=0, X], and column 3 = E[Y|A=1, X].
#' a = 0 (column 1) and a = 1 (column 2).
estimate_mu_nonpar <- function(Y, A, X, sl_lib, nsplits=NULL, split_inds=NULL, ...) {
  if (is.null(nsplits) & is.null(split_inds)) {
    stop("One of nsplits or split_inds must be provided.")
  }
  if (is.null(split_inds)) {
    split_inds <- generate_folds(length(A), nsplits)
  }
  else {
    nsplits = length(unique(split_inds))
  }
  XA <- as.data.frame(cbind(X, A))
  X <- as.data.frame(X)
  n <- length(A)
  muhat <- as_tibble(matrix(0, nrow = n, ncol = 3,
                            dimnames = list(c(), c("mu0", "mu1", "mu"))))
  coefs <- rep(0, length(sl_lib))

  for (vfold in 1:nsplits) {
    train <- split_inds != vfold
    test <- split_inds == vfold
    if (nsplits == 1) {
      train <- test
    }
    model0 <- SuperLearner(Y[(A == 0) & train], X[(A == 0) & train, ],
                           newX = X[test, ], family = gaussian, SL.library = sl_lib, ...)
    muhat[test, "mu0"] <- model0$SL.predict
    coefs <- coefs + model0$coef
    model1 <- SuperLearner(Y[(A == 1) & train], X[(A == 1) & train, ],
                           newX = X[test, ], family = gaussian, SL.library = sl_lib, ...)
    muhat[test, "mu1"] <- model1$SL.predict
    coefs <- coefs + model1$coef
  }
  coefs <- coefs/(2*nsplits)
  muhat <- muhat %>%
    mutate(mu = A*mu1 + (1 - A)*mu0)

  return(list(preds = muhat, coefs = coefs))
}


#' Bound estimated propensity scores away from 0 and 1.
#'
#' @param pi Vector of estimated propensity scores in [0, 1].
#' @param lower Lower bound.
#' @param upper Upper bound.
#' @return The bounded scores.
bound_pi <- function(pi, lower=0.025, upper=0.975) {
  pi[pi < lower] <- lower
  pi[pi > upper] <- upper
  return(pi)
}


#' Compute regression, aka g-computation, estimate of the ATE.
#'
#' @param muhat 2-column data frame of estimated regression values, where
#' column 1 = E[Y|A=0, X] and column 2 = E[Y|A=1, X]. The columns must be named
#' mu0 and mu1.
#' @return A named vector with elements "est": the estimated ATE using the
#' g-computation estimator; and "se", which is always NA and which is a
#' placeholder for the bootstrap-estimated standard error. This element is
#' included to make the output of this function equivalent to the output from
#' other estimator functions.
est_gcomp <- function(muhat) {
  est <- mean(muhat$mu1 - muhat$mu0)
  return(c("est" = est, "se" = NA))
}


#' Compute regression, aka g-computation, estimate, via linear regression.
#'
#' This function is specifically designed to be passed to boot::boot() when
#' estimating the standard error via bootstrapping.
#' @param dat Data frame containing covariates, treatments, and outcomes.
#' @param rows Which rows of the data to include. If NULL, then all rows are
#' included. This argument is required in order to pass the function to boot().
#' @param covs String representing the names of the covariate columns, which
#' should be numbered starting with 1. For example, the default "X" means that
#' the covariate columns are named X1, X2, etc.
#' @param treat The name of the treatment column.
#' @param outcome The name of the outcome column.
#' @return A named vector with elements "est": the estimated ATE using the
#' g-computation estimator; and "se", which is always NA and which is a
#' placeholder for the bootstrap-estimated standard error. This element is
#' included to make the output of this function equivalent to the output from
#' other estimator functions.
est_gcomp_parametric <- function(dat, rows=NULL, covs="X", treat="A", outcome="Y") {
  if (!(is.null(rows))) {dat <- dat[rows, ]}
  covcols <- colnames(dat)[startsWith(colnames(dat), covs)]
  if (length(covcols) == 0) {
    stop(paste("There are no columns in the data that start with", covs))
  }
  form <- as.formula(paste(outcome, "~", treat, "+", paste(covcols, collapse = "+")))
  mumod <- glm(form, data=dat)
  dat[[treat]] <- 1; muhat1 <- predict(mumod, newdata=dat)
  dat[[treat]] <- 0; muhat0 <- predict(mumod, newdata=dat)
  est <- mean(muhat1 - muhat0)

  return(c("est" = est, "se" = NA))
}


#' Compute g-computation estimate of ATE with nonparametric outcome model.
#'
#' This function is specifically designed to be passed to boot::boot() when
#' estimating the standard error via bootstrapping.
#' @param dat Data frame containing covariates, treatments, and outcomes.
#' @param rows Which rows of the data to include. If NULL, then all rows are
#' included. This argument is required in order to pass the function to boot().
#' @param sl_lib A character vector of learner functions created via
#' SuperLearner, as in sl_lib <- SuperLearner::create.Learner(...)$names, where
#' ... are the learner type and learner parameters. These are passed to
#' SuperLearner() to estimate the propensity score function.
#' @param covs String representing the names of the covariate columns, which
#' should be numbered starting with 1. For example, the default "X" means that
#' the covariate columns are named X1, X2, etc.
#' @param treat The name of the treatment column.
#' @param outcome The name of the outcome column.
#' @param nsplits The number of splits to use for sample splitting, if fold
#' indices aren't provided to the split_inds argument.
#' @param split_inds Fold indices for sample splitting, if nsplits is NULL. One
#' of split_inds or nsplits must be provided. split_inds takes precedence if
#' it's not NULL.
#' @param ... Additional arguments passed to SuperLearner.
#' @return A named vector with elements "est": the estimated ATE using the
#' g-computation estimator; and "se", which is always NA and which is a
#' placeholder for the bootstrap-estimated standard error. This element is
#' included to make the output of this function equivalent to the output from
#' other estimator functions.
est_gcomp_nonparametric <- function(dat, rows=NULL, sl_lib, covs="X",
                                    treat="A", outcome="Y", nsplits=NULL,
                                    split_inds=NULL, ...) {
  if (!(is.null(rows))) {dat <- dat[rows, ]}
  covcols <- colnames(dat)[startsWith(colnames(dat), covs)]
  if (length(covcols) == 0) {
    stop(paste("There are no columns in the data that start with", covs))
  }
  if (is.null(nsplits) & is.null(split_inds)) {
    stop("One of nsplits or split_inds must be provided.")
  }
  if (is.null(split_inds)) {
    split_inds <- generate_folds(length(A), nsplits)
  }
  else {
    nsplits = length(unique(split_inds))
  }

  X <- dat[covcols]
  A <- dat[[treat]]
  Y <- dat[[outcome]]
  n <- length(A)
  muhat0 <- rep(NA, n)
  muhat1 <- rep(NA, n)

  for (vfold in 1:nsplits) {
    train <- split_inds != vfold
    test <- split_inds == vfold
    if (nsplits == 1) {
      train <- test
    }
    model0 <- SuperLearner(Y[(A == 0) & train], X[(A == 0) & train, ],
                           newX = X[test, ], family = gaussian, SL.library = sl_lib, ...)
    muhat0[test] <- model0$SL.predict
    model1 <- SuperLearner(Y[(A == 1) & train], X[(A == 1) & train, ],
                           newX = X[test, ], family = gaussian, SL.library = sl_lib, ...)
    muhat1[test] <- model1$SL.predict
  }

  est <- mean(muhat1 - muhat0)

  return(c("est" = est, "se" = NA))
}


#' Compute IPW estimate of ATE.
#'
#' @param Y A vector of continuous outcomes of length n.
#' @param A A vector of treatment values, in {0, 1}.
#' @param pihat A vector of estimated propensity scores. Estimation can be poor
#' when values are close to 0 or 1, so it may be advisable to bound these in,
#' say, [0.025, 0.975].
#' @return A named vector with elements "est", the estimated ATE, and "se", the
#' standard error, estimated analytically.
est_ipw <- function(Y, A, pihat) {
  weights <- A * (mean(A) / pihat) + (1 - A) * ((1 - mean(A)) / (1 - pihat))
  mod <- lm(Y ~ A, weights = weights)
  est <- coef(mod)[2]
  se <- sqrt(vcovHC(mod, type = "HC")[2,2])

  return(c("est" = unname(est), "se" = se))
}


#' Compute AIPW estimate of ATE.
#'
#' @param Y A vector of continuous outcomes of length n.
#' @param A A vector of treatment values, in {0, 1}.
#' @param pihat A vector of estimated propensity scores. Estimation can be poor
#' when values are close to 0 or 1, so it may be advisable to bound these in,
#' say, [0.025, 0.975].
#' @param muhat 2-column data frame of estimated regression values, where
#' column 1 = E[Y|A=0, X] and column 2 = E[Y|A=1, X]. The columns must be named
#' mu0 and mu1.
#' @return A named vector with elements "est", the estimated ATE, and "se", the
#' standard error, estimated analytically.
est_aipw <- function(Y, A, pihat, muhat) {
  vec <- (((2*A - 1)*(Y - muhat$mu))/
            ((2*A - 1)*pihat + (1 - A)) + muhat$mu1 - muhat$mu0)
  est <- mean(vec)
  se <- sd(vec)/sqrt(length(Y))

  return(c("est" = unname(est), "se" = se))
}


#' Compute tMLE estimate of ATE using tmle::tmle().
#'
#' @param Y A vector of continuous outcomes of length n.
#' @param A A vector of treatment values, in {0, 1}.
#' @param X n-by-p matrix of covariate data.
#' @param pihat A vector of estimated propensity scores. Estimation can be poor
#' when values are close to 0 or 1, so it may be advisable to bound these in,
#' say, [0.025, 0.975].
#' @param muhat 2-column data frame of estimated regression values, where
#' column 1 = E[Y|A=0, X] and column 2 = E[Y|A=1, X]. The columns must be named
#' mu0 and mu1.
#' @return A named vector with elements "est", the estimated ATE, and "se", the
#' standard error, estimated analytically.
est_tmle <- function(Y, A, X, muhat, pihat) {
  Q <- cbind(muhat$mu0, muhat$mu1)
  mod <- tmle(Y, A, X, Q = Q, g1W = pihat)
  est <- mod$estimates$ATE$psi
  se <- sqrt(mod$estimates$ATE$var.psi)

  return(c("est" = unname(est), "se" = se))
}


#' Calculate performance metrics for an estimate of the ATE.
#'
#' The metrics calculated are the bias, squared error, and the coverage and
#' width of a 95% confidence interval.
#' @param est A named vector with elements "est", the estimated ATE, and "se",
#' the estimated standard error.
#' @param effect_size The true ATE, a number.
#' @return A named vector with elements "bias," "mse," "cov" (coverage), and
#' "width". The confidence interval is computed as est +/- 1.96*se. The coverage
#' is 0 or 1, indicating whether the true ATE falls in the interval.
calc_metrics <- function(est, effect_size) {
  effect_est <- est["est"]
  se_est <- est["se"]
  bias <- effect_est - effect_size
  mse <- (effect_est - effect_size)^2
  coverage <- (effect_est - 1.96*se_est < effect_size) & (effect_size < effect_est + 1.96*se_est)
  coverage <- +coverage
  width <- (effect_est + 1.96*se_est) - (effect_est-1.96*se_est)

  return(c(bias = unname(bias), mse = unname(mse),
           cov = unname(coverage), width = unname(width)))
}


#' Compare performance of ATE estimators via simulation.
#'
#' @param counter Monte Carlo iteration number. Used to set up randomization,
#' via set.seed(), and to name output files.
#' @param N The number of Monte Carlo samples to generate.
#' @param beta_a Vector of logistic coefficients of length p, relating the
#' propensity scores to the covariates. Passed to generate_covariates().
#' @param beta_y Vector of linear coefficients of length p, relating the outcome
#' Y to the covariates. Passed to generate_covariates().
#' @param effect_size The treatment effect, which is constant.
#' @param interactions Logical, whether or not to include interaction terms in
#' the  design matrices for both the baseline and transformed covariates. FALSE
#' by default.
#' @param estimators Which estimators to use for nonparametric estimation of the
#' nuisance parameters, some subset of c("aipw", "gcomp", "ipw", "tmle"). If
#' nonpar == FALSE, then this argument is ignored.
#' @param nonpar Logical, whether to estimate the nuisance parameters
#' nonparametrically (if TRUE, via specified learners) or parametrically (if
#' FALSE, via logistic regression).
#' @param nsplits The number of splits to use for sample splitting if
#' estimating nuisance parameters nonparametrically.
#' @param sl_lib_mu Learners to pass to SuperLearner to estimate the outcome
#' regression. Should generally be a vector of strings that refer to functions
#' created by SuperLearner::create.Learner(). For example, this argument could
#' be c(ranger_learner$names) after having called create.Learner("SL.ranger").
#' If nonpar == FALSE, then this argument is ignored.
#' @param sl_lib_pi Learners to pass to SuperLearner to estimate the propensity
#' regression. If nonpar == FALSE, then this argument is ignored.
#' @param cv_folds The number of folds for Superlearner to use internally. Each
#' split from 1 to nsplits is passed to SuperLearner, which uses cv_folds for
#' its learning procedure.
#' @param bootstrap_reg Whether to estimate the variance of the
#' g-computation-based estimator. The other estimators all have closed-form
#' asymptotic variances. This is computationally expensive, especially if
#' nonpar == TRUE, so it is set to FALSE by default.
#' @param bootNum The number of bootstrap iterations to use when estimating the
#' variance of the g-computation-based estimator. If bootstrap_reg == FALSE,
#' then this argument is ignored.
#' @param oracle_estimators Which estimators to use to generate oracle
#' estimates, some subset of c("aipw", "gcomp", "ipw", "tmle"). Oracle estimates
#' use the true propensity scores and outcome regression values. This is useful
#' for benchmarking the performance of real estimators. If NULL or "", then no
#' oracle estimates are computed.
#' @param remainder Whether to compute (an empirical estimate of) the remainder
#' terms for the doubly robust estimators. Two remainders are computed, one for
#' the baselines covariates X and one for the transformed covariates Z. The
#' doubly robust estimators are based on a von Mises expansion that consists of
#' an influence function term plus a remainder term. The remainder term
#' determines how quickly the doubly robust estimators converge to the truth.
#' The remainder term is "second order," meaning that it involves a product of
#' the error in the propensity estimator and the error in the outcome regression
#' estimator, so the doubly robust estimators are root-n consistent if, for
#' example, both the nuisance parameter estimators are consistent at faster than
#' n^(1/4) rates. Checking whether the remainder term is shrinking at the
#' expected rate can therefore be useful for diagnostic purposes.
#' @return A list of simulation results, which includes the data generated,
#' a dataframe of the true and estimated nuisance parameter values, the
#' remainders if remainder == TRUE, metrics for the estimators (the bias and
#' mse, as well as coverage and width if CIs were computed), and, if
#' nonpar == TRUE, the coefficients for the learners in the SuperLearner
#' ensemble for both the propensity and outcome regressions.
simulate_inference <- function(counter, N, beta_a, beta_y, effect_size,
                               interactions=FALSE,
                               estimators=c("aipw", "gcomp", "ipw", "tmle"),
                               nonpar=FALSE, nsplits=1, sl_lib_mu=NULL, sl_lib_pi=NULL,
                               cv_folds=5, bootstrap_reg=FALSE, bootNum=100,
                               oracle_estimators=NULL, remainder=FALSE) {
  set.seed(counter)
  estimators <- sort(tolower(estimators))
  print(paste("...Estimating effects for these estimators:",
              paste(estimators, collapse = ", ")))
  est <- list()

  ## Generate covariates X and transformed covariates Z, which induce misspecification
  covs <- generate_covariates(N)
  X_base <- covs[["X"]]
  Z_base <- covs[["Z"]]
  if (interactions) {
    X <- model.matrix(~(X1 + X2 + X3 + X4)^2, data = as.data.frame(X_base))[, -1]
    Z <- model.matrix(~(Z1 + Z2 + Z3 + Z4)^2, data = as.data.frame(Z_base))[, -1]
  }
  else {
    X <- model.matrix(~X1 + X2 + X3 + X4, data = as.data.frame(X_base))[, -1]
    Z <- model.matrix(~Z1 + Z2 + Z3 + Z4, data = as.data.frame(Z_base))[, -1]
  }
  colnames(X) <- paste0("X", 1:ncol(X))
  colnames(Z) <- paste0("Z", 1:ncol(Z))

  ## Generate treatments and outcomes
  obs <- generate_outcomes(X_base, beta_a, beta_y)
  A <- obs[["A"]]
  Y <- obs[["Y"]]

  ## Create combined dataset
  dat <- bind_cols(tbl_df(X), tbl_df(Z), obs)

  ## Estimate nuisance functions using X ("_T") and Z ("_F")
  pihat_coefs <- NULL
  muhat_coefs <- NULL
  if (nonpar) {
    # Make sure pihat and muhat use the same folds for cross-fitting
    inds <- generate_folds(N, nsplits)

    pihat_T_list <- estimate_pi_nonpar(A, X, sl_lib_pi, split_inds = inds,
                                       cvControl = list(V = cv_folds))
    pihat_T <- pihat_T_list$preds %>% bound_pi
    pihat_F_list <- estimate_pi_nonpar(A, Z, sl_lib_pi, split_inds = inds,
                                       cvControl = list(V = cv_folds))
    pihat_F <- pihat_F_list$preds %>% bound_pi
    pihat_coefs <- bind_cols(learner = names(pihat_T_list$coefs),
                             pihat_T = pihat_T_list$coefs, pihat_F = pihat_F_list$coefs)

    muhat_T_list <- estimate_mu_nonpar(Y, A, X, sl_lib_mu, split_inds = inds,
                                       cvControl = list(V = cv_folds))
    muhat_T <- muhat_T_list$preds
    muhat_F_list <- estimate_mu_nonpar(Y, A, Z, sl_lib_mu, split_inds = inds,
                                       cvControl = list(V = cv_folds))
    muhat_F <- muhat_F_list$preds
    muhat_coefs <- bind_cols(learner = names(muhat_T_list$coefs),
                             muhat_T = muhat_T_list$coefs, muhat_F = muhat_F_list$coefs)
  }
  else {
    pihat_T <- estimate_pi_par(A, X) %>% bound_pi
    pihat_F <- estimate_pi_par(A, Z) %>% bound_pi
    muhat_T <- estimate_mu_par(Y, A, X)
    muhat_F <- estimate_mu_par(Y, A, Z)
  }
  pi_true <- 1 - expit(cbind(1, X_base) %*% beta_a)
  mu0 <- cbind(1, X_base) %*% beta_y
  mu1 <- mu0 + effect_size
  mu = A*mu1 + (1 - A)*mu0
  mu_true <- data.frame(mu0 = mu0, mu1 = mu1, mu = mu)
  predictions = data.frame(pi_true = pi_true, pihat_T = pihat_T, pihat_F = pihat_F,
                           mu0_true = mu0, mu1_true = mu1, mu_true = mu,
                           muhat0_T = muhat_T$mu0, muhat1_T = muhat_T$mu1, muhat_T = muhat_T$mu,
                           muhat0_F = muhat_F$mu0, muhat1_F = muhat_F$mu1, muhat_F = muhat_F$mu1) %>%
    as_tibble

  ## Estimate ATE and se for each chosen estimator
  if ("gcomp" %in% estimators) {
    est$gcomp_T <- est_gcomp(muhat_T)
    est$gcomp_F <- est_gcomp(muhat_F)
  }
  if ("ipw" %in% estimators) {
    est$ipw_T <- est_ipw(Y, A, pihat_T)
    est$ipw_F <- est_ipw(Y, A, pihat_F)
  }
  if ("aipw" %in% estimators) {
    est$aipw_T <- est_aipw(Y, A, pihat_T, muhat_T)
    est$aipw_F <- est_aipw(Y, A, pihat_F, muhat_F)
  }
  if ("tmle" %in% estimators) {
    est$tmle_T <- est_tmle(Y, A, X, muhat_T, pihat_T)
    est$tmle_F <- est_tmle(Y, A, Z, muhat_F, pihat_F)
  }

  ## Compute bootstrap SEs for regression estimator
  if (("gcomp" %in% estimators) && bootstrap_reg) {
    datT <- bind_cols(as_tibble(X), A = A, Y = Y)
    datF <- bind_cols(as_tibble(Z), A = A, Y = Y)
    if (nonpar) {
      bs_T <- boot(datT, est_gcomp_nonparametric, R=bootNum, covs = "X",
                   split_inds = inds, sl_lib = sl_lib_mu,
                   cvControl = list(V = cv_folds))
      bs_F <- boot(datF, est_gcomp_nonparametric, R=bootNum, covs = "Z",
                   split_inds = inds, sl_lib = sl_lib_mu,  cvControl = list(V = cv_folds))
    }
    else {
      bs_T <- boot(datT, est_gcomp_parametric, R=bootNum, covs="X")
      bs_F <- boot(datF, est_gcomp_parametric, R=bootNum, covs="Z")
    }
    est[["gcomp_T"]]["se"] <- sd(bs_T$t[, 1])
    est[["gcomp_F"]]["se"] <- sd(bs_F$t[, 1])
  }

  ## Add suffix to indicate parametric or nonparametric
  if (nonpar) {
    names(est) <- paste0(names(est), "NP")
  }
  else {
    names(est) <- paste0(names(est), "PM")
  }

  ## Compute oracle predictions
  if ("gcomp" %in% oracle_estimators) {
    est$gcomp_oracle <- est_gcomp(mu_true)
  }
  if ("ipw" %in% oracle_estimators) {
    est$ipw_oracle <- est_ipw(Y, A, pi_true)
  }
  if ("aipw" %in% oracle_estimators) {
    est$aipw_oracle <- est_aipw(Y, A, pi_true, mu_true)
  }
  if ("tmle" %in% oracle_estimators) {
    est$tmle_oracle <- est_tmle(Y, A, X, mu_true, pi_true)
  }

  ## Calculate bias, MSE, coverage, and width of confidence intervals
  performance <- sapply(est, calc_metrics, effect_size = effect_size)
  metrics <- rbind(data.frame(est), performance)
  metrics <- cbind(counter = counter, N = N, metrics)
  metrics <- rownames_to_column(metrics, var = "metric")

  ## Compute remainder term for diagnostic purposes
  if (remainder) {
    R2_T <- mean((pi_true - pihat_T)/pihat_T*(muhat_T$mu1 - mu1)) -
      mean(((1 - pi_true) - (1 - pihat_T))/(1 - pihat_T)*(muhat_T$mu0 - mu0))
    R2_F <- mean((pi_true - pihat_F)/pihat_F*(muhat_F$mu1 - mu1)) -
      mean(((1 - pi_true) - (1 - pihat_F))/(1 - pihat_F)*(muhat_F$mu0 - mu0))
    remainder <- data.frame(N = N, R2_T = R2_T, R2_F = R2_F)
  }
  else {
    remainder <- NULL
  }
  out <- list(data = dat, predictions = predictions, remainder = remainder,
              metrics = metrics, coefs = list(pihat_coefs = pihat_coefs,
                                              muhat_coefs = muhat_coefs))
  invisible(out)
}


#' Save simulation results.
#'
#' @param results The output of simulate_inference().
#' @param counter Monte Carlo iteration number. Used to name output files.
#' @param estimators The estimators that are included in the results, used to
#' name output files.
#' @param outdir Directory to write results to.
#' @param include Which elements of the results to save, either "all" or some
#' subset of c("data", "predictions", "remainder", "metrics", "coefs"), which
#' are the components of the results. Nothing will be written for components
#' that are NULL.
#' @param suffix An additional optional suffix to add to the filenames.
#' @return NULL. As a side effect, this function writes results to outdir. Each
#' included component gets its own subdirectory.
write_results <- function(results, counter, estimators, outdir, include="all",
                          suffix="") {
  if (!dir.exists(outdir)) {dir.create(outdir)}
  if (include[1] == "all") {
    include <- c("data", "predictions", "remainder", "metrics", "coefs")
  }
  # Write data
  if (("data" %in% include) && (!is.null(results$data))) {
    subdir <- file.path(outdir, "data")
    if (!dir.exists(subdir)) {dir.create(subdir)}
    path <- file.path(subdir, paste0("iteration", counter, "_data", suffix, ".txt"))
    write_tsv(results[["data"]], path)
  }
  if (("predictions" %in% include) && (!is.null(results$predictions))) {
    subdir <- file.path(outdir, "predictions")
    if (!dir.exists(subdir)) {dir.create(subdir)}
    path <- file.path(subdir, paste0("iteration", counter, "_predictions", suffix, ".txt"))
    write_tsv(results[["predictions"]], path)
  }
  if (("remainder" %in% include) && (!is.null(results$remainder))) {
    subdir <- file.path(outdir, "remainder")
    if (!dir.exists(subdir)) {dir.create(subdir)}
    path = file.path(subdir, paste0("iteration", counter, "_remainder", suffix, ".txt"))
    write_tsv(results[["remainder"]], path)
  }
  if (("metrics" %in% include) && (!is.null(results$metrics))) {
    subdir <- file.path(outdir, "metrics")
    if (!dir.exists(subdir)) {dir.create(subdir)}
    path = file.path(subdir, paste0("iteration", counter, "_",
                                    paste(estimators, collapse = "_"), suffix, ".txt"))
    write_tsv(results[["metrics"]], path)
  }
  if (("coefs" %in% include) && (!is.null(results$coefs[[1]]))) {
    subdir <- file.path(outdir, "superlearner_coefs")
    if (!dir.exists(subdir)) {dir.create(subdir)}
    pi_path = file.path(subdir, paste0("iteration", counter, "_pihat", suffix, ".txt"))
    write_tsv(results[["coefs"]][["pihat_coefs"]], pi_path)
    mu_path = file.path(subdir, paste0("iteration", counter, "_muhat", suffix, ".txt"))
    write_tsv(results[["coefs"]][["muhat_coefs"]], mu_path)
  }
}
