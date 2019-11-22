#' Expit function, aka inverse logit.
expit <- function(x){exp(x)/(1+exp(x))}

#' Generate indices for CV folds.
#'
#' @param n Number of rows in the dataset.
#' @param cv_folds Number of folds to create.
#' @param probs Probabilities associated with each fold. If NULL, then data
#' points are assigned to each fold with equal probability. Fold memberships are
#' iid sampled for each data point.
#' @return A list of row indices corresponding to each fold.
generate_folds <- function(n, cv_folds, probs=NULL) {
  inds <- sample(1:cv_folds, size = n, replace = TRUE, prob = probs)
  out <- lapply(1:cv_folds, function(x) which(inds == x))
  return(out)
}

#' Generate raw and transformed covariates.
#'
#' The raw covariates are a 4-dimensional multivariate normal with mean 0 and
#' diagonal covariance. They are transformed following Kang \& Schafer, 2007,
#' to induce misspecification in the parametric propensity and outcome models.
#' @param n Number of samples.
#' @return A list containing two n-by-4 matrices: X, the raw covariates, and Z,
#' the transformed version of X.
generate_covariates <- function(n) {
  X <- rmvnorm(n, mean = rep(0, 4), sigma = diag(4))
  colnames(X) <- paste0("X", 1:4)
  Z <- X
  Z[,1] <- exp(X[,1]/2)
  Z[,2] <- X[,2]/(1 + exp(X[,1])) + 10
  Z[,3] <- (X[,1]*X[,3]/25 + 0.6)^3
  Z[,4] <- (X[,2] + X[,4] + 20)^2
  colnames(Z) <- paste0("Z", 1:4)

  return(list(X = X, Z = Z))
}

#' Generate treatment levels and outcomes from covariates.
#'
#' Propensity scores \eqn{\pi(X) := P(A = 1|X)} are a logistic function of the
#' covariates X. The treatment A is binary, 0 or 1. The outcomes
#' \eqn{\mu(X, a) := E[Y|X, A = a]}, for \eqn{a \in {0, 1}}, are linear in
#' the covariates, plus independent Gaussian noise, plus a constant treatment
#' effect for those who received treatment. Given this structure, the treatment
#' and the potential outcomes are independent conditional on X, but not
#' marginally because X confounds the treatment and the outcome.
#' @param X n-by-p matrix of covariate data.
#' @param beta_a Vector of logistic coefficients of length p, relating the
#' propensity scores to X by the equation \eqn{\pi(X) = expit(X^T \beta_a)}.
#' @param beta_y Vector of linear coefficients of length p, relating the outcome
#' Y to X by the equation \eqn{\mu(X, a) = X^T \beta_y}.
#' @param effect_size The treatment effect, which is constant.
#' @return A list containing the treatments A and the outcomes Y, each of length
#' n.
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
#' @return A vector of estimated propensities in [0, 1].
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


#' Estimate propensity scores nonparametrically using a superlearner.
#'
#' @param A A vector of treatment values, in {0, 1}.
#' @param X A matrix or data frame of covariate values.
#' @param sl A superlearner, an object of class Lrnr_sl from the \code{sl3}
#' library.
#' @param split Logical, whether to perform sample splitting. If TRUE, then the
#' sample is randomly split in 2; learning is performed one part, and
#' propensities are estimated on the other. The process is repeated with the
#' samples switched so that a full vector of estimated propensities is
#' generated. If FALSE, then learning and estimation are both performed on the
#' entire sample.
#' @param ... Additional arguments passed to SuperLearner.
#' @return A vector of estimated propensities in [0, 1].
estimate_pi_nonpar <- function(A, X, sl, ...) {
  dat <- bind_cols(tbl_df(X), A = A)
  task <- make_sl3_Task(data = dat, covariates = colnames(X), outcome = "A",
                        outcome_type = "binomial", ...)
  sl_fit <- sl$train(task)
  out <- list(preds = sl_fit$predict(), coefs = sl_fit$coefficients)
  return(out)
}


#' Estimate outcome regression function nonparametrically using a superlearner.
#'
#' Outcomes are treated as continuous.
#' @param Y A vector of continuous outcomes of length n.
#' @param A A vector of treatment values, in {0, 1}, of length n.
#' @param X An n-by-p matrix or data frame of covariate values.
#' @param sl A superlearner, an object of class Lrnr_sl from the \code{sl3}
#' library.
#' @param ... Additional arguments passed to \code{sl3::make_sl3_Task}.
#' @return An n-by-3 data frame of estimated regression values, where column
#' 1 = \eqn{E[Y|A, X]}, column 2 = \eqn{E[Y|A=0, X]}, and column 3 =
#' \eqn{E[Y|A=1, X]}.
estimate_mu_nonpar <- function(Y, A, X, sl, ...) {
  dat <- bind_cols(tbl_df(X), A = A, Y = Y)
  task <- make_sl3_Task(data = dat, covariates = c(colnames(X), "A"),
                        outcome = "Y", outcome_type = "continuous", ...)
  sl_fit <- sl$train(task)
  mu <- sl_fit$predict()
  task0 <- make_sl3_Task(data = dat %>% mutate(A = 0),
                         covariates = c(colnames(X), "A"), outcome = "Y",
                         outcome_type = "continuous", ...)
  mu0 <- sl_fit$predict(task = task0)
  task1 <- make_sl3_Task(data = dat %>% mutate(A = 1),
                         covariates = c(colnames(X), "A"), outcome = "Y",
                         outcome_type = "continuous", ...)
  mu1 <- sl_fit$predict(task = task1)
  out <- list(preds = data.frame(mu = mu, mu0 = mu0, mu1 = mu1),
              coefs = sl_fit$coefficients)
  return(out)
}

#' Bound estimated propensity scores away from 0 and 1.
#'
#' @param pi Vector of estimated propensity scores in [0, 1].
#' @param lower Lower bound.
#' @param upper Upper bound.
#' @return The scores, bounded in [lower, upper].
bound_pi <- function(pi, lower=0.025, upper=0.975) {
  pi[pi < lower] <- lower
  pi[pi > upper] <- upper
  return(pi)
}


#' Compute regression, aka g-computation, estimate of ATE.
#'
#' @param muhat 2-column data frame of estimated regression values, where
#' column 1 = \eqn{\hat{E}[Y|A=0, X]} and column 2 = \eqn{\hat{E}[Y|A=1, X]}. The
#' columns must be named \code{mu0} and \code{mu1}.
#' @return Estimated ATE, the mean of \code{mu1} - \code{mu0}.
estimate_gcomp <- function(muhat) {
  est <- mean(muhat$mu1 - muhat$mu0)
  return(c("est" = est, "se" = NA))
}


#' Compute regression, aka g-computation, estimate, via linear regression.
#'
#' This function is specifically designed to be passed to \code{boot::boot()}
#' when estimating the standard error via bootstrapping.
#' @param dat Data frame containing covariates, treatments, and outcomes.
#' @param rows Which rows of the data to include. If NULL, then all rows are
#' included. This argument is required in order to pass the function to
#' \code{boot()}.
#' @param covs String representing the names of the covariate columns, which
#' should be numbered starting with 1. For example, the default "X" means that
#' the covariate columns are named \code{X1}, \code{X2}, etc.
#' @param treat The name of the treatment column.
#' @param outcome The name of the outcome column.
#' @return A named vector with elements "est": the estimated ATE using the
#' g-computation estimator; and "se", which is always NA and which is a
#' placeholder for the bootstrap-estimated standard error. This element is
#' included to make the output of this function equivalent to the output from
#' other estimator functions. The se must be estimated by bootstrapping, using
#' this function, and then computing the empirical standard error over the
#' bootstrap samples.
estimate_gcomp_par <- function(dat, rows=NULL, covs="X", treat="A", outcome="Y") {
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
#' @param SL.library Vector of names of learners to use when estimating the
#' outcome model, passed to SuperLearner().
#' @param cv_folds Number of cross-validation folds to use, passed to the V
#' argument of the cvControl list in SuperLearner().
#' @param treat The name of the treatment column.
#' @param outcome The name of the outcome column.
#' @return A named vector with elements "est": the estimated ATE using the
#' g-computation estimator; and "se", which is always NA and which is a
#' placeholder for the bootstrap-estimated standard error. This element is
#' included to make the output of this function equivalent to the output from
#' other estimator functions. The se must be estimated by bootstrapping, using
#' this function, and then computing the empirical standard error over the
#' bootstrap samples.
est_gcomp_nonparametric <- function(dat, rows=NULL, sl, treat="A",
                                    outcome="Y", ...) {
  if (!(is.null(rows))) {dat <- dat[rows, ]}
  covs <- names(dat)[names(dat) != outcome]
  task <- make_sl3_Task(data = dat, covariates = covs, outcome = outcome,
                        outcome_type = "continuous", ...)
  sl_fit <- sl$train(task)
  dat0 <- dat; dat0[[treat]] <- 0
  task0 <- make_sl3_Task(data = dat0, covariates = covs, outcome = outcome,
                         outcome_type = "continuous", ...)
  mu0 <- sl_fit$predict(task = task0)

  dat1 <- dat; dat1[[treat]] <- 1
  task1 <- make_sl3_Task(data = dat1, covariates = covs, outcome = outcome,
                         outcome_type = "continuous", ...)
  mu1 <- sl_fit$predict(task = task1)

  est <- mean(mu1 - mu0)

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
estimate_ipw <- function(Y, A, pihat) {
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
estimate_aipw <- function(Y, A, pihat, muhat) {
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
#' @param pihat A vector of estimated propensity scores. Estimation can be poor
#' when values are close to 0 or 1, so it may be advisable to bound these in,
#' say, [0.025, 0.975].
#' @param muhat 2-column data frame of estimated regression values, where
#' column 1 = E[Y|A=0, X] and column 2 = E[Y|A=1, X]. The columns must be named
#' mu0 and mu1.
#' @return A named vector with elements "est", the estimated ATE, and "se", the
#' standard error, estimated analytically.
estimate_tmle <- function(Y, A, X, muhat, pihat) {
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
#' @param counter Integer used to set up randomization, via set.seed(), and to
#' name output files.
#' @param N The number of samples to generate.
#' @param beta_a A length 5 vector of coefficients, where the first coefficient
#' is the intercept.
#' @export
simulate_inference <- function(counter, N, beta_a, beta_y, effect_size,
                 estimators=c("aipw", "gcomp", "ipw", "tmle"),
                 nonpar=FALSE, sl_lib_mu=NULL, sl_lib_pi=NULL,
                 cv_folds=5, bootstrap_reg=FALSE, bootNum=100,
                 oracle_estimators=NULL, remainder=FALSE) {
  # set.seed(counter)
  estimators <- sort(tolower(estimators))
  print(paste("...Estimating effects for these estimators:",
              paste(estimators, collapse = ", ")))
  est <- list()
  # est <- as.list(rep(NA, 8))
  # names(est) <- c("aipw_F", "aipw_T", "gcomp_F", "gcomp_T", "iwp_F", "ipw_T",
  #                 "tmle_F", "tmle_T")

  ## Generate covariates X and transformed covariates Z, which induce misspecification
  covs <- generate_covariates(N)
  X_base <- covs[["X"]]
  X <- X_base
  X <- model.matrix(~(-1 + X1 + X2 + X3 + X4)^2, data.frame(X_base))  # 2-way interactions, no intercept
  colnames(X) <- paste0("X", 1:ncol(X))
  Z <- covs[["Z"]]
  Z <- model.matrix(~(-1 + Z1 + Z2 + Z3 + Z4)^2, data.frame(Z))
  colnames(Z) <- paste0("Z", 1:ncol(Z))
  # Z <- as_tibble(Z_base) %>%
  #   mutate(Z1sq_Z2 = Z1^2*Z2, Z5_new = Z3^(1/3)/log(Z1)) %>%
  #   dplyr::select(c(Z1, Z2, Z4, Z1sq_Z2, Z5_new)) %>%
  #   as.matrix
  # colnames(Z_base) <- c("W1", "W2", "W3", "W4")
  # Z <- as_tibble(Z_base) %>%
  #   mutate(Z1 = log(W1), Z2 = 1/log(W1), Z3 = W1^2, Z4 = W2, Z5 = sqrt(W4),
  #          Z6 = W1^2*W2, Z7 = W3^(1/3)/log(W1))

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
    # Make sure pihat and muhat use the same folds
    inds <- generate_folds(N, cv_folds)
    cvControl = list(V = cv_folds, validRows = inds)

    pihat_T_list <- estimate_pi_nonpar(A, X, sl_lib_pi, cvControl = cvControl)
    pihat_T <- pihat_T_list$preds %>% bound_pi
    pihat_F_list <- estimate_pi_nonpar(A, Z, sl_lib_pi, cvControl = cvControl)
    pihat_F <- pihat_F_list$preds %>% bound_pi
    pihat_coefs <- bind_cols(learner = names(pihat_T_list$coefs),
                             pihat_T = pihat_T_list$coefs, pihat_F = pihat_F_list$coefs)

    muhat_T_list <- estimate_mu_nonpar(Y, A, X, sl_lib_mu, cvControl = cvControl)
    muhat_T <- muhat_T_list$preds
    muhat_F_list <- estimate_mu_nonpar(Y, A, Z, sl_lib_mu, cvControl = cvControl)
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
    est$gcomp_T <- estimate_gcomp(muhat_T)
    est$gcomp_F <- estimate_gcomp(muhat_F)
  }
  if ("ipw" %in% estimators) {
    est$ipw_T <- estimate_ipw(Y, A, pihat_T)
    est$ipw_F <- estimate_ipw(Y, A, pihat_F)
  }
  if ("aipw" %in% estimators) {
    est$aipw_T <- estimate_aipw(Y, A, pihat_T, muhat_T)
    est$aipw_F <- estimate_aipw(Y, A, pihat_F, muhat_F)
  }
  if ("tmle" %in% estimators) {
    est$tmle_T <- estimate_tmle(Y, A, X, muhat_T, pihat_T)
    est$tmle_F <- estimate_tmle(Y, A, Z, muhat_F, pihat_F)
  }

  ## Compute bootstrap SEs for regression estimator
  if (("gcomp" %in% estimators) && bootstrap_reg) {
    datT <- bind_cols(as_tibble(X), A = A, Y = Y)
    datF <- bind_cols(as_tibble(Z), A = A, Y = Y)
    if (nonpar) {
      bs_T <- boot(datT, estimate_gcomp_nonpar, R=bootNum,
                   SL.library=sl_mu, folds = folds)
      bs_F <- boot(datF, estimate_gcomp_nonpar, R=bootNum,
                   SL.library=sl_mu,  folds = folds)
    }
    else {
      bs_T <- boot(datT, estimate_gcomp_par, R=bootNum, covs="X")
      bs_F <- boot(datF, estimate_gcomp_par, R=bootNum, covs="Z")
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
    est$gcomp_oracle <- estimate_gcomp(mu_true)
  }
  if ("ipw" %in% oracle_estimators) {
    est$ipw_oracle <- estimate_ipw(Y, A, pi_true)
  }
  if ("aipw" %in% oracle_estimators) {
    est$aipw_oracle <- estimate_aipw(Y, A, pi_true, mu_true)
  }
  if ("tmle" %in% oracle_estimators) {
    est$tmle_oracle <- estimate_tmle(Y, A, X, mu_true, pi_true)
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


#' Write results.
write_results <- function(results, counter, estimators, outdir, include="all",
                          suffix="") {
  if (!dir.exists(outdir)) {dir.create(outdir)}
  if (include[1] == "all") {
    include <- c("data", "predictions", "remainder", "metrics", "coefs")
  }
  # Write data
  if (("data" %in% include) && ("data" %in% names(results))) {
    subdir <- file.path(outdir, "data")
    if (!dir.exists(subdir)) {dir.create(subdir)}
    path <- file.path(subdir, paste0("iteration", counter, "_data", suffix, ".txt"))
    write_tsv(results[["data"]], path)
  }
  if (("predictions" %in% include) && ("predictions" %in% names(results))) {
    subdir <- file.path(outdir, "predictions")
    if (!dir.exists(subdir)) {dir.create(subdir)}
    path <- file.path(subdir, paste0("iteration", counter, "_predictions", suffix, ".txt"))
    write_tsv(results[["predictions"]], path)
  }
  if (("remainder" %in% include) && ("remainder" %in% names(results))) {
    subdir <- file.path(outdir, "remainder")
    if (!dir.exists(subdir)) {dir.create(subdir)}
    path = file.path(subdir, paste0("iteration", counter, "_remainder", suffix, ".txt"))
    write_tsv(results[["remainder"]], path)
  }
  if (("metrics" %in% include) && ("metrics" %in% names(results))) {
    subdir <- file.path(outdir, "metrics")
    if (!dir.exists(subdir)) {dir.create(subdir)}
    path = file.path(subdir, paste0("iteration", counter, "_",
                     paste(estimators, collapse = "_"), suffix, ".txt"))
    write_tsv(results[["metrics"]], path)
  }
  if (("coefs" %in% include) && ("coefs" %in% names(results))) {
    subdir <- file.path(outdir, "superlearner_coefs")
    if (!dir.exists(subdir)) {dir.create(subdir)}
    pi_path = file.path(subdir, paste0("iteration", counter, "_pihat", suffix, ".txt"))
    write_tsv(results[["coefs"]][["pihat_coefs"]], pi_path)
    mu_path = file.path(subdir, paste0("iteration", counter, "_muhat", suffix, ".txt"))
    write_tsv(results[["coefs"]][["muhat_coefs"]], mu_path)
  }
}
