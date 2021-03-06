% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_analysis.R
\name{estimate_pi_nonpar}
\alias{estimate_pi_nonpar}
\title{Estimate propensity scores nonparametrically using SuperLearner.}
\usage{
estimate_pi_nonpar(A, X, sl_lib, nsplits = NULL, split_inds = NULL, ...)
}
\arguments{
\item{A}{A vector of treatment values, in {0, 1}.}

\item{X}{A matrix or data frame of covariate values.}

\item{sl_lib}{A character vector of learner functions created via
SuperLearner, as in sl_lib <- SuperLearner::create.Learner(...)$names, where
... are the learner type and learner parameters. These are passed to
SuperLearner() to estimate the propensity score function.}

\item{nsplits}{The number of splits to use for sample splitting, if fold
indices aren't provided to the split_inds argument.}

\item{split_inds}{Fold indices for sample splitting, if nsplits is NULL. One
of split_inds or nsplits must be provided. split_inds takes precedence if
it's not NULL.}

\item{...}{Additional arguments passed to SuperLearner.}
}
\value{
A vector of estimated treatment propensities P(A = 1|X) in [0, 1].
}
\description{
Estimate propensity scores nonparametrically using SuperLearner.
}
