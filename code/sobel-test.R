#' @title Sobel's test
#'
#' @description Compute Sobel's test statistic for the case of fitting two covariates, where
#'   one covariate is constant across genes and the other covariate is gene-specific
#'
#' @param exprs gene by sample expression matrix (G by N)
#' @param fixed_covariate length-N vector of sample condition labels (numeric vector,
#'   such as 0, 1 or 1,2).
#' @param varing_covariate gene by sample covariate matrix (G by N)
#'
#' @return a list of estimates
#'   d: estimated indirect effect
#' @export
#'
# exprs=exprs_pair
# fixed_covariates=list(tissue=tissue)
# varying_covariate=methyl_pair

sobel <- function(exprs, fixed_covariates=list(), varying_covariate) {

  source("../code/lm-varying-covariate.R")
  library(limma)

  # Fitting model 1
  design_1 <- model.matrix(~fixed_covariates[[1]])
  colnames(design_1)[-1] <- names(fixed_covariates)
  model_1 <- lmFit(exprs, design_1)

  # Fitting model 2
  design_2 <- model.matrix(~fixed_covariates[[1]])
  colnames(design_2)[-1] <- names(fixed_covariates)[1]
  model_2 <- lmFit(varying_covariate, design_2)

  # Fitting model 3
  model_3 <- lm_varying_covariate(exprs,
                                  fixed_covariates,
                                  varying_covariate)

  tau <- coef(model_1)[,2]
  tau_prime <- model_3$coefs[,2]
  alpha <- coef(model_2)[,2]
  beta <- model_3$coefs[,3]
  sigma2_alpha <- (model_2$stdev.unscaled[,2]*model_2$sigma)^2
  sigma2_beta <- (model_3$stdev.unscaled[,3]*model_3$sigma)^2

  d <- tau-tau_prime
  se <- sqrt((alpha^2)*sigma2_beta + (beta^2)*sigma2_alpha)
  return(list(tau=tau,
              tau_prime=tau_prime,
              alpha=alpha,
              beta=beta,
              sigma2_alpha=sigma2_alpha,
              sigma2_beta=sigma2_beta,
              d=d,
              se=se))
}


# sobel test applying to model 1 data after voom transformation
# exprs=exprs_pair
# fixed_covariates=list(tissue=tissue)
# varying_covariate=methyl_pair
sobel_voom_model1 <- function(exprs, fixed_covariates=list(), varying_covariate) {

  source("../code/lm-varying-covariate.R")
  library(limma)

  # Fitting model 1
  design_1 <- model.matrix(~fixed_covariates[[1]])
  colnames(design_1)[-1] <- names(fixed_covariates)

  exprs_counts <- 2^exprs
  exprs_voom <- voom(exprs_counts, design=design_1, normalize.method = "none")

  model_1 <- lmFit(exprs_voom, design_1)
  model_1 <- eBayes(model_1)

  # Fitting model 2
  design_2 <- model.matrix(~fixed_covariates[[1]])
  colnames(design_2)[-1] <- names(fixed_covariates)[1]
  model_2 <- lmFit(varying_covariate, design_2)

  # Fitting model 3
  model_3 <- lm_varying_covariate(exprs,
                                  fixed_covariates,
                                  varying_covariate)

  tau <- coef(model_1)[,2]
  tau_prime <- model_3$coefs[,2]
  alpha <- coef(model_2)[,2]
  beta <- model_3$coefs[,3]
  sigma2_alpha <- (model_2$stdev.unscaled[,2]*model_2$sigma)^2
  sigma2_beta <- (model_3$stdev.unscaled[,3]*model_3$sigma)^2

  d <- tau-tau_prime
  se <- sqrt((alpha^2)*sigma2_beta + (beta^2)*sigma2_alpha)
  return(list(tau=tau,
              tau_prime=tau_prime,
              alpha=alpha,
              beta=beta,
              sigma2_alpha=sigma2_alpha,
              sigma2_beta=sigma2_beta,
              d=d,
              se=se))
}




# sobel test applying to model 1 and 3 data after voom transformation
# exprs=df$exprs_pair
# fixed_covariates=list(tissue=df$tissue)
# varying_covariate=df$methyl_pair
sobel_voom_model13 <- function(exprs, fixed_covariates=list(), varying_covariate) {

  source("code/lm-varying-covariate.R")
  library(limma)

  # Fitting model 1
  design_1 <- model.matrix(~fixed_covariates[[1]])
  colnames(design_1)[-1] <- names(fixed_covariates)

  exprs_counts <- 2^exprs
  exprs_voom <- voom(exprs_counts, design=design_1, normalize.method = "none")

  model_1 <- lmFit(exprs_voom, design_1)
  model_1 <- eBayes(model_1)

  # Fitting model 2
  design_2 <- model.matrix(~fixed_covariates[[1]])
  colnames(design_2)[-1] <- names(fixed_covariates)[1]
  model_2 <- lmFit(varying_covariate, design_2)

  # Fitting model 3
  model_3 <- lm_varying_covariate(exprs,
                                  fixed_covariates,
                                  varying_covariate)

  tau <- coef(model_1)[,2]
  tau_prime <- model_3$coefs[,2]
  alpha <- coef(model_2)[,2]
  beta <- model_3$coefs[,3]
  sigma2_alpha <- (model_2$stdev.unscaled[,2]*model_2$sigma)^2
  sigma2_beta <- (model_3$stdev.unscaled[,3]*model_3$sigma)^2

  d <- tau-tau_prime
  se <- sqrt((alpha^2)*sigma2_beta + (beta^2)*sigma2_alpha)
  return(list(tau=tau,
              tau_prime=tau_prime,
              alpha=alpha,
              beta=beta,
              sigma2_alpha=sigma2_alpha,
              sigma2_beta=sigma2_beta,
              d=d,
              se=se))
}

# version with two covariates
#
# exprs=exprs_pair
# fixed_covariates=list(tissue=tissue,RIN=RIN)
# varying_covariate=methyl_pair
#
sobel_2 <- function(exprs, fixed_covariates=list(), varying_covariate) {

  source("../code/lm-varying-covariate.R")
  library(limma)

  # Fitting model 1
  design_1 <- model.matrix(~fixed_covariates[[1]] + fixed_covariates[[2]])
  colnames(design_1)[-1] <- names(fixed_covariates)
  model_1 <- lmFit(exprs, design_1)

  # Fitting model 2
  design_2 <- model.matrix(~fixed_covariates[[1]])
  colnames(design_2)[-1] <- names(fixed_covariates)[1]
  model_2 <- lmFit(varying_covariate, design_2)

  # Fitting model 3
  model_3 <- lm_varying_covariate_2(exprs,
                                  fixed_covariates,
                                  varying_covariate)

  tau <- coef(model_1)[,2]
  tau_prime <- model_3$coefs[,2]
  alpha <- coef(model_2)[,2]
  beta <- model_3$coefs[,4]
  sigma2_alpha <- (model_2$stdev.unscaled[,2]*model_2$sigma)^2
  sigma2_beta <- (model_3$stdev.unscaled[,4]*model_3$sigma)^2

  d <- tau-tau_prime
  se <- sqrt((alpha^2)*sigma2_beta + (beta^2)*sigma2_alpha)
  return(list(tau=tau,
              tau_prime=tau_prime,
              alpha=alpha,
              beta=beta,
              sigma2_alpha=sigma2_alpha,
              sigma2_beta=sigma2_beta,
              d=d,
              se=se))
}
