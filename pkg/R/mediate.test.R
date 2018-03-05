#' @title Mediation test
#'
#' @description Compute statistics for evaluating mediating effect in genomewide expression data.
#'
#' @param Y gene by sample expression matrix (G by N). Required to be normalized and log transformed.
#' @param X sample condition labels, assumed to be binary for now.
#' @param M data for the mediating variable, a gene by sample matrix (G by N).
#'
#' @export
#'
mediate.test <- function(Y, X, M) {

  library(limma)
  library(assertthat)

  # evaluate argument conditions
  assert_that(all.equal(dim(Y), dim(M)))
  assert_that(all.equal(length(unique(X)), 2))
  assert_that(all.equal(length(X), dim(Y)[2]))

  G <- dim(Y)[1]
  N <- dim(Y)[2]

  Y <- as.matrix(Y)
  M <- as.matrix(M)
  X <- factor(X)

  # Fitting model 1: Y_g = \gamma_{1g} + \tau_g X + \epsilon_{1g}
  design_1 <- model.matrix(~X)
  model_1 <- lmFit(Y, design_1)

  # Fitting model 2: M_g = \gamma_{2g} + \alpha X + \epsilon_{2g}
  design_2 <- model.matrix(~X)
  model_2 <- lmFit(M, design_2)

  # Fitting model 3: M_g = \gamma_{3g} + \tau^{\prime}_g X + \beta_g M_g \epsilon_{3g}
  model_3 <- lm.varyingCovariates(Y=Y, X=X, M=M)

  tau <- coef(model_1)[,2]
  tau_prime <- model_3$coefs[,2]
  alpha <- coef(model_2)[,2]
  beta <- model_3$coefs[,3]
  sigma_alpha <- model_2$stdev.unscaled[,2]*model_2$sigma
  sigma_beta <- model_3$stdev.unscaled[,3]*model_3$sigma
  sigma_tau <- model_1$stdev.unscaled[,2]*model_1$sigma
  sigma_tau_prime <- model_3$stdev.unscaled[,2]*model_3$sigma
  phi_x.m <- sapply(1:nrow(M), function(i) cor(as.integer(X), unlist(M[i,])))

  d <- tau-tau_prime
  se.sobel <- sqrt((alpha*sigma_beta)^2 + (beta*sigma_alpha)^2)
  se.fs <- sqrt(sigma_tau^2 + sigma_tau_prime^2 - 2*sigma_tau*sigma_tau_prime*sqrt(1-phi_x.m^2))

  return(data.frame(d=d,
              se.sobel=se.sobel,
              tau=tau,
              tau_prime=tau_prime,
              ab =alpha*beta,
              alpha=alpha,
              beta=beta,
              sigma_alpha=sigma_alpha,
              sigma_beta=sigma_beta,
              sigma_tau=sigma_tau,
              sigma_tau_prime=sigma_tau_prime,
              phi_x.m = phi_x.m
              ))
}


