#' @title Mediation test for data after voom transformation
#'
#' @description Compute statistics for evaluating mediating effect in genomewide expression data. Apply voom transform to log2 expression data.
#'
#' @param Y gene by sample expression matrix (G by N). Required to be normalized and log transformed.
#' @param X sample condition labels, assumed to be binary for now.
#' @param M data for the mediating variable, a gene by sample matrix (G by N).
#'
#' @export
#'
mediate.test.voom <- function(Y, X, M) {

  library(limma)
  library(assertthat)
  library(ltm)

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
  # First transform the log2 expression to counts and
  # compute voom weights
  design_1 <- model.matrix(~X)
  wts <- vooma(Y, design = design_1)
  model_1 <- lmFit(Y, design_1, weights= wts$weights)
  model_1 <- eBayes(model_1)

  # Fitting model 2: M_g = \gamma_{2g} + \alpha X + \epsilon_{2g}
  design_2 <- model.matrix(~X)
  model_2 <- lmFit(M, design_2)
  model_2 <- eBayes(model_2)

  # Fitting model 3: M_g = \gamma_{3g} + \tau^{\prime}_g X + \beta_g M_g \epsilon_{3g}
  model_3 <- lm.varyingCovariates.voom(Y=Y, X=X, M=M, weights = NULL)

  tau <- coef(model_1)[,2]
  tau_prime <- model_3$coefs[,2]
  alpha <- coef(model_2)[,2]
  beta <- model_3$coefs[,3]
  sigma_alpha <- model_2$stdev.unscaled[,2]*sqrt(model_2$s2.post)
  sigma_beta <- model_3$stdev.unscaled[,3]*sqrt(model_3$s2.post)
  sigma_tau <- model_1$stdev.unscaled[,2]*sqrt(model_1$s2.post)
  sigma_tau_prime <- model_3$stdev.unscaled[,2]*sqrt(model_3$s2.post)

  corr.xm <- sapply(1:nrow(Y), function(g) {
    biserial.cor(M[g,], as.numeric(X))
  })

  # library(CorShrink)
  # corr.xm.shrink <- CorShrinkVector(corr.xm, rep(ncol(Y),nrow(Y)))

  se.sobel <- sqrt((alpha*sigma_beta)^2 + (beta*sigma_alpha)^2)
  se.fs <- sqrt(sigma_tau^2 + sigma_tau_prime^2 - 2*sigma_tau*sigma_tau_prime*sqrt(1-corr.xm^2))
  # se.fs.shrink <- sqrt(sigma_tau^2 + sigma_tau_prime^2 - 2*sigma_tau*sigma_tau_prime*sqrt(1-corr.xm.shrink^2))

  se.unbiased <- sqrt((alpha*sigma_beta)^2 + (beta*sigma_alpha)^2 - (sigma_alpha*sigma_beta)^2)
 # lik.fun <- function(x, mu, sigma) {
 #  sum(-((x-mu)^2)/(2*sigma^2) - log(sigma) - 0.5*log(2*pi))
 # }
 #
 # mu.model_1 <- model_1$coefficients%*%t(model.matrix(~X))
 # lik.model_1 <- sapply(1:nrow(Y), function(g) {
 #   lik.fun(Y[g,], mu=mu.model_1[g,], model_1$s2.post[g]) })
 #
 # mu.model_3 <- do.call(rbind, lapply(1:nrow(Y), function(g) {
 #   model_3$coefs[g,]%*%t(model.matrix(~X+M[g,])) }) )
 # lik.model_3 <- sapply(1:nrow(Y.de), function(g) {
 #   lik.fun(Y.de[g,], mu=mu.model_3[g,], model_3$s2.post[g]) })
 #
 #  lrt <- 2*(lik.model_3-lik.model_1)
 #  lrt.pval <- 1-pchisq(lrt, df = 1, lower.tail = F)

  return(data.frame(d=tau-tau_prime,
              se.sobel=se.sobel,
              se.fs=se.fs,
              # se.fs.shrink=se.fs.shrink,
              se.unbiased=se.unbiased,
              se.ab.z=sqrt((alpha/sigma_alpha)^2*(beta/sigma_beta)^2+1),
              tau=tau,
              tau_prime=tau_prime,
              ab=alpha*beta,
              ab.z=(alpha/sigma_alpha)*(beta/sigma_beta),
              alpha=alpha,
              beta=beta,
              sigma_alpha=sigma_alpha,
              sigma_beta=sigma_beta,
              sigma_tau = sigma_tau,
              sigma_tau_prime=sigma_tau_prime))
              # corr.xm=corr.xm,
              # corr.xm.shrink=corr.xm.shrink,
              #lrt.pval=lrt.pval,
              # cor.tau.tau_prime=cor.tau.tau_prime,
              # cor.tau.tau_prime.shrink=cor.tau.tau_prime.shrink,
#              df=min(model_1$df.total, model_3$df.total)[1]))
}


