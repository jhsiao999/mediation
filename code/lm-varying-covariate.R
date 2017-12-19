#' @title linear model for gene expression datasets with gene-specific covariates
#'
#' @param exprs sample by gene expressin matrix.
#' @param fixed_covariate length-N covariate vector constant across genes (such as phenotype).
#' @param varying_covariate sample by gene covariate measurement matrix.
#'
#' @export

lm_varying_covariate <- function(exprs, fixed_covariates=list(), varying_covariate) {

  assertthat::assert_that(all.equal(dim(exprs), dim(varying_covariate)))

  G <- dim(exprs)[1]

  # assume the number of fixed covariate is 1
  est <- lapply(1:G, function(g) {
    cov <- unlist(varying_covariate[g,])
    y <- unlist(exprs[g,])
#    design <- model.matrix(~fixed_covariates[[1]]+fixed_covariates[[2]] + cov)
#    colnames(design)[c(2,3)] <- names(fixed_covariates)
    design <- model.matrix(~fixed_covariates[[1]] + cov)
    colnames(design)[2] <- names(fixed_covariates)

    fit <- lm.fit(y=y, x=design)
    return(fit)
  })

  coefs <- do.call(rbind, lapply(est, "coef"))
  rownames(coefs) <- rownames(exprs)

  cov.coefficient <- lapply(est, function(x) {
    chol2inv(x$qr$qr, size = x$qr$rank) })

  stdev.unscaled <- do.call(rbind, lapply(cov.coefficient, function(x) {
    sqrt(diag(x)) }))
  colnames(stdev.unscaled) <- colnames(coefs)

  sigma <- sapply(est, function(x) {
    sqrt(sum(x$residuals^2)/x$df.residual)
  })
  return(list(coefs=coefs,
              stdev.unscaled=stdev.unscaled,
              sigma=sigma))
}




#
# version with two covariates
#
lm_varying_covariate_2 <- function(exprs, fixed_covariates=list(), varying_covariate) {

  assertthat::assert_that(all.equal(dim(exprs), dim(varying_covariate)))

  G <- dim(exprs)[1]

  # assume the number of fixed covariate is 1
  est <- lapply(1:G, function(g) {
    cov <- unlist(varying_covariate[g,])
    y <- unlist(exprs[g,])
    design <- model.matrix(~fixed_covariates[[1]]+fixed_covariates[[2]] + cov)
    colnames(design)[c(2,3)] <- names(fixed_covariates)

    fit <- lm.fit(y=y, x=design)
    return(fit)
  })

  coefs <- do.call(rbind, lapply(est, "coef"))
  rownames(coefs) <- rownames(exprs)

  cov.coefficient <- lapply(est, function(x) {
    chol2inv(x$qr$qr, size = x$qr$rank) })

  stdev.unscaled <- do.call(rbind, lapply(cov.coefficient, function(x) {
    sqrt(diag(x)) }))
  colnames(stdev.unscaled) <- colnames(coefs)

  sigma <- sapply(est, function(x) {
    sqrt(sum(x$residuals^2)/x$df.residual)
  })
  return(list(coefs=coefs,
              stdev.unscaled=stdev.unscaled,
              sigma=sigma))
}



