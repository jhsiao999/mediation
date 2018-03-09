#' @title linea model for genomewide expression data with gene-specific covariates
#'
#' @param Y gene by sample expression matrix (G by N). Required to be normalized and log transformed.
#' @param X sample condition labels, assumed to be binary for now.
#' @param M data for the mediating variable, a gene by sample matrix (G by N).
#'
#' @export

lm.varyingCovariates <- function(Y, X, M) {

  require(assertthat)
  assert_that(all.equal(dim(Y), dim(M)))
  assert_that(all.equal(length(unique(X)), 2))

  Y <- as.matrix(Y)
  M <- as.matrix(M)
  X <- factor(X)

  G <- dim(Y)[1]
  N <- dim(Y)[2]

  fits <- lapply(1:G, function(g) {
    M_g <- M[g,]
    Y_g <- Y[g,]
    design_g <- model.matrix(~X + M_g)
    fit <- lm.fit(y=Y_g, x=design_g)
    return(fit)
  })

  coefs <- do.call(rbind, lapply(fits, "coef"))
  rownames(coefs) <- rownames(Y)

  cov.coefficient <- lapply(fits, function(x) {
    chol2inv(x$qr$qr, size = x$qr$rank) })

  stdev.unscaled <- do.call(rbind, lapply(cov.coefficient, function(x) {
    sqrt(diag(x)) }))
  colnames(stdev.unscaled) <- colnames(coefs)

  sigma <- sapply(fits, function(x) {
    sqrt(sum(x$residuals^2)/x$df.residual)
  })

  predicted <- do.call(rbind, lapply(fits, function(x) {
    x$fitted.values
  }))

  return(list(coefs=coefs,
              stdev.unscaled=stdev.unscaled,
              sigma=sigma,
              predicted=predicted))
}




