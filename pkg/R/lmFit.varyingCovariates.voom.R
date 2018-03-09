#' @title linea model for genomewide expression data with gene-specific covariates including voom weights
#'
#' @param Y gene by sample expression matrix (G by N). Required to be normalized and log transformed.
#' @param X sample condition labels, assumed to be binary for now.
#' @param M data for the mediating variable, a gene by sample matrix (G by N).
#'
#' @export
lm.varyingCovariates.voom <- function(Y, X, M, weights) {

  require(assertthat)
  assert_that(all.equal(dim(Y), dim(M)))
  assert_that(all.equal(length(unique(X)), 2))

  Y <- as.matrix(Y)
  M <- as.matrix(M)
  X <- factor(X)

  G <- dim(Y)[1]
  N <- dim(Y)[2]

  ncoef <- ncol(model.matrix(~X))+1
  beta <- matrix(0, nrow=G, ncol=ncoef)
  stdev.unscaled <- matrix(0, nrow=G, ncol=ncoef)
  df.residual <- vector("numeric", G)
  sigma <- vector("numeric", G)

  if (is.null(weights)) {
    # compute voom weights
    # adapted from limma::voom
    fits <- lm.varyingCovariates(Y=Y, X=X, M=M)

    Ameans <- fits$coefs[,1]
    sigma <- fits$sigma
    lib.size <- colSums(2^Y)
    sx <- Ameans + mean(log2(lib.size + 1)) - log2(1e+06)
    sy <- sqrt(sigma)
    l <- lowess(sx, sy, f = .5)
    f <- approxfun(l, rule = 2)
    fitted.values <- fits$predicted
    fitted.cpm <- 2^fitted.values
    fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + 1))
    fitted.logcount <- log2(fitted.count)
    weights <- 1/f(fitted.logcount)^4
    dim(weights) <- dim(fitted.logcount)
  }

  # adapted from lm.series
  for (g in 1:G) {
    M_g <- M[g,]
    Y_g <- as.vector(Y[g,])
    design_g <- model.matrix(~X + M_g)

    W_g <- as.vector(weights[g, ])
    out <- lm.wfit(design_g, Y_g, W_g)
    est <- !is.na(out$coef)
    beta[g, est] <- out$coef
    stdev.unscaled[g, est] <- sqrt(diag(chol2inv(out$qr$qr,
                                                 size = out$rank)))
    df.residual[g] <- out$df.residual
    sigma[g] <- sqrt(mean(out$effects[-(1:out$rank)]^2))
    QR_g <- qr(design_g)
    cov.coef_g <- chol2inv(QR_g$qr, size = QR_g$rank)
    est_g <- QR_g$pivot[1:QR_g$rank]
  }
  #dimnames(cov.coef) <- list(coef.names[est], coef.names[est])

  return(list(coefs=beta,
              stdev.unscaled=stdev.unscaled,
              sigma=sigma))
}




