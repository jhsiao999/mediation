#' @title Simulate data for mediational analysis
#'
#' @param N_each number of samples each group
#' @param Y_means, Y_sds dependent variable means and standard deviations for each group
#' @param M_means, M_sds mediating variable means and standard deviations for each group
#'
#' @export
#'
#' @examples
#' df.sim <- mediate.simulate(N_each=5, G=50, Y_means = c(0,0), Y_sd=c(1,1), M_means = c(0,1), M_sds = c(1,1))
#' test.sobel <- mediate.test(Y=df.sim$Y, X=df.sim$X, M=df.sim$M)

mediate.simulate <- function(N_each, G, Y_means, Y_sds,
                           M_means, M_sds) {
  #  n_each <- 40
  condition <- factor(rep(c(0,1), each = N_each), levels=c(0,1))

  set.seed(17)
  if (!is.null(G) & G > 1) {
    df.list <- lapply(1:G, function(g) {
      yy <- c(rnorm(N_each, Y_means[1], Y_sds[1]),
            rnorm(N_each, Y_means[2], Y_sds[2]))
      mm <- c(rnorm(N_each, M_means[1], M_sds[1]),
              rnorm(N_each, M_means[2], M_sds[2]))
      return(list(yy=yy, mm=mm))
    })
  }
  Y <- do.call(rbind, lapply(df.list, "[[", 1))
  M <- do.call(rbind, lapply(df.list, "[[", 2))

  return(list(Y=Y, X=condition, M=M))
}


