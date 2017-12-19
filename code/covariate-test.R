# # <---- subset data to the tissue pair
# # compare human heart and kidney
# index_sample <- which(samples_sub$Species == "human"& (samples_sub$Tissue == "heart"|samples_sub$Tissue == "kidney"))
# exprs_pair <- exprs_subset[,index_sample]
# methyl_pair <- methyl_subset[,index_sample]
# species <- droplevels.factor(samples$Species[index_sample])
# tissue <- droplevels.factor(samples$Tissue[index_sample])
# RIN <- samples$RIN[index_sample]

df <- get(load(file="data/example-kidney-v-heart-human.rda"))

# Extract parameters
N <- ncol(df$exprs_pair)
ngenes <- nrow(df$methyl_pair)
tissue <- df$tissue

# <---- Model 1: exprs ~ tissue
# specify tissue coding
design_1 <- model.matrix(~tissue)
design_1[design_1[,2]==0,2] <- -1

model_1 <- lmFit(exprs_pair, design_1)

# <---- Model 3: exprs corrected for methylation ~ tissue
# specify design matrix
resid_exprs <- array(0, dim = dim(exprs_pair))
for (index in 1:nrow(exprs_pair)){
  resid_exprs[index,] <- lm(t(exprs_pair[index, ]) ~ t(methyl_pair[index, ]))$resid
}
rownames(resid_exprs) <- rownames(exprs_pair)

model_3 <- lmFit(resid_exprs, design_1)


# get effect sizes
beta1 <- coef(model_1[,2])
beta3 <- coef(model_3[,2])

se_beta1 <- se_beta2 <- se_beta3 <- cov_beta13 <- vector("numeric", ngenes)

# <---- get variances of beta1^S tissue effect
se_beta1 <- model_1$sigma*sqrt((solve(t(design_1)%*%design_1))[2,2])

head(cbind(se_beta1,model_1$sigma*sqrt(model_1$cov.coefficients[2,2])))

# <---- get variances of beta2 methylation effect
for (g in 1:length(se_beta2)) {
  design_2g <- model.matrix(~ as.numeric(methyl_pair[g,]))
  sigma_g <- model_1$sigma[g]
  se_beta2[index] <- sigma_g*sqrt((solve(t(design_2g)%*%design_2g))[2,2])
}

# <---- get variances of beta3 tissue effect
A <- solve(t(design_1)%*%design_1)%*%t(design_1)
contr.vector <- array(c(0,1), dim=c(2,1))

# compute beta3 by hand
for (g in 1:length(se_beta3)) {
  M_g <- t(methyl_pair[g,])
  design_2g <- model.matrix(~ as.numeric(M_g))
  A_2g <- solve(t(design_2g)%*%design_2g)%*%t(design_2g)
  sigma_g <- model_1$sigma[g]
  var_beta2g <- (se_beta2^2)[g]
  var_part1 <- (se_beta1[g])^2
  var_part2 <- ( A%*%M_g%*%var_beta2g%*%t(M_g)%*%t(A) )[2,2]
  var_part3 <- ( 2*(sigma_g^2)*A%*%t(A_2g)%*%contr.vector%*%t(M_g)%*%t(A) )[2,2]
  se_beta3[g] <- sqrt(var_part1 + var_part2 + var_part3)
}

# cov(beta1,beta3)
for (g in 1:length(cov_beta13)) {
  M_g <- t(methyl_pair[g,])
  design_2g <- model.matrix(~ as.numeric(M_g))
  A_2g <- solve(t(design_2g)%*%design_2g)%*%t(design_2g)
  sigma_g <- model_1$sigma[g]
  var_part1 <- (se_beta1[g])^2
  cov_beta13[g] <- var_part1 - (sigma_g^2)*((A %*% t(A_2g) %*%contr.vector%*%t(M_g)%*%t(A))[2,2])
}

# cov(beta1-beta3)
cov_diff_sqrt <- sqrt(se_beta1^2 + se_beta3^2 - 2*cov_beta13)
beta_diff <- beta1-beta3




#' #' @param beta1 effect size from model1 (expressin ~ species)
#' #' @param beta2 effect size from model2 (expression|methylation ~ species)
#' #' @param M methylation matrix gene by sample
#' #' @param RIN RIN scores
#' #' @param n1 sample size for human
#' #' @param n2 sample size for chimp
#' #' @param ngenes number of genes
#'
#' library(data.table)
#' df1 <- get(load("Limma_output_fit1.rda"))
#' df2 <- get(load("Limma_output_fit2.rda"))
#'
#' beta1 <- coef(df1[,2])
#' beta2 <- coef(df2[,2])
#' N <- n1+n2
#' RIN <- df1$design[,3]
#' ngenes <- dim(coef(df1))[1]
#'
#' sigma2 <- df1$sigma^2
#' design <- cbind(int=rep(1,N),species=c(rep(1,n1), rep(-1,n2)), RIN=RIN)
#'
#' # compute some matrices to simply calculations
#' H <- design%*%solve(t(design)%*%design)%*%t(design)
#' I <- diag(N)
#' A <- solve(t(design)%*%design)%*%t(design)
#'
#' var_beta1 <- var_beta2 <- cov_beta12 <- vector("numeric", ngenes)
#'
#' for (g in 1:ngenes) {
#'   M_g <- M[g,]
#'   D <-  1/(2*n*t(M_g)%*%(I-H)%*%M_g)
#'
#'   # standard error in model 1
#'   var_beta1[g] <- sigma2*t(solve(t(design)%*%design))
#'
#'   # standard error of methylation
#'   var_M <- sigma2/((t(M_g)%*%(I-H)%*%M_g)^2)*t(M_g)%*%(I-H)%*%t(I-H)%*%M_g
#'
#'   # standard error in model 2
#'   var_beta2[g] <- var_beta1[g] + A%*%M_g*var_M*t(M_g)t(A)-2*sigma2*A%*%t(I-H)%*%M_g%*%t(M_g)%*%t(A)
#'
#'
#'   # covariance of effect sizes
#'   cov_beta12[g] <- var_beta1[g] - sigma2/(t(M_g)%*%(I-H)%*%M_g)*(A%*%t(I-H)%*%M_g%*%t(M_g)%*%t(A))
#' }
#'
#' cov_diff <- var_beta1 + var_beta2 - 2*cov_beta12
#'
#'
#'
#' library(ashr)
#' ash(beta12, cof_diff)
#'
#'
#'









