# # <---- subset data to the pheno pair
# import sample labels
samples <- read.table("data/metadata.txt", sep = ",", header = TRUE)

# expression
exprs <- read.table("data/expression_vals.txt", sep=",", header = TRUE)

# hi-c data
hic <- read.table("data/hic_vals.txt", sep=",", header = TRUE)


# Extract parameters
exprs_pair <- exprs
methyl_pair <- hic
pheno <- samples$SP
N <- ncol(exprs_pair)
ngenes <- nrow(exprs_pair)

# <---- Model 1: exprs ~ pheno
# specify pheno coding
design_1 <- model.matrix(~pheno)

exprs_counts <- 2^exprs_pair
exprs_voom <- voom(exprs_counts, design=design_1, normalize.method = "none")

model_1 <- lmFit(exprs_voom, design_1)
model_1 <- eBayes(model_1)

# <---- Model 3: exprs corrected for methylation ~ pheno
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

# <---- get variances of beta1^S pheno effect
se_beta1 <- model_1$sigma*sqrt((solve(t(design_1)%*%design_1))[2,2])

# checking the validity of the results
head(cbind(se_beta1,model_1$sigma*sqrt(model_1$cov.coefficients[2,2])))

# <---- get variances of beta2 methylation effect
for (g in 1:length(se_beta2)) {
  design_2g <- model.matrix(~ as.numeric(methyl_pair[g,]))
  sigma_g <- model_1$sigma[g]
  se_beta2[index] <- sigma_g*sqrt((solve(t(design_2g)%*%design_2g))[2,2])
}

# <---- get variances of beta3 pheno effect
A <- solve(t(design_1)%*%design_1)%*%t(design_1)
contr.vector <- array(c(0,1), dim=c(2,1))

# compute beta3 by hand
# can't directly apply limma because limma requires
# that sample covariates have the same values across genes
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

cov_diff_sqrt <- sqrt(se_beta1^2 + se_beta3^2 - 2*cov_beta13)
beta_diff <- beta1-beta3

# putting results in a list
# beta_diff: differences between species effect on expression
# and species effect on expression data after regressing out hi-C
df_reg <- list(beta_diff=beta_diff, cov_diff_sqrt=cov_diff_sqrt)

# applying ash
# df: degrees of freedom = total number of samples - 2
ash_reg <- ash(as.vector(beta_diff), cov_diff_sqrt, df= length(pheno)-2)

# save(df_reg, ash_reg,
#      reg3, reg1,
#      file="output/sobeltest.Rmd/results-reg.rda")



par(mfrow=c(2,2), mar=c(4,4,4,1))
plot(x=df_reg$beta_diff,y=df_reg$cov_diff_sqrt, pch=16, cex=.6,
     xlab="Effect size", ylab="Standard error",
     main = "Effect size vs SE",
     col=ifelse(ash_reg$result$svalue<.05, "red", "black"))
abline(v=0, lty=2, col="grey20")
plot(x=df_reg$beta_diff, y=ash_reg$result$PosteriorMean, pch=16, cex=.6,
     xlab="Effect size", ylab="Posterior mean",
     main = "Effect size vs Posterior Mean",
     col=ifelse(ash_reg$result$svalue<.05, "red", "black"))
plot(x=df_reg$cov_diff_sqrt, y=ash_reg$result$PosteriorSD, pch=16, cex=.6,
     xlab="Standard error", ylab="Posterior standard deviation",
     main = "SE vs Posterior SD",
     col=ifelse(ash_reg$result$svalue<.05, "red", "black"))
plot(x=df_reg$beta_diff, y=-log10(ash_reg$result$svalue), pch=16, cex=.6,
     xlab="Effect size", ylab="-log10(s-value)",
     main = "Effect size vs -log10(s-value)",
     col=ifelse(ash_reg$result$svalue<.05, "red", "black"))






