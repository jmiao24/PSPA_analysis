rm(list = ls())
require(data.table)
source("Fun.R")
library(optparse)

options(stringsAsFactors = F)
option_list = list(
  make_option("--tissue", action = "store", default = NA, type = "character")
)
opt = parse_args(OptionParser(option_list=option_list))

tissue <- opt$tissue
# tissue <- "Brain_Cerebellar_Hemisphere"

# Input covariates
X_labeled_cov <- fread(paste0("./sexde/result/sexde/", tissue, ".lab_x.covs.txt"))
X_unlabeled_cov <- fread(paste0("./sexde/result/sexde/", tissue, ".unlab_x.covs.txt"))
# Variables for SVA
X_labeled_Sur_Y <- fread(paste0("./sexde/result/sexde/aux/covs.svs.", tissue, ".lab_y.v8.FINAL.txt"))
X_labeled_Sur_Yhat <- fread(paste0("./sexde/result/sexde/aux/covs.svs.", tissue, ".lab_yhat.v8.FINAL.txt"))
X_unlabeled_Sur <- fread(paste0("./sexde/result/sexde/aux/covs.svs.", tissue, ".unlab_yhat.v8.FINAL.txt"))

# Gene expression
Y_labeled <- fread(paste0("./sexde/result/hyfa/", tissue, ".lab_y.csv"))
Yhat_labeled <- fread(paste0("./sexde/result/hyfa/", tissue, ".lab_yhat.csv"))
Yhat_unlabeled <- fread(paste0("./sexde/result/hyfa/", tissue, ".unlab_yhat.csv"))

# Match the order
X_labeled_Sur_Y <- X_labeled_Sur_Y[match(X_labeled_cov$V1, X_labeled_Sur_Y$V1), ]
X_labeled_Sur_Yhat <- X_labeled_Sur_Yhat[match(X_labeled_cov$V1, X_labeled_Sur_Yhat$V1), ]
X_unlabeled_Sur <- X_unlabeled_Sur[match(X_unlabeled_cov$V1, X_unlabeled_Sur$V1), ]
keep_n <- min(ncol(X_labeled_Sur_Y), ncol(X_labeled_Sur_Yhat))

Y_labeled <- Y_labeled[match(X_labeled_cov$V1, Y_labeled$"Participant ID"), ]
Yhat_labeled <- Yhat_labeled[match(X_labeled_cov$V1, Yhat_labeled$"Participant ID"), ]
Yhat_unlabeled <- Yhat_unlabeled[match(X_unlabeled_cov$V1, Yhat_unlabeled$"Participant ID"), ]

X_labeled_Y <- cbind(X_labeled_cov[, -1], X_labeled_Sur_Y[, -1])
X_labeled_Yhat <- cbind(X_labeled_cov[, -1], X_labeled_Sur_Yhat[, -1])
X_unlabeled <-  cbind(X_unlabeled_cov[, -1], X_unlabeled_Sur[, -1])

Y_labeled <- as.data.frame(Y_labeled[, -c(1:2)])
Yhat_labeled <- as.data.frame(Yhat_labeled[, -c(1:2)])
Yhat_unlabeled <- as.data.frame(Yhat_unlabeled[, -c(1:2)])

Sex_labeled_Y <- X_labeled_Y[, 1]
Cov_labeled_Y <- X_labeled_Y[, -1]
Sex_labeled_Yhat <- X_labeled_Yhat[, 1]
Cov_labeled_Yhat <- X_labeled_Yhat[, -1]
Sex_unlabeled_Yhat <- X_unlabeled[, 1]
Cov_unlabeled_Yhat <- X_unlabeled[, -1]

# Get residual for Sex
Sex_labeled_Y_res <- residuals(lm(SEX ~ ., data= X_labeled_Y))
Sex_labeled_Yhat_res <- residuals(lm(SEX ~ ., data= X_labeled_Yhat))
Sex_unlabeled_Yhat_res <- residuals(lm(SEX ~ ., data= X_unlabeled))

# Get the residual for Sex
# Initialize Y_res with the same dimensions as Y
top <- ncol(Yhat_labeled)
Yhat_labeled_res <- matrix(nrow = nrow(Yhat_labeled), ncol = top)
for (i in 1:top) {
  # Linear regression of Y's column on X
  model <- lm(Yhat_labeled[, i] ~ ., data = Cov_labeled_Yhat[, c(1:6)])
  
  # Extract residuals and store in Y_res
  Yhat_labeled_res[, i] <- residuals(model)
}
fwrite(Yhat_labeled_res, paste0("./GE_residuals/", tissue, "_Yhat_labeled.txt.gz"), sep = "\t", quote = F, row.names = F, col.names = T)

# Get the residual for Sex
# Initialize Y_res with the same dimensions as Y
top <- ncol(Y_labeled)
Y_labeled_res <- matrix(nrow = nrow(Y_labeled), ncol = top)
for (i in 1:top) {
  # Linear regression of Y's column on X
  model <- lm(Y_labeled[, i] ~ ., data = Cov_labeled_Y[, c(1:6)])
  
  # Extract residuals and store in Y_res
  Y_labeled_res[, i] <- residuals(model)
}
fwrite(Yhat_labeled_res, paste0("./GE_residuals/", tissue, "_Y_labeled.txt.gz"), sep = "\t", quote = F, row.names = F, col.names = T)

# Get the residual for Sex
# Initialize Y_res with the same dimensions as Y
top <- ncol(Yhat_unlabeled)
Yhat_unlabeled_res <- matrix(nrow = nrow(Yhat_unlabeled), ncol = top)
for (i in 1:top) {
  # Linear regression of Y's column on X
  model <- lm(Yhat_unlabeled[, i] ~ ., data = Cov_unlabeled_Yhat[, c(1:6)])
  
  # Extract residuals and store in Y_res
  Yhat_unlabeled_res[, i] <- residuals(model)
}
fwrite(Yhat_labeled_res, paste0("./GE_residuals/", tissue, "_Yhat_unlabeled.txt.gz"), sep = "\t", quote = F, row.names = F, col.names = T)

# Compute column-wise correlation
column_correlations <- mapply(function(x, y) cor(x, y, use = "complete.obs"), as.data.frame(Yhat_labeled_res), as.data.frame(Y_labeled_res))
fwrite(data.frame(GENE = colnames(Yhat_labeled), cor = column_correlations), paste0("./corr/", tissue, ".txt.gz"), sep = "\t", quote = F, row.names = F, col.names = T)


num_cols <- ncol(Yhat_labeled_res)

# Pre-allocate vectors with the final size needed
beta_vec_classic <- numeric(num_cols)
beta_vec_PopInf <- numeric(num_cols)
beta_vec_PP <- numeric(num_cols)
beta_vec_eff <- numeric(num_cols)

se_vec_classic <- numeric(num_cols)
se_vec_PopInf <- numeric(num_cols)
se_vec_PP <- numeric(num_cols)
se_vec_eff <- numeric(num_cols)

P_vec_classic <- numeric(num_cols)
P_vec_PopInf <- numeric(num_cols)
P_vec_PP <- numeric(num_cols)
P_vec_eff <- numeric(num_cols)

# Iterate over columns
for (i in 1:num_cols) {
  # print(i)
  Y_labeled_tmp <- Y_labeled_res[, i]
  Yhat_labeled_tmp <- Yhat_labeled_res[, i]
  Yhat_unlabeled_tmp <- Yhat_unlabeled_res[, i]
  
  res_PP <- pp_ols_asymptotic(Sex_labeled_Y_res, Sex_unlabeled_Yhat_res, Y_labeled_tmp, Yhat_labeled_tmp, Yhat_unlabeled_tmp, alpha = 0.05, sandwich = F)
  res_eff <- eff_ols_asymptotic(Sex_labeled_Y_res, Sex_unlabeled_Yhat_res, Y_labeled_tmp, Yhat_labeled_tmp, Yhat_unlabeled_tmp, alpha = 0.05, sandwich = F)
  res_PopInf <- PopInf_ols_asymptotic(Sex_labeled_Y_res, Sex_unlabeled_Yhat_res, Y_labeled_tmp, Yhat_labeled_tmp, Yhat_unlabeled_tmp, alpha = 0.05, sandwich = F)
  res_classic <- classic_ols_asymptotic(Sex_labeled_Y_res, Y_labeled_tmp, alpha = 0.05, sandwich = F)
  
  beta_vec_classic[i] <- res_classic$theta[2]
  beta_vec_PopInf[i] <- res_PopInf$theta[2]
  beta_vec_PP[i] <- res_PP$theta[2]
  beta_vec_eff[i] <- res_eff$theta[2]
  
  se_vec_classic[i] <- res_classic$se[2]
  se_vec_PopInf[i] <- res_PopInf$se[2]
  se_vec_PP[i] <- res_PP$se[2]
  se_vec_eff[i] <- res_eff$se[2]
  
  P_vec_classic[i] <- res_classic$P[2]
  P_vec_PopInf[i] <- res_PopInf$P[2]
  P_vec_PP[i] <- res_PP$P[2]
  P_vec_eff[i] <- res_eff$P[2]
}
df_out <- data.frame(beta_vec_classic, beta_vec_PopInf, beta_vec_PP, beta_vec_eff, se_vec_classic, se_vec_PopInf, se_vec_PP, se_vec_eff, P_vec_classic, P_vec_PopInf, P_vec_PP, P_vec_eff)
fwrite(df_out, 
       paste0("./SexBSE/", tissue, ".txt.gz"),
       sep = "\t", quote = F, row.names = F, col.names = T)
