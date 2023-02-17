# Script that performs calculations for comparing performance of updated prognostic models with original models
source("nb_diff.R")
source("dca.R")
library(pROC)

######
# Comparison of areas under the ROC curves
# MELD
roc_meld_new <- roc(test.data$D90_surv, test.data$meld.surv.updated)
roc_meld_old <- roc(test.data$D90_surv, test.data$MELD.surv)
roc.test(roc_meld_old, roc_meld_new) 

# CLIF-C ACLF
roc_clif_new <- roc(test.data$D90_surv, test.data$clif.surv.updated)
roc_clif_old <- roc(test.data$D90_surv, test.data$CLIF.surv)
roc.test(roc_clif_new, roc_clif_old) 

# Lille
roc_lille_new <- roc(test.data.c$D90_surv, test.data.c$lille.surv.updated)
roc_lille_old <- roc(test.data.c$D90_surv, test.data.c$Lille.surv)
roc.test(roc_lille_new, roc_lille_old)

######
# Comparison of NB
# First define the variables corresponding to mortality
test.data.c$clif.mort <- 1 - test.data.c$CLIF.surv
test.data.c$meld.mort <- 1 - test.data.c$MELD.surv
test.data.c$meld.mort2 <- 1 - test.data.c$MELD.surv2
test.data.c$lille.mort <- 1 - test.data.c$Lille.surv

test.data.c$clif.mort.new <- 1 - test.data.c$clif.surv.updated
test.data.c$meld.mort.new <- 1 - test.data.c$meld.surv.updated
test.data.c$lille.mort.new <- 1 - test.data.c$lille.surv.updated

# Initiate the bootstrap procedure and calculate p-values
library(boot)
set.seed(34)
R <- 500 # Number of bootstrap iterations

# Calculate p-values for comparison of MELD scores
boot.diff.meld <- boot(data=test.data.c, statistic = nb_diff,
                        R = R, outcome = "D90_DTH", pred1 = "meld.mort",
                        pred2 = "meld.mort.new", xstart = 0.25, xstop = 0.75,
                        step = 0.05)
pvalue.meld <- NULL

for(i in 1:length(boot.diff.meld$t0)){
  pvalue.meld <- c(pvalue.meld, mean(abs(boot.diff.meld$t[,i] - boot.diff.meld$t0[i]) > abs(boot.diff.meld$t0[i])))
}

# Calculate p-values for comparison of MELD2 scores
boot.diff.meld2 <- boot(data=test.data.c, statistic = nb_diff,
                       R = R, outcome = "D90_DTH", pred1 = "meld.mort2",
                       pred2 = "meld.mort.new", xstart = 0.25, xstop = 0.75,
                       step = 0.05)
pvalue.meld2 <- NULL

for(i in 1:length(boot.diff.meld2$t0)){
  pvalue.meld2 <- c(pvalue.meld2, mean(abs(boot.diff.meld2$t[,i] - boot.diff.meld2$t0[i]) > abs(boot.diff.meld2$t0[i])))
}

# Calculate p-values for comparison of CLIF-scores
boot.diff.clif <- boot(data=test.data.c, statistic = nb_diff,
                       R = R, outcome = "D90_DTH", pred1 = "clif.mort",
                       pred2 = "clif.mort.new", xstart = 0.25, xstop = 0.75,
                       step = 0.05)
pvalue.clif <- NULL

for(i in 1:length(boot.diff.clif$t0)){
  pvalue.clif <- c(pvalue.clif, mean(abs(boot.diff.clif$t[,i] - boot.diff.clif$t0[i]) > abs(boot.diff.clif$t0[i])))
}

# Calculate p-values for comparison of Lille scores
boot.diff.lille <- boot(data=test.data.c, statistic = nb_diff,
                       R = R, outcome = "D90_DTH", pred1 = "lille.mort",
                       pred2 = "lille.mort.new", xstart = 0.25, xstop = 0.75,
                       step = 0.05)
pvalue.lille <- NULL

for(i in 1:length(boot.diff.lille$t0)){
  pvalue.lille <- c(pvalue.lille, mean(abs(boot.diff.lille$t[,i] - boot.diff.lille$t0[i]) > abs(boot.diff.lille$t0[i])))
}

# Append p-values together in a dataframe for easy accessibility
pvalues_oldnew <- data.frame("threshold probability" = seq(from = 0.25, to = 0.75, by = 0.05),
                              "p-values MELD" = pvalue.meld,
                              "p-values MELD2" = pvalue.meld2,
                              "p-values CLIF" = pvalue.clif,
                              "p-values Lille" = pvalue.lille)

# Also show plots that illustrate how the models compare (to interpret p-values)
meld_dca <- dca(data = test.data.c, outcome = "D90_DTH", predictors = c("meld.mort", "meld.mort2", "meld.mort.new"), xstop = 0.75)
clif_dca <- dca(data = test.data.c, outcome = "D90_DTH", predictors = c("clif.mort", "clif.mort.new"), xstop = 0.75)
lille_dca <- dca(data = test.data.c, outcome = "D90_DTH", predictors = c("lille.mort", "lille.mort.new"), xstop = 0.75)

nb_data_meld <- meld_dca$net.benefit
nb_data_clif <- clif_dca$net.benefit
nb_data_lille <- lille_dca$net.benefit

# NB plot for the old and new MELD scores
plot(nb_data_meld$threshold, nb_data_meld$none, type = "l", lwd = 2, xlab = "Threshold mortality probability", ylab = "Net benefit", ylim = c(-0.05, 0.25))
lines(nb_data_meld$threshold, nb_data_meld$all, type = "l", col = 8, lwd = 2)
lines(nb_data_meld$threshold, nb_data_meld$meld.mort, type = "l", col = "darkred", lwd = 2)
lines(nb_data_meld$threshold, nb_data_meld$meld.mort2, type = "l", col = "orange", lwd = 2)
lines(nb_data_meld$threshold, nb_data_meld$meld.mort.new, type = "l", col = "darkgreen", lwd = 2)
# Add a legend
legend("topright", cex = 0.8, legend = c("Treat none", "Treat all", "Original MELD_1", "Original MELD_2", "Updated MELD"),
       col = c(17, 8, "darkred", "orange", "darkgreen"), lwd = c(2, 2, 2, 2, 2, 2, 2))

# NB plot for the old and new CLIF-C ACLF scores
plot(nb_data_clif$threshold, nb_data_clif$none, type = "l", lwd = 2, xlab = "Threshold mortality probability", ylab = "Net benefit", ylim = c(-0.05, 0.25))
lines(nb_data_clif$threshold, nb_data_clif$all, type = "l", col = 8, lwd = 2)
lines(nb_data_clif$threshold, nb_data_clif$clif.mort, type = "l", col = "darkred", lwd = 2)
lines(nb_data_clif$threshold, nb_data_clif$clif.mort.new, type = "l", col = "darkgreen", lwd = 2)
# Add a legend
legend("topright", cex = 0.8, legend = c("Treat none", "Treat all", "Original CLIF-C ACLF", "Updated CLIF-C ACLF"),
       col = c(17, 8, "darkred", "darkgreen"), lwd = c(2, 2, 2, 2, 2, 2))

# NB plot for the old and new Lille scores
plot(nb_data_lille$threshold, nb_data_lille$none, type = "l", lwd = 2, xlab = "Threshold mortality probability", ylab = "Net benefit", ylim = c(-0.05, 0.25))
lines(nb_data_lille$threshold, nb_data_lille$all, type = "l", col = 8, lwd = 2)
lines(nb_data_lille$threshold, nb_data_lille$lille.mort, type = "l", col = "darkred", lwd = 2)
lines(nb_data_lille$threshold, nb_data_lille$lille.mort.new, type = "l", col = "darkgreen", lwd = 2)
# Add a legend
legend("topright", cex = 0.8, legend = c("Treat none", "Treat all", "Original Lille", "Updated Lille"),
       col = c(17, 8, "darkred", "darkgreen"), lwd = c(2, 2, 2, 2, 2, 2))


