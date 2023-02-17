# R Script for sub-group analysis
library(rms)
library(pROC)
library(runway)
source("nb_diff.R")
source("dca.R")

######
# Comparison of ROC between MELD3.0 and original MELD
#roc_meld_u <- roc(test.meld$D90_surv, test.meld$meld.surv.updated)
#roc_meld3 <- roc(stph.meld$D90_DTH, stph.meld$MELD3.surv)

#roc.test(roc_meld3, roc_meld)
#roc.test(roc_meld3, roc_meld_u)


######
# Sub-group based on sex 
# Calculations on the MELD3.0 score - first create two separate dfs per sex
stph.meld.males <- stph.meld[stph.meld$Gender == 0,]
stph.meld.females <- stph.meld[stph.meld$Gender == 1,]

# Calibration plots and calculations for males/females separately
val.prob(stph.meld.males$MELD3.surv, stph.meld.males$D90_surv, 
         xlab = "Predicted survival probability", ylab = "Actual survival probability", 
         pl = TRUE, smooth = FALSE, logistic.cal = TRUE, legendloc = FALSE, statloc = FALSE)

val.prob(stph.meld.females$MELD3.surv, stph.meld.females$D90_surv, 
         xlab = "Predicted survival probability", ylab = "Actual survival probability", 
         pl = TRUE, smooth = FALSE, logistic.cal = TRUE, legendloc = FALSE, statloc = FALSE)

# Discrimination for males/females separately
# Create ROC plot and calculate AUC for MELD3.0 per gender
roc_meld3_m <- roc(stph.meld.males$D90_DTH, stph.meld.males$MELD3.surv)
auc_meld3_m <- auc(roc_meld3_m)
auc_meld3_ci_m <- ci.auc(roc_meld3_m) # Confidence intervals

roc_plot(stph.meld.males, "D90_surv", "MELD3.surv", ci = TRUE, plot_title = "ROC curve for the MELD score")

roc_meld3_f <- roc(stph.meld.females$D90_DTH, stph.meld.females$MELD3.surv)
auc_meld3_f <- auc(roc_meld3_f)
auc_meld3_ci_f <- ci.auc(roc_meld3_f) # Confidence intervals

roc_plot(stph.meld.females, "D90_surv", "MELD3.surv", ci = TRUE, plot_title = "ROC curve for the MELD score")

######
# Calculations on the updated MELD score
test.males <- test.data.c[test.data.c$Gender == 0,]
test.females <- test.data.c[test.data.c$Gender == 1,]

# Calibration 
val.prob(test.males$meld.surv.updated, test.males$D90_surv, 
         pl = TRUE, smooth = FALSE, logistic.cal = TRUE, 
         xlab = "Predicted survival probability", ylab = "Actual survival probability", 
         legendloc = FALSE, statloc = FALSE)

val.prob(test.females$meld.surv.updated, test.females$D90_surv, 
         pl = TRUE, smooth = FALSE, logistic.cal = TRUE,
         xlab = "Predicted survival probability", ylab = "Actual survival probability", 
         legendloc = FALSE, statloc = FALSE)

# Discrimination
roc_meld_m <- roc(test.males$D90_DTH, test.males$meld.surv.updated)
auc_meld_m <- auc(roc_meld_m)
auc_meld_ci_m <- ci.auc(roc_meld_m) # Confidence intervals

roc_plot(test.males, "D90_surv", "meld.surv.updated", ci = TRUE, plot_title = "ROC curve for the MELD score (males)")

roc_meld_f <- roc(test.females$D90_DTH, test.females$meld.surv.updated)
auc_meld_f <- auc(roc_meld_f)
auc_meld_ci_f <- ci.auc(roc_meld_f) # Confidence intervals

roc_plot(test.females, "D90_surv", "meld.surv.updated", ci = TRUE, plot_title = "ROC curve for the MELD score (females)")

######
# Calculations on the updated CLIF-C ACLF score
# Calibration 
val.prob(test.males$clif.surv.updated, test.males$D90_surv, 
         xlab = "Predicted survival probability", ylab = "Actual survival probability", 
         pl = TRUE, smooth = FALSE, logistic.cal = TRUE, legendloc = FALSE, statloc = FALSE)

val.prob(test.females$clif.surv.updated, test.females$D90_surv, 
         xlab = "Predicted survival probability", ylab = "Actual survival probability", 
         pl = TRUE, smooth = FALSE, logistic.cal = TRUE, legendloc = FALSE, statloc = FALSE)

# Discrimination
roc_clif_m <- roc(test.males$D90_DTH, test.males$clif.surv.updated)
auc_clif_m <- auc(roc_clif_m)
auc_clif_ci_m <- ci.auc(roc_clif_m) # Confidence intervals

roc_plot(test.males, "D90_surv", "clif.surv.updated", ci = TRUE, plot_title = "ROC curve for the MELD score (males)")

roc_clif_f <- roc(test.females$D90_DTH, test.females$clif.surv.updated)
auc_clif_f <- auc(roc_clif_f)
auc_clif_ci_f <- ci.auc(roc_clif_f) # Confidence intervals

roc_plot(test.females, "D90_surv", "meld.surv.updated", ci = TRUE, plot_title = "ROC curve for the MELD score (females)")

######
# Calculations on the updated Lille score
# Calibration 
val.prob(test.males$lille.surv.updated, test.males$D90_surv, 
         xlab = "Predicted survival probability", ylab = "Actual survival probability", 
         pl = TRUE, smooth = FALSE, logistic.cal = TRUE, legendloc = FALSE, statloc = FALSE)

val.prob(test.females$lille.surv.updated, test.females$D90_surv,
         xlab = "Predicted survival probability", ylab = "Actual survival probability", 
         pl = TRUE, smooth = FALSE, logistic.cal = TRUE, legendloc = FALSE, statloc = FALSE)

# Discrimination
roc_lille_m <- roc(test.males$D90_DTH, test.males$lille.surv.updated)
auc_lille_m <- auc(roc_lille_m)
auc_lille_ci_m <- ci.auc(roc_lille_m) # Confidence intervals

roc_lille_f <- roc(test.females$D90_DTH, test.females$lille.surv.updated)
auc_lille_f <- auc(roc_lille_f)
auc_lille_ci_f <- ci.auc(roc_lille_f) # Confidence intervals

######
# p-value comparisons between male and female sub-samples

# Formally compare the c-statistics across models using bootstrap method
compareroc.meld3 <- roc.test(roc_meld3_f, roc_meld3_m) 
compareroc.meld <- roc.test(roc_meld_f, roc_meld_m) 
compareroc.clif <- roc.test(roc_clif_f, roc_clif_m)
compareroc.lille <- roc.test(roc_lille_f, roc_lille_m)

# Tabulate the p-values
roc_pvalues_mf <- c(compareroc.meld3$p.value, compareroc.meld$p.value, compareroc.clif$p.value, compareroc.lille$p.value)
names(roc_pvalues_mf) <- c("MELD 3.0", "MELD", "CLIF", "Lille")
roc_pvalues_mf

######
# Age stratification: first find the mean age and the age distribution of the data
hist(stph$Age.at.randomisation..calc.)
summary(stph$Age.at.randomisation..calc.)

summary(stph[stph$D90_DTH == 0,]$Age.at.randomisation..calc.) # 47.4 median age
summary(stph[stph$D90_DTH == 1,]$Age.at.randomisation..calc.) # 53.1 median age

#####
# MELD
# stratify by the median age of 49 years
test.y <- test.data.c[test.data.c$Age.at.randomisation..calc. < 49,]
test.o <- test.data.c[test.data.c$Age.at.randomisation..calc. >= 49,]

# Calibration of the age-stratified MELD
val.prob(test.y$meld.surv.updated, test.y$D90_surv, 
         xlab = "Predicted survival probability", ylab = "Actual survival probability", 
         pl = TRUE, smooth = FALSE, logistic.cal = TRUE, legendloc = FALSE, statloc = FALSE)

val.prob(test.o$meld.surv.updated, test.o$D90_surv, 
         xlab = "Predicted survival probability", ylab = "Actual survival probability", 
         pl = TRUE, smooth = FALSE, logistic.cal = TRUE, legendloc = FALSE, statloc = FALSE)

# Discrimination for the MELD
roc_meld_y <- roc(test.y$D90_DTH, test.y$meld.surv.updated)
auc_meld_y <- auc(roc_meld_y)
auc_meld_ci_y <- ci.auc(roc_meld_y) # Confidence intervals

roc_plot(test.y, "D90_surv", "meld.surv.updated", ci = TRUE, plot_title = "ROC curve for the MELD score (males)")

roc_meld_o <- roc(test.o$D90_DTH, test.o$meld.surv.updated)
auc_meld_o <- auc(roc_meld_o)
auc_meld_ci_o <- ci.auc(roc_meld_o) # Confidence intervals

roc_plot(test.o, "D90_surv", "meld.surv.updated", ci = TRUE, plot_title = "ROC curve for the MELD score (females)")

#####
# LILLE
# Calibration of the age-stratified Lille
val.prob(test.y$lille.surv.updated, test.y$D90_surv, 
         xlab = "Predicted survival probability", ylab = "Actual survival probability", 
         pl = TRUE, smooth = FALSE, logistic.cal = TRUE, legendloc = FALSE, statloc = FALSE)

val.prob(test.o$lille.surv.updated, test.o$D90_surv, 
         xlab = "Predicted survival probability", ylab = "Actual survival probability", 
         pl = TRUE, smooth = FALSE, logistic.cal = TRUE, legendloc = FALSE, statloc = FALSE)

# Discrimination for the Lille
roc_lille_y <- roc(test.y$D90_DTH, test.y$lille.surv.updated)
auc_lille_y <- auc(roc_lille_y)
auc_lille_ci_y <- ci.auc(roc_lille_y) # Confidence intervals

roc_lille_o <- roc(test.o$D90_DTH, test.o$lille.surv.updated)
auc_lille_o <- auc(roc_lille_o)
auc_lille_ci_o <- ci.auc(roc_lille_o) # Confidence intervals

#####
# CLIF-C ACLF
# Calibration of the age-stratified CLIF
val.prob(test.y$clif.surv.updated, test.y$D90_surv, 
         xlab = "Predicted survival probability", ylab = "Actual survival probability", 
         pl = TRUE, smooth = FALSE, logistic.cal = TRUE, legendloc = FALSE, statloc = FALSE)

val.prob(test.o$clif.surv.updated, test.o$D90_surv, 
         xlab = "Predicted survival probability", ylab = "Actual survival probability", 
         pl = TRUE, smooth = FALSE, logistic.cal = TRUE, legendloc = FALSE, statloc = FALSE)

# Discrimination for the CLIF
roc_clif_y <- roc(test.y$D90_DTH, test.y$clif.surv.updated)
auc_clif_y <- auc(roc_clif_y)
auc_clif_ci_y <- ci.auc(roc_clif_y) # Confidence intervals

roc_clif_o <- roc(test.o$D90_DTH, test.o$clif.surv.updated)
auc_clif_o <- auc(roc_clif_o)
auc_clif_ci_o <- ci.auc(roc_clif_o) # Confidence intervals

######
# Formally compare the c-statistics across models using bootstrap method
compareroc.meld <- roc.test(roc_meld_y, roc_meld_o) 
compareroc.clif <- roc.test(roc_clif_y, roc_clif_o)
compareroc.lille <- roc.test(roc_lille_y, roc_lille_o)

# Tabulate the p-values
roc_pvalues_yo <- c(compareroc.meld$p.value, compareroc.clif$p.value, compareroc.lille$p.value)
names(roc_pvalues_yo) <- c("MELD", "CLIF", "Lille")
roc_pvalues_yo

######
# Clinical utility
# First make dataframes for males and females separately for overlapping observations
test.males$clif.mort <- 1 - test.males$clif.surv.updated
test.males$meld.mort <- 1 - test.males$meld.surv.updated
test.males$lille.mort <- 1 - test.males$lille.surv.updated
test.males$meld3.mort <- 1 - test.males$MELD3.surv

test.females$clif.mort <- 1 - test.females$clif.surv.updated
test.females$meld.mort <- 1 - test.females$meld.surv.updated
test.females$lille.mort <- 1 - test.females$lille.surv.updated
test.females$meld3.mort <- 1 - test.females$MELD3.surv

male_dca <- dca(data = test.males, outcome = "D90_DTH", predictors = c("clif.mort", "meld.mort", "lille.mort", "meld3.mort"), xstop = 0.75)
nb_data_m <- male_dca$net.benefit

# Create figure
plot(nb_data_m$threshold, nb_data_m$none, type = "l", lwd = 2, xlab = "Threshold mortality probability", ylab = "Net benefit", ylim = c(-0.05, 0.25))
lines(nb_data_m$threshold, nb_data_m$all, type = "l", col = 8, lwd = 2)
lines(nb_data_m$threshold, nb_data_m$meld.mort, type = "l", col = "darkblue", lwd = 2)
lines(nb_data_m$threshold, nb_data_m$clif.mort, type = "l", col = "darkgreen", lwd = 2)
lines(nb_data_m$threshold, nb_data_m$lille.mort, type = "l", col = "darkred", lwd = 2)
lines(nb_data_m$threshold, nb_data_m$meld3.mort, type = "l", col = "orange", lwd = 2)
# Add a legend
legend("topright", cex = 0.8, legend = c("Treat none", "Treat all", "MELD", "CLIF-C ACLF", "Lille", "MELD 3.0"),
       col = c(17, 8, "darkblue", "darkgreen","darkred", "orange"), lwd = c(2, 2, 2, 2, 2, 2))

female_dca <- dca(data = test.females, outcome = "D90_DTH", predictors = c("clif.mort", "meld.mort", "lille.mort", "meld3.mort"), xstop = 0.75)
nb_data_f <- female_dca$net.benefit

# Create figure
plot(nb_data_f$threshold, nb_data_f$none, type = "l", lwd = 2, xlab = "Threshold mortality probability", ylab = "Net benefit", ylim = c(-0.05, 0.25))
lines(nb_data_f$threshold, nb_data_f$all, type = "l", col = 8, lwd = 2)
lines(nb_data_f$threshold, nb_data_f$meld.mort, type = "l", col = "darkblue", lwd = 2)
lines(nb_data_f$threshold, nb_data_f$clif.mort, type = "l", col = "darkgreen", lwd = 2)
lines(nb_data_f$threshold, nb_data_f$lille.mort, type = "l", col = "darkred", lwd = 2)
lines(nb_data_f$threshold, nb_data_f$meld3.mort, type = "l", col = "orange", lwd = 2)
# Add a legend
legend("topright", cex = 0.8, legend = c("Treat none", "Treat all", "MELD", "CLIF-C ACLF", "Lille", "MELD 3.0"),
       col = c(17, 8, "darkblue", "darkgreen","darkred", "orange"), lwd = c(2, 2, 2, 2, 2, 2))

######
# Same procedure for age-based sub-groups
test.y$clif.mort <- 1 - test.y$clif.surv.updated
test.y$meld.mort <- 1 - test.y$meld.surv.updated
test.y$lille.mort <- 1 - test.y$lille.surv.updated

test.o$clif.mort <- 1 - test.o$clif.surv.updated
test.o$meld.mort <- 1 - test.o$meld.surv.updated
test.o$lille.mort <- 1 - test.o$lille.surv.updated

young_dca <- dca(data = test.y, outcome = "D90_DTH", predictors = c("clif.mort", "meld.mort", "lille.mort"), xstop = 0.75)
nb_data_y <- young_dca$net.benefit

# Create figure
plot(nb_data_y$threshold, nb_data_y$none, type = "l", lwd = 2, xlab = "Threshold mortality probability", ylab = "Net benefit", ylim = c(-0.10, 0.25))
lines(nb_data_y$threshold, nb_data_y$all, type = "l", col = 8, lwd = 2)
lines(nb_data_y$threshold, nb_data_y$meld.mort, type = "l", col = "darkblue", lwd = 2)
lines(nb_data_y$threshold, nb_data_y$clif.mort, type = "l", col = "darkgreen", lwd = 2)
lines(nb_data_y$threshold, nb_data_y$lille.mort, type = "l", col = "darkred", lwd = 2)
# Add a legend
legend("topright", cex = 0.8, legend = c("Treat none", "Treat all", "MELD", "CLIF-C ACLF", "Lille"),
       col = c(17, 8, "darkblue", "darkgreen","darkred"), lwd = c(2, 2, 2, 2, 2, 2))

old_dca <- dca(data = test.o, outcome = "D90_DTH", predictors = c("clif.mort", "meld.mort", "lille.mort"), xstop = 0.75)
nb_data_o <- old_dca$net.benefit

# Create figure
plot(nb_data_o$threshold, nb_data_o$none, type = "l", lwd = 2, xlab = "Threshold mortality probability", ylab = "Net benefit", ylim = c(-0.10, 0.25))
lines(nb_data_o$threshold, nb_data_o$all, type = "l", col = 8, lwd = 2)
lines(nb_data_o$threshold, nb_data_o$meld.mort, type = "l", col = "darkblue", lwd = 2)
lines(nb_data_o$threshold, nb_data_o$clif.mort, type = "l", col = "darkgreen", lwd = 2)
lines(nb_data_o$threshold, nb_data_o$lille.mort, type = "l", col = "darkred", lwd = 2)
# Add a legend
legend("topright", cex = 0.8, legend = c("Treat none", "Treat all", "MELD", "CLIF-C ACLF", "Lille"),
       col = c(17, 8, "darkblue", "darkgreen","darkred"), lwd = c(2, 2, 2, 2, 2, 2))





