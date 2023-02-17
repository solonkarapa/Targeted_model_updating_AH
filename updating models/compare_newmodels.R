# This script performs model comparisons for the re-calibrated models

# H-L test for the updated models
library(ResourceSelection)

hoslem.test(test.data.c$D90_surv, test.data.c$meld.surv.updated)
hoslem.test(test.data.c$D90_surv, test.data.c$lille.surv.updated)
hoslem.test(test.data.c$D90_surv, test.data.c$clif.surv.updated)

# Sensitivity analysis
hoslem.test(test.meld$D90_surv, test.meld$meld.surv.updated)
hoslem.test(test.lille$D90_surv, test.lille$lille.surv.updated)
hoslem.test(test.clif$D90_surv, test.clif$clif.surv.updated)

#####
# Comparison of c-statistics using p-values
library(pROC)
roc_meld <- roc(test.data.c$D90_surv, test.data.c$meld.surv.updated)
roc_lille <- roc(test.data.c$D90_surv, test.data.c$lille.surv.updated)
roc_clif <- roc(test.data.c$D90_surv, test.data.c$clif.surv.updated)

compareroc.mc.new <- roc.test(roc_clif, roc_meld) # comparison between MELD and CLIF
compareroc.ml.new <- roc.test(roc_meld, roc_lille) # comparison between MELD and Lille
compareroc.cl.new <- roc.test(roc_clif, roc_lille) # comparison between CLIF and Lille

rocnew_pvalues <- c(compareroc.mc.new$p.value, compareroc.ml.new$p.value, compareroc.cl.new$p.value)
names(rocnew_pvalues) <- c("p-value MELD-CLIF", "p-value MELD-Lille", "p-value CLIF-Lille")

#####
# NRI on re-calibrated scores
# Define events and probability vectors
event <- test.data.c$D90_surv
p.MELD <- test.data.c$meld.surv.updated
p.LILLE <- test.data.c$lille.surv.updated
p.CLIF <- test.data.c$clif.surv.updated

# Define cut-off points
cut_lille <- 0.45
cut_meld <- 0.707^(exp(2.5 - 1.127))
cut_meld2 <- 0.98465^(exp(0.1635*(25 - 10)))  
cut_clif <- exp(-0.0079 * exp(0.0869*51))

# Calculate NRIs
library(nricens)
# MELD and Lille
NRI_ML <- nribin(event = event, p.std = p.MELD, p.new = p.LILLE, cut = cut_meld, niter = 0, updown = 'category')
# CLIF-C ACLF and Lille
NRI_CL <- nribin(event = event, p.std = p.CLIF, p.new = p.LILLE, cut = cut_clif, niter = 0, updown = 'category')
# MELD and CLIF-C ACLF
NRI_MC <- nribin(event = event, p.std = p.MELD, p.new = p.CLIF, cut = cut_meld, niter = 0, updown = 'category')
# Lille and MELD
NRI_LM <- nribin(event = event, p.std = p.LILLE, p.new = p.MELD, cut = cut_lille, niter = 0, updown = 'category')
# Lille and CLIF-C ACLF
NRI_LC <- nribin(event = event, p.std = p.LILLE, p.new = p.CLIF, cut = cut_lille, niter = 0, updown = 'category')
# CLIF-C ACLF and MELD
NRI_CM <- nribin(event = event, p.std = p.CLIF, p.new = p.MELD, cut = cut_clif, niter = 0, updown = 'category')
# MELD_2 and Lille
NRI_M2L <- nribin(event = event, p.std = p.MELD, p.new = p.LILLE, cut = cut_meld2, niter = 0, updown = 'category')
# MELD_2 and CLIF-C ACLF
NRI_M2C <- nribin(event = event, p.std = p.MELD, p.new = p.CLIF, cut = cut_meld2, niter = 0, updown = 'category')

#######
# Clinical utility
source("nb_diff.R")
source("dca.R")

test.data.c$clif.mort <- 1 - test.data.c$clif.surv.updated
test.data.c$meld.mort <- 1 - test.data.c$meld.surv.updated
test.data.c$lille.mort <- 1 - test.data.c$lille.surv.updated

updated_dca <- dca(data = test.data.c, outcome = "D90_DTH", predictors = c("clif.mort", "meld.mort", "lille.mort"), xstop = 0.75)
nb_data <- updated_dca$net.benefit

# Create figure
plot(nb_data$threshold, nb_data$none, type = "l", lwd = 2, xlab = "Threshold mortality probability", ylab = "Net benefit", ylim = c(-0.05, 0.25))
lines(nb_data$threshold, nb_data$all, type = "l", col = 8, lwd = 2)
lines(nb_data$threshold, nb_data$meld.mort, type = "l", col = "darkblue", lwd = 2)
lines(nb_data$threshold, nb_data$clif.mort, type = "l", col = "darkred", lwd = 2)
lines(nb_data$threshold, nb_data$lille.mort, type = "l", col = "orange", lwd = 2)
# Add a legend
legend("topright", cex = 0.8, legend = c("Treat none", "Treat all", "MELD", "CLIF-C ACLF", "Lille"),
       col = c(17, 8, "darkblue", "darkred", "orange"), lwd = c(2, 2, 2, 2, 2, 2))

# Tabulate NB values
table_output <- dca(data = stph.cnew, outcome = "D90_DTH", 
                    predictors = c("clif.mort", "meld.mort", "lille.mort"), 
                    xstart = 0.25, xstop = 0.75, xby = 0.10, graph = F)
table_output

# Formally compare the NB values for different values of the threshold probability using a bootstrap approach
library(boot)
set.seed(34)
R <- 500

# First the difference between clif and meld
boot.diff.cm <- boot(data=test.data.c, statistic = nb_diff,
                     R = R, outcome = "D90_DTH", pred1 = "clif.mort",
                     pred2 = "meld.mort", xstart = 0.25, xstop = 0.75,
                     step = 0.05)
pvalue.cm <- NULL

for(i in 1:length(boot.diff.cm$t0)){
  pvalue.cm <- c(pvalue.cm, mean(abs(boot.diff.cm$t[,i] - boot.diff.cm$t0[i]) > abs(boot.diff.cm$t0[i])))
}

# Now the difference between MELD and Lille
boot.diff.ml <- boot(data=test.data.c, statistic = nb_diff,
                     R = R, outcome = "D90_DTH", pred1 = "lille.mort",
                     pred2 = "meld.mort", xstart = 0.25, xstop = 0.75,
                     step = 0.05)
pvalue.ml <- NULL

for(i in 1:length(boot.diff.ml$t0)){
  pvalue.ml <- c(pvalue.ml, mean(abs(boot.diff.ml$t[,i] - boot.diff.ml$t0[i]) > abs(boot.diff.ml$t0[i])))
}

# And finally the difference between Lille and CLIF-C ACLF
boot.diff.cl <- boot(data=test.data.c, statistic = nb_diff,
                     R = R, outcome = "D90_DTH", pred1 = "clif.mort",
                     pred2 = "lille.mort", xstart = 0.25, xstop = 0.75,
                     step = 0.05)
pvalue.cl <- NULL

for(i in 1:length(boot.diff.cl$t0)){
  pvalue.cl <- c(pvalue.cl, mean(abs(boot.diff.cl$t[,i] - boot.diff.cl$t0[i]) > abs(boot.diff.cl$t0[i])))
}


# append p-values together in a df
pvalues_updated <- data.frame("threshold probability" = seq(from = 0.25, to = 0.75, by = 0.05),
                              "p-values Meld-CLIF" = pvalue.cm,
                              "p-values MELD-Lille" = pvalue.ml,
                              "p-values CLIF-Lille" = pvalue.cl)

#####
# NB sensitivity analysis: perform DCA for the other data splitting method
test.clif$clif.mort <- 1 - test.clif$clif.surv.updated
test.meld$meld.mort <- 1 - test.meld$meld.surv.updated
test.lille$lille.mort <- 1 - test.lille$lille.surv.updated

test.sets <- merge(test.clif, test.meld, by = "Subject")
test.sets <- merge(test.sets, test.lille, by = "Subject")

updated_dca <- dca(data = test.sets, outcome = "D90_DTH", predictors = c("clif.mort", "meld.mort", "lille.mort"), xstop = 0.75)
nb_data <- updated_dca$net.benefit

# Create figure
plot(nb_data$threshold, nb_data$none, type = "l", lwd = 2, xlab = "Threshold mortality probability", ylab = "Net benefit", ylim = c(-0.05, 0.25))
lines(nb_data$threshold, nb_data$all, type = "l", col = 8, lwd = 2)
lines(nb_data$threshold, nb_data$meld.mort, type = "l", col = "darkblue", lwd = 2)
lines(nb_data$threshold, nb_data$clif.mort, type = "l", col = "darkred", lwd = 2)
lines(nb_data$threshold, nb_data$lille.mort, type = "l", col = "orange", lwd = 2)
# Add a legend
legend("topright", cex = 0.8, legend = c("Treat none", "Treat all", "MELD", "CLIF-C ACLF", "Lille"),
       col = c(17, 8, "darkblue", "darkred", "orange"), lwd = c(2, 2, 2, 2, 2, 2))





