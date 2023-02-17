
#############################################   
############### Sensitivity Analysis  #######
#############################################

# load full data 
path_data1 <- "/Users/work/IDrive-Sync/Projects/MIMAH/code/AH_code/pre-analysis/"
load(paste0(path_data1, "full_sample.Rdata"))

# load data with original models 
path_data2 <- "/Users/work/IDrive-Sync/Projects/MIMAH/code/AH_code/original models/"
load(paste0(path_data2, "original_models.Rdata"))

#split data on full stph sample and then select the complete cases per model

# Split full sample in training and test observations
## make sure to use the same seed and fraction as in recalibrate.R script 
fraction <- 0.8

#dt <- sort(sample(nrow(stph), nrow(stph) * fraction))
#test.data <- stph[dt,]
#train.data <- stph[-dt,]

#############################################
############### MELD score ##################
#############################################
# MELD score re-calibration - subset data
#test.meld <- stph.meld[stph.meld$Subject %in% test.data$Subject,]
#train.meld <- stph.meld[stph.meld$Subject %in% train.data$Subject,]

set.seed(111) # Set random seed
dt <- sort(sample(nrow(stph.meld), nrow(stph.meld) * fraction))
test.meld <- stph.meld[dt,]
train.meld <- stph.meld[-dt,]

# Fit a logistic regression model on the training data
meld_regr <- glm(D90_surv ~ MELD.calc, family = binomial("logit"), data = train.meld)

# Save the regression coefficients (alpha and beta)
ic_meld <- meld_regr$coefficients[1]
slope_meld <- meld_regr$coefficients[2]

# Update scores on training set
train.meld$updated.meld <- ic_meld + slope_meld*train.meld$MELD.calc
train.meld$meld.surv.updated <- 1/(1 + exp(-train.meld$updated.meld))

# Update scores on test set
test.meld$updated.meld <- ic_meld + slope_meld*test.meld$MELD.calc
test.meld$meld.surv.updated <- 1/(1 + exp(-test.meld$updated.meld))

#############################################
############### Lille score #################
#############################################
# Lille score re-calibration - subset data
#test.lille <- stph.lille[stph.lille$Subject %in% test.data$Subject,]
#train.lille <- stph.lille[stph.lille$Subject %in% train.data$Subject,]

set.seed(111)
dt <- sort(sample(nrow(stph.lille), nrow(stph.lille) * fraction))
test.lille <- stph.lille[dt,]
train.lille <- stph.lille[-dt,]

# Run logistic regression and save regression coefficients
lille_regr <- glm(D90_surv ~ LILLE, family = binomial("logit"), data = train.lille)

# Save the regression coefficients (alpha and beta)
ic_lille <- lille_regr$coefficients[1]
slope_lille <- lille_regr$coefficients[2]

# Update scores on training set
train.lille$updated.lille <- ic_lille + slope_lille*train.lille$LILLE
train.lille$lille.surv.updated <- 1/(1 + exp(-train.lille$updated.lille))

# Update scores on test set
test.lille$updated.lille <- ic_lille + slope_lille*test.lille$LILLE
test.lille$lille.surv.updated <- 1/(1 + exp(-test.lille$updated.lille))

#############################################
############### CLIF-C ACLF score ###########
#############################################
# CLIF-C ACLF score re-calibration -  subset data
#test.clif <- stph.clif[stph.clif$Subject %in% test.data$Subject,]
#train.clif <- stph.clif[stph.clif$Subject %in% train.data$Subject,]

set.seed(111)
dt <- sort(sample(nrow(stph.clif), nrow(stph.clif) * fraction))
test.clif <- stph.clif[dt,]
train.clif <- stph.clif[-dt,]

# Run logistic regression and save regression coefficients
clif_regr <- glm(D90_surv ~ CLIF.C, family = binomial("logit"), data = train.clif)

# Save the regression coefficients (alpha and beta)
ic_clif <- clif_regr$coefficients[1]
slope_clif <- clif_regr$coefficients[2]

# Update scores on training set
train.clif$updated.clif <- ic_clif + slope_clif*train.clif$CLIF.C
train.clif$clif.surv.updated <- 1/(1 + exp(-train.clif$updated.clif))

# Assess performance on the test set
test.clif$updated.clif <- ic_clif + slope_clif*test.clif$CLIF.C
test.clif$clif.surv.updated <- 1/(1 + exp(-test.clif$updated.clif))

# overall sample size per model 
nrow(stph.meld); nrow(stph.lille); nrow(stph.clif)

# train/test sample size
model <- c("MELD", "Lille", "CLIF")
train_size <- c(nrow(train.meld), nrow(train.lille), nrow(train.clif))
test_size <- c(nrow(test.meld), nrow(test.lille), nrow(test.clif))
df <- data.frame(model, train_size, test_size)
df

library(xtable)
xtable(t(df))

#path_to_save <- "/Users/work/IDrive-Sync/Projects/MIMAH/code/AH_code/updating models/"
#setwd(path_to_save)
#save(test.meld, test.lille, test.clif, file = "recalibrated_models_sens.Rdata")

