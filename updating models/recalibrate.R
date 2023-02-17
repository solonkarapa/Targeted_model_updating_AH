# This script performs recalibration of models using the method proposed by Steyerberg et al. (2004).
# For each prognostic model (MELD, Lille and CLIF-C ACLF), a logistic regression model is used to 
# find recalibration parameters and the survival probability is re-calculated using these parameters. 
library(rms)
library(pROC)

# load data
path_data <- "/Users/work/IDrive-Sync/Projects/MIMAH/code/AH_code/original models/"
load(paste0(path_data, "complete_cases_models.Rdata"))

# Set random seed
set.seed(111)

######
# Split complete-case sample in training and test observations
fraction <- 0.8
dt <- sort(sample(nrow(stph.c), nrow(stph.c)*fraction))
test.data <- stph.c[dt,]
train.data <- stph.c[-dt,]

#############################################
############### MELD score ##################
#############################################
# Fit a logistic regression model on the training data
meld_regr <- glm(D90_surv ~ MELD.calc, family = binomial("logit"), data = train.data)

# Save the regression coefficients (alpha and beta)
ic_meld_c <- meld_regr$coefficients[1]
slope_meld_c <- meld_regr$coefficients[2]

# Calculate adjusted MELD score and survival probability on training set
train.data$updated.meld <- ic_meld_c + slope_meld_c*train.data$MELD.calc
train.data$meld.surv.updated <- 1/(1 + exp(-train.data$updated.meld))

# Calculate adjusted scores on the test set and assess performance
test.data$updated.meld <- ic_meld_c + slope_meld_c*test.data$MELD.calc
test.data$meld.surv.updated <- 1/(1 + exp(-test.data$updated.meld))

#############################################
################### Lille score #############
#############################################
# Run logistic regression and save regression coefficients
lille_regr <- glm(D90_surv ~ LILLE, family = binomial("logit"), data = train.data)

ic_lille_c <- lille_regr$coefficients[1]
slope_lille_c <- lille_regr$coefficients[2]

# Update scores in training set
train.data$updated.lille <- ic_lille_c + slope_lille_c*train.data$LILLE
train.data$lille.surv.updated <- 1/(1 + exp(-train.data$updated.lille))

# Assess performance on the test set
test.data$updated.lille <- ic_lille_c + slope_lille_c*test.data$LILLE
test.data$lille.surv.updated <- 1/(1 + exp(-test.data$updated.lille))

#############################################
######### CLIF-C ACLF Score #################
#############################################
# Run logistic regression and save regression coefficients
clif_regr <- glm(D90_surv ~ CLIF.C, family = binomial("logit"), data = train.data)

ic_clif_c <- clif_regr$coefficients[1]
slope_clif_c <- clif_regr$coefficients[2]

# Update scores in training set
train.data$updated.clif <- ic_clif_c + slope_clif_c*train.data$CLIF.C
train.data$clif.surv.updated <- 1/(1 + exp(-train.data$updated.clif))

# Assess performance on the test set
test.data$updated.clif <- ic_clif_c + slope_clif_c*test.data$CLIF.C
test.data$clif.surv.updated <- 1/(1 + exp(-test.data$updated.clif))

#############################################
######### collect all params ################
#############################################
df_meld <- data.frame(intercept = ic_meld_c, slope = slope_meld_c)
df_meld$Score <- "MELD" 

df_lille <- data.frame(intercept = ic_lille_c, slope = slope_lille_c)
df_lille$Score <- "Lille" 

df_clif <- data.frame(intercept = ic_clif_c, slope = slope_clif_c)
df_clif$Score <- "CLIF-C ACLF" 

df <- rbind(df_meld, df_lille, df_clif)
rownames(df) <- NULL

library(xtable)
xtable(df[c(3, 1, 2)])



