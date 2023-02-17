######
# This script calculates the prognostic scores (MELD, Lille, CLIF-C ACLF) and corresponding survival probabilities.
library(tidyverse)
library(dplyr)
library(ggplot2)

#############################################
############### MELD score ##################
#############################################
#####
# First constrain INR, Creatinine, and Bilirubin to be within certain bandwidth
stph.meld$INR <- ifelse(stph.meld$INR < 1, 1, stph.meld$INR)
stph.meld$Creatinine.mg.dl.MELD <- ifelse(stph.meld$Creatinine.mg.dl < 1, 1, stph.meld$Creatinine.mg.dl)
stph.meld$Bilirubin.mg.dl.MELD <- ifelse(stph.meld$Bilirubin.mg.dl < 1, 1, stph.meld$Bilirubin.mg.dl)
stph.meld$Creatinine.mg.dl.MELD <- ifelse(stph.meld$Creatinine.mg.dl > 4, 4, stph.meld$Creatinine.mg.dl)

# Calculate basic MELD score
stph.meld$MELD.calc <- 0.378*log(stph.meld$Bilirubin.mg.dl.MELD) + 
    1.120*log(stph.meld$INR) + 
    0.957*log(stph.meld$Creatinine.mg.dl.MELD) 

# Round score to the nearest tenth
stph.meld$MELD.calc <- round(stph.meld$MELD.calc, 1)

# Multiply score by ten
stph.meld$MELD.calc <- 10*stph.meld$MELD.calc

# Constrain serum sodium to be within certain bandwidth
stph.meld$Sodium <- ifelse(stph.meld$Sodium < 125, 125, stph.meld$Sodium)
stph.meld$Sodium <- ifelse(stph.meld$Sodium > 137, 137, stph.meld$Sodium)

# Recalculate to get the MELD-Na score
for(i in 1:nrow(stph.meld)){
  if(stph.meld$MELD.calc[i] > 11){
    temp <- stph.meld$MELD.calc[i] + 1.32*(137 - stph.meld$Sodium[i]) - 
      0.033*stph.meld$MELD.calc[i]*(137 - stph.meld$Sodium[i])
    stph.meld$MELD.calc[i] <- temp
  } else {
    stph.meld$MELD.calc[i] <- stph.meld$MELD.calc[i]
  }
}

# Round to nearest integer
stph.meld$MELD.calc <- round(stph.meld$MELD.calc, 0)

# Constrain MELD to be between 6 and 40
stph.meld$MELD.calc <- ifelse(stph.meld$MELD.calc < 6, 6, stph.meld$MELD.calc)
stph.meld$MELD.calc <- ifelse(stph.meld$MELD.calc > 40, 40, stph.meld$MELD.calc)

# Calculate 90-day survival probability based on MELD-Na score (2 different functions)
stph.meld$MELD.surv <- 0.707^(exp((stph.meld$MELD.calc/10) - 1.127)) 
stph.meld$MELD.surv2 <- 0.98465^(exp(0.1635*(stph.meld$MELD.calc - 10)))  

# Extract 90-day survival probability from VanDerwerke et al, 2021 
library(readxl)
MELD_VanDerwerken <- read_excel("~/IDrive-Sync/Projects/MIMAH/data/MELD_VanDerwerken.xlsx")

for(i in 1:nrow(stph.meld)){
    stph.meld$MELD_Van[i] <- MELD_VanDerwerken[stph.meld$MELD.calc[i] - 5, ]$SURV
}

####
# Calculate the MELD 3.0 score (sex-adjusted version of the MELD)
# Constrain Albumin 
stph.meld$Albumin.MELD <- ifelse(stph.meld$Albumin < 1.5, 1.5, stph.meld$Albumin)
stph.meld$Albumin.MELD <- ifelse(stph.meld$Albumin > 3.5, 3.5, stph.meld$Albumin)

# Calculate MELD 3.0, round to nearest integer, and calculate corresponding survival
stph.meld$MELD3 <- 1.33*stph.meld$Gender + 4.56*log(stph.meld$Bilirubin.mg.dl.MELD) +
  0.82*(137 - stph.meld$Sodium) - (0.24*(137 - stph.meld$Sodium)*log(stph.meld$Bilirubin.mg.dl.MELD)) +
  9.09*log(stph.meld$INR) + 11.14*log(stph.meld$Creatinine.mg.dl.MELD) + 1.85*(3.5 - stph.meld$Albumin.MELD) -
  (1.83*(3.5 - stph.meld$Albumin.MELD)*log(stph.meld$Creatinine.mg.dl.MELD)) + 6
  
stph.meld$MELD3 <- round(stph.meld$MELD3, 0)
stph.meld$MELD3.surv <- 0.946^(exp(0.17698*stph.meld$MELD3 - 3.56))

#############################################
######### CLIF-C ACLF Score #################
#############################################
#####
# Start by calculating organ-failure sub-scores
stph.clif$liver.score <- ifelse(stph.clif$Bilirubin.mg.dl < 6, 1, 2)
stph.clif$liver.score <- ifelse(stph.clif$Bilirubin.mg.dl < 12, stph.clif$liver.score, 3)

stph.clif$kidney.score <- ifelse(stph.clif$Creatinine.mg.dl < 2, 1, 2)
stph.clif$kidney.score <- ifelse(stph.clif$Creatinine.mg.dl < 3.5, stph.clif$kidney.score, 3) 

stph.clif$brain.score <- ifelse(stph.clif$HE == 0, 1, 2)
stph.clif$brain.score <- ifelse(stph.clif$HE > 2, 3, stph.clif$brain.score)

stph.clif$coag.score <- ifelse(stph.clif$INR < 2, 1, 2)
stph.clif$coag.score <- ifelse(stph.clif$INR < 2.5, stph.clif$coag.score, 3)

stph.clif$circ.score <- ifelse(stph.clif$MAP >= 70, 1, 2)

# Calculate CLIF-OF score, which is the sum of sub-scores
stph.clif$CLIF.OF <- stph.clif$liver.score + stph.clif$kidney.score + 
  stph.clif$brain.score + stph.clif$coag.score + stph.clif$circ.score + 1 # respiratory score equal to 1 for all patients

# Calculate the CLIF-C ACLF score and corresponding survival
stph.clif$CLIF.C <- 10*(0.33*stph.clif$CLIF.OF + 0.04*stph.clif$Age + 0.63*log(stph.clif$WBC) - 2)
stph.clif$CLIF.surv <- exp(-0.0079 * exp(0.0869*stph.clif$CLIF.C))

#############################################
################### Lille score #############
#############################################
#####
# First create renal insufficiency dummy
stph.lille$ren.insuf <- ifelse(stph.lille$Creatinine.mg.dl < 1.3, 0, 1)

# Calculate the change in bilirubin
stph.lille$delta.bili <- stph.lille$Bilirubin.Merged - stph.lille$Bilirubin.day.7

# Calculate the Lille score and corresponding survival probability
stph.lille$LILLE <- 3.19 - 0.101*stph.lille$Age + 0.147*stph.lille$Albumin + 0.0165*stph.lille$delta.bili - 
  0.206*stph.lille$ren.insuf - 0.0065*stph.lille$Bilirubin.Merged - 0.0096*stph.lille$protime
stph.lille$Lille.surv <- 1 - (exp(-stph.lille$LILLE)/(1 + exp(-stph.lille$LILLE)))

#####
# Create dataframe of all complete cases including prognostic scores (used later)
stph.c <- merge(stph.meld, stph.clif, by = "Subject")
stph.c <- merge(stph.c, stph.lille, by = "Subject")


# 
#path <- "/Users/work/IDrive-Sync/Projects/MIMAH/code/AH_code/original models/"
#setwd(path)
#save(stph.meld, stph.lille, stph.clif, file = "original_models.Rdata")
#save(stph.c, file = "complete_cases_models.Rdata")
# save(stph, file = "full_sample.Rdata")

# Tabulate the calculated prognostic scores and survival probabilities
#library(table1)
#table1::table1(~MELD.calc + MELD.surv + MELD.surv2 + MELD3.surv, data = stph.meld)
#table1::table1(~LILLE + Lille.surv, data = stph.lille)
#table1::table1(~CLIF.C + CLIF.surv, data = stph.clif)


