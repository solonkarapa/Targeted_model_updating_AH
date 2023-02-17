######
# This script reads in the stph data and creates workable data frames
library(tidyverse)
library(dplyr)

# Load data 
path <- "/Users/work/IDrive-Sync/Projects/MIMAH/data"
load(paste0(path, "/stopah_plus_scores.Rdata"))
stph <- rename(data_stph_prelim3)

# Rename relevant variables
stph <- rename(stph, Bilirubin.mg.dl = Bilirubin.Merged..mg.dL..Merged..calc.)
stph <- rename(stph, Creatinine = Creatinine...Merged)
stph <- rename(stph, Albumin = Albumin...Merged)
stph <- rename(stph, WBC = WBC...Merged) # White blood count
stph <- rename(stph, INR = INR...Merged.clinical.and.calc) 
stph <- rename(stph, protime = Prothrombin.Time..patient....Merged) # Prothrombin time
stph <- rename(stph, HE = Hepatic.Encephalopathy...Merged) # Brain function
stph <- rename(stph, Sodium = Sodium...Merged)
    
# Create survival variable
stph$D90_surv <- 1 - stph$D90_DTH

# Transform creatinine into mg/dl
stph$Creatinine.mg.dl <- 0.0113*stph$Creatinine

# Handle missing data and create data frames for each of the prognostic scores
stph.meld <- stph[complete.cases(stph$Bilirubin.mg.dl, stph$Creatinine, stph$INR, stph$Sodium),]
stph.clif <- stph[complete.cases(stph$Bilirubin.mg.dl, stph$Creatinine, stph$INR, stph$WBC, 
                                 stph$HE, stph$MAP),]
stph.lille <- stph[complete.cases(stph$Bilirubin.day.7, stph$Bilirubin.Merged, stph$Creatinine,
                                  stph$protime, stph$Albumin),]


#path <- "/Users/work/IDrive-Sync/Projects/MIMAH/code/AH_code/pre-analysis/"
#setwd(path)
#save(stph, file = "full_sample.Rdata")

# Create table of descriptive statistics and number of missing values per variable
library(table1)
# Some factor variables
stph$HE_f <- as.factor(stph$HE) 
stph$D90_surv_f <- as.factor(stph$D90_surv)
stph$Gender_f <- as.factor(stph$Gender)
tb1 <- table1::table1(~ Bilirubin.mg.dl + Creatinine.mg.dl + Albumin + WBC + protime + 
                 INR + HE_f + Bilirubin.day.7 + Gender_f + D90_surv_f + Sodium, data = stph)

library(xtable)
xtable(as_tibble(tb1))

