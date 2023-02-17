# This script performs all calculations of model calibration
library(ggplot2)
library(survminer) # for plotting theme

path_funs <- "/Users/work/IDrive-Sync/Projects/MIMAH/code/funs"
source(paste0(path_funs, "/calibration_fun.R"))

path_data <- "/Users/work/IDrive-Sync/Projects/MIMAH/code/AH_code/updating models"
load(paste0(path_data, "/recalibrated_models_default.Rdata"))

#############################################   
############### Calculate calibration  ######
#############################################   
# MELD
cal_MELD.surv <- calibration(test.data$meld.surv.updated, y = test.data$D90_surv)
cal_MELD.surv$Score <- "MELD"

# Lille
cal_Lille <- calibration(test.data$lille.surv.updated, y = test.data$D90_surv)
cal_Lille$Score <- "Lille"

# CLIF-C ACLF
cal_CLIF <- calibration(test.data$clif.surv.updated, y = test.data$D90_surv)
cal_CLIF$Score <- "CLIF-C ACLF"

# combine dfs
df_cal <- rbind(cal_MELD.surv, cal_Lille, cal_CLIF)

#############################################   
###################### Plots  ###############
############################################# 

# plot without ribbon and without MELD 3.0
df_cal %>%
    ggplot(., aes(x = pred, y = obs, col = Score)) +
    geom_line(lwd = 1)  + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    #geom_ribbon(aes(ymin = lower, ymax = upper, linetype = NA), 
    #            alpha = 0.3, show.legend = F) + 
    #scale_fill_manual("", values = col) + 
    #scale_color_manual(name = "Score", values = col) + 
    xlim(0, 1) + 
    ylim(0, 1) + 
    ylab("Observed proportion") + 
    xlab("Predicted probability") + 
    theme_classic() 

# plot with ribbon 
df_cal %>%
    ggplot(., aes(x = pred, y = obs, col = Score)) +
    geom_line(lwd = 1)  + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = Score, linetype = NA),  
                alpha = 0.3, show.legend = F) + 
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    #scale_fill_manual("", values = col) + 
    #scale_color_manual(name = "Score", values = col) + 
    facet_grid(. ~ Score) +
    coord_equal() +
    xlim(0, 1) + 
    ylim(0, 1) + 
    ylab("Observed survival proportion") + 
    xlab("Predicted survival probability") + 
    theme_classic2() 
 
