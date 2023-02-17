
# This script performs all calculations of model calibration
library(ggplot2)

path_funs <- "/Users/work/IDrive-Sync/Projects/MIMAH/code/funs"
source(paste0(path_funs, "/calibration_fun.R"))

path_data <- "/Users/work/IDrive-Sync/Projects/MIMAH/code/AH_code/updating models"
load(paste0(path_data, "/recalibrated_models_sens.Rdata"))

#############################################   
############### Calculate calibration  ######
#############################################   
# MELD
cal_MELD.surv <- calibration(test.meld$meld.surv.updated, y = test.meld$D90_surv)
cal_MELD.surv$Score <- "MELD"

# Lille
cal_Lille <- calibration(test.lille$lille.surv.updated, y = test.lille$D90_surv)
cal_Lille$Score <- "Lille"

# CLIF-C ACLF
cal_CLIF <- calibration(test.clif$clif.surv.updated, y = test.clif$D90_surv)
cal_CLIF$Score <- "CLIF-C ACLF"

#############################################   
###################### Plots  ###############
############################################# 
# plot with ribbon 
p1 <- cal_MELD.surv %>%
    ggplot(., aes(x = pred, y = obs)) +
    geom_line(linewidth = 1)  + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_ribbon(aes(ymin = lower, ymax = upper, linetype = NA),  
                alpha = 0.3, show.legend = F) + 
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    #scale_fill_manual("", values = col) + 
    #scale_color_manual(name = "Score", values = col) + 
    facet_grid(. ~ Score) +
    coord_equal() +
    xlim(0, 1) + 
    ylim(0, 1) + 
    ylab("Observed proportion") + 
    xlab("Predicted probability") + 
    theme_classic2() 

p2 <- cal_Lille %>%
    ggplot(., aes(x = pred, y = obs)) +
    geom_line(lwd = 1)  + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_ribbon(aes(ymin = lower, ymax = upper, linetype = NA),  
                alpha = 0.3, show.legend = F) + 
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    #scale_fill_manual("", values = col) + 
    #scale_color_manual(name = "Score", values = col) + 
    facet_grid(. ~ Score) +
    coord_equal() +
    xlim(0, 1) + 
    ylim(0, 1) + 
    ylab("Observed proportion") + 
    xlab("Predicted probability") + 
    theme_classic2() 

p3 <- cal_CLIF %>%
    ggplot(., aes(x = pred, y = obs)) +
    geom_line(lwd = 1)  + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_ribbon(aes(ymin = lower, ymax = upper, linetype = NA),  
                alpha = 0.3, show.legend = F) + 
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") + 
    #scale_fill_manual("", values = col) + 
    #scale_color_manual(name = "Score", values = col) + 
    facet_grid(. ~ Score) +
    coord_equal() +
    xlim(0, 1) + 
    ylim(0, 1) + 
    ylab("Observed proportion") + 
    xlab("Predicted probability") + 
    theme_classic2() 

ggarrange(p1, p2, p3, nrow = 1, ncol = 3, common.legend = TRUE)
