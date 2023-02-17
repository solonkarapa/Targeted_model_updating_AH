# This script performs all calculations for assessing discrimination

library(pROC)
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(purrr)

# data
path_data <- "/Users/work/IDrive-Sync/Projects/MIMAH/code/AH_code/original models/"
load(paste0(path_data, "original_models.Rdata"))

#############################################   
############### Calculate AUCs  #############
#############################################  

# MELD (survival function 1 and 2 do not matter here)
roc_meld <- roc(stph.meld$D90_DTH, stph.meld$MELD.surv)

# MELD 3.0
roc_meld3 <- roc(stph.meld$D90_DTH, stph.meld$MELD3.surv)

# MELD from VanDerwerken et al 2021
roc_meld.VanDerwerken <- roc(stph.meld$D90_DTH, stph.meld$MELD_Van)

# Lille
roc_lille <- roc(stph.lille$D90_DTH, stph.lille$Lille.surv)

# CLIF-C ACLF
roc_clif <- roc(stph.clif$D90_DTH, stph.clif$CLIF.surv)

#############################################   
###################### Plots  ###############
############################################# 
roc.list <- list("CLIF-C ACLF" = roc_clif,
                 "Lille" = roc_lille,
                 "MELD" = roc_meld, 
                 "MELD VanDerwerken" = roc_meld.VanDerwerken)

#save(roc.list, file = "ROC_original.Rdata")

g.list <- ggroc(roc.list)

# ROC plot of all models combined
g.list +
    geom_line(lwd = 1) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed") + 
    theme_classic()

# faceting 
g.list + 
    facet_grid(. ~ name) + 
    geom_line(lwd = 1) + 
    geom_abline(slope = 1, intercept = 1, linetype = "dashed") +
    theme_classic() + 
    theme(legend.position="none") 

# add confidence bands
ci.list <- lapply(roc.list, ci.se, specificities = seq(0, 1, l = 25))

dat.ci.list <- lapply(ci.list, function(ciobj) 
    data.frame(x = as.numeric(rownames(ciobj)),
               lower = ciobj[, 1],
               upper = ciobj[, 3]))

df <- plyr::ldply(dat.ci.list, data.frame, .id = "name")

# rorder based on AUC scores
#df_cal$Score <- factor(df_cal$Score, labels = levels(data_wide$condition))

# To have all the curves of the same color, use aes="group":
#g.group <- ggroc(roc.list, aes="group")
#g.group
#g.group + facet_grid(.~name)

#############################################   
###################### AUC  ###############
#############################################
# CI
#auc_meld <- auc(roc_meld)
#auc_meld_ci <- ci.auc(roc_meld) # Confidence intervals
#roc_plot(stph.meld, "D90_surv", "MELD.surv", ci = TRUE, plot_title = "ROC curve for the MELD score")

df_AUC <- as.data.frame(map_dfr(roc.list, ci.auc))
rownames(df_AUC) <- c("low_CL", "mean", "upper_CL")

df_AUC2 <- tibble::rownames_to_column(df_AUC, var = "AUC")

df3 <- gather(df_AUC2, condition, measurement, `CLIF-C ACLF`:`MELD VanDerwerken`, factor_key = TRUE)
data_wide <- spread(df3, AUC, measurement) %>% arrange(mean)

# reorder factor levels
data_wide$condition <- fct_reorder(data_wide$condition, data_wide$mean)

p_auc <- ggplot(data_wide, aes(x = mean, y = condition, col = condition)) +
    geom_point(lwd = 2)  + 
    coord_cartesian(xlim = c(0.5, 0.86)) +
    geom_errorbar(aes(xmin = low_CL, xmax = upper_CL), 
                  alpha = 1, show.legend = F, lwd = 1, width = 0.5) + 
    labs(y = "Score", col = "Score", x = "AUC with 95% limits") +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_classic2() +
    theme(legend.position = "none") 

#####
p_roc <- ggroc(roc.list, lwd = 1.1) + 
    facet_grid(. ~ name) + 
    geom_ribbon(data = df, aes(x = x, ymin = lower, ymax = upper, fill = name), alpha = 0.3, inherit.aes = F) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed") + 
    labs(x = "Specificity", y = "Sensitivity") + 
    coord_equal() +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_classic2() +
    theme(legend.position = "none") 

data_wide$name <- data_wide$condition
data_wide$low_CL <- round(data_wide$low_CL, 2)
data_wide$mean <- round(data_wide$mean, 2)
data_wide$upper_CL <- round(data_wide$upper_CL, 2)

p_roc + geom_text(data = data_wide, 
                  aes(0.02, 0.15, label = paste0("AUC (95% CI): ", mean, " (", low_CL, "-", upper_CL, ")" ), 
                                     hjust = 1), col = "black")


