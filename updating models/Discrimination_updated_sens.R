
library(pROC)
library(dplyr)
library(ggplot2)
#library(ggROC)

# load data
path_data <- "/Users/work/IDrive-Sync/Projects/MIMAH/code/AH_code/updating models"
load(paste0(path_data, "/recalibrated_models_sens.Rdata"))

#############################################   
############### Calculate AUCs  #############
#############################################  

# MELD (survival function 1 and 2 do not matter here)
roc_meld <- roc(test.meld$D90_surv, test.meld$meld.surv.updated)

# MELD 3.0
#roc_meld3 <- roc(stph.meld$D90_DTH, stph.meld$MELD3.surv)

# MELD from VanDerwerken et al 2021
#roc_meld.VanDerwerken <- roc(stph.meld$D90_DTH, stph.meld$MELD_Van)

# Lille
roc_lille <- roc(test.lille$D90_surv, test.lille$lille.surv.updated)

# CLIF-C ACLF
roc_clif <- roc(test.clif$D90_surv, test.clif$clif.surv.updated)

#############################################   
###################### Plots  ###############
############################################# 
roc.list <- list("MELD" = roc_meld, 
                 #"MELD VanDerwerken" = roc_meld.VanDerwerken,
                 "CLIF-C ACLF" = roc_clif, 
                 "Lille" = roc_lille)
#save(roc.list, file = "ROC_updated.Rdata")

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


pl <- ggroc(roc.list) + 
    facet_grid(. ~ name) +
    theme_minimal() + 
    geom_ribbon(data = df, aes(x = x, ymin = lower, ymax = upper, fill = name), alpha = 0.3, inherit.aes = F) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed") + 
    labs(x = "Specificity", y = "Sensitivity") + 
    coord_equal() +
    theme_classic() +
    theme(legend.position = "none") 

data_wide$name <- data_wide$condition
data_wide$low_CL <- round(data_wide$low_CL, 2)
data_wide$mean <- round(data_wide$mean, 2)
data_wide$upper_CL <- round(data_wide$upper_CL, 2)
pl + geom_text(data = data_wide, aes(0.05, 0.25, label = paste0("AUC (95% CI): ", mean, " (", low_CL, "-", upper_CL, ")" ), 
                                     hjust = 1), col = "black")

#############################################   
###################### AUC  #################
#############################################
# CI
#auc_meld <- auc(roc_meld)
#auc_meld_ci <- ci.auc(roc_meld) # Confidence intervals
#roc_plot(stph.meld, "D90_surv", "MELD.surv", ci = TRUE, plot_title = "ROC curve for the MELD score")

df_AUC <- as.data.frame(map_dfr(roc.list, ci.auc))
rownames(df_AUC) <- c("low_CL", "mean", "upper_CL")

df_AUC2 <- tibble::rownames_to_column(df_AUC, var = "AUC")

df3 <- gather(df_AUC2, condition, measurement, MELD:Lille, factor_key = TRUE)
data_wide <- spread(df3, AUC, measurement) %>% arrange(mean)

# reorder factor levels
data_wide$condition <- fct_reorder(data_wide$condition, data_wide$mean)

ggplot(data_wide, aes(x = mean, y = condition, col = condition)) +
    geom_point(lwd = 2)  + 
    coord_cartesian(xlim = c(0.5, 0.86)) +
    geom_errorbar(aes(xmin = low_CL, xmax = upper_CL), 
                  alpha = 1, show.legend = F, lwd = 1, width = 0.5) + 
    labs(y = "Score", col = "Score", x = "AUC with 95% limits") +
    theme_classic() 







