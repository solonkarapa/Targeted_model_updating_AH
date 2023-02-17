
library(pROC)
library(ggplot2)
library(dplyr)

# data original models
#path_data <- "/Users/work/IDrive-Sync/Projects/MIMAH/code/AH_code/original models/"
#load(paste0(path_data, "original_models.Rdata"))

# data recalibrated models
path_data_rec <- "/Users/work/IDrive-Sync/Projects/MIMAH/code/AH_code/updating models"
load(paste0(path_data_rec, "/recalibrated_models_default.Rdata"))

# funs
path_funs <- "/Users/work/IDrive-Sync/Projects/MIMAH/code/funs"
source(paste0(path_funs, "/calibration_fun.R"))

age_cut_off <- 49 # stratify by the median age of 49 years

test.y <- test.data[test.data$Age.at.randomisation..calc. < age_cut_off,]
test.o <- test.data[test.data$Age.at.randomisation..calc. >= age_cut_off,]

# choose age-group to calculate statistics 
age_group <- "o" # "y" or "o"

if(age_group == "y"){
    dataset <- test.y 
}else if(age_group == "o"){
    dataset <- test.o
}

#### Discrimination ###############
# MELD updated
roc_meld <- roc(dataset$D90_DTH, dataset$meld.surv.updated)
# CLIF - updated
roc_clif <- roc(dataset$D90_DTH, dataset$clif.surv.updated)
# Lille - updated 
roc_lille <- roc(dataset$D90_DTH, dataset$lille.surv.updated)

## Plots  ###############
roc.list <- list("MELD updated" = roc_meld, 
                 #"MELD 3.0" = roc_meld3,
                 "CLIF-C ACLF updated" = roc_clif, 
                 "Lille updated" = roc_lille)

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
    ggtitle(paste0("age_group ", age_group)) +
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

ggroc(roc.list) + 
    facet_grid(. ~ name) +
    theme_minimal() + 
    ggtitle(paste0("age_group ", age_group)) +
    geom_ribbon(data = df, aes(x = x, ymin = lower, ymax = upper, fill = name), alpha = 0.3, inherit.aes = F) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed") + 
    labs(x = "Specificity", y = "Sensitivity") + 
    coord_equal() +
    theme_classic() +
    theme(legend.position = "none") 

#############################################   
###################### AUC  #################
#############################################
df_AUC <- as.data.frame(map_dfr(roc.list, ci.auc))
rownames(df_AUC) <- c("low_CL", "mean", "upper_CL")

df_AUC2 <- tibble::rownames_to_column(df_AUC, var = "AUC")

df3 <- gather(df_AUC2, condition, measurement, `MELD updated`:`Lille updated`, factor_key = TRUE)
data_wide <- spread(df3, AUC, measurement) %>% arrange(mean)

# reorder factor levels
data_wide$condition <- fct_reorder(data_wide$condition, data_wide$mean)

p_auc_young <- ggplot(data_wide, aes(x = mean, y = condition, col = condition)) +
    #ggtitle(age_group) +
    geom_point(lwd = 2)  + 
    coord_cartesian(xlim = c(0.5, 0.86)) +
    geom_errorbar(aes(xmin = low_CL, xmax = upper_CL), 
                  alpha = 1, show.legend = F, lwd = 1, width = 0.5) + 
    labs(y = "Score", col = "Score", x = "AUC with 95% limits") +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") + 
    theme_classic() 

#### Calibration ###############
# MELD updated 
cal_meld <- calibration(dataset$meld.surv.updated, dataset$D90_surv)
cal_meld$Score <- "MELD updated"
# CLIF - updated
cal_clif <-  calibration(dataset$clif.surv.updated, dataset$D90_surv)
cal_clif$Score <- "CLIF-C ACLF updated"
# Lille - updated 
cal_lille <- calibration(dataset$lille.surv.updated, dataset$D90_surv)
cal_lille$Score <- "Lille updated"

# combine dfs
df_cal <- rbind(cal_meld, cal_clif, cal_lille)

df_cal$Score <- factor(df_cal$Score, labels = levels(data_wide$condition))

# plot without ribbon 
df_cal %>%
    ggplot(., aes(x = pred, y = obs, col = Score)) +
    ggtitle(paste0("age_group ", age_group)) +
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
p_calibration_young <- df_cal %>%
    ggplot(., aes(x = pred, y = obs, col = Score)) +
    #ggtitle(paste0("age_group ", age_group)) +
    geom_line(lwd = 1)  + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = Score, linetype = NA),  
                alpha = 0.3, show.legend = F) + 
    #scale_fill_manual("", values = col) + 
    #scale_color_manual(name = "Score", values = col) + 
    facet_grid(. ~ Score) +
    coord_equal() +
    xlim(0, 1) + 
    ylim(0, 1) + 
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") + 
    ylab("Observed survival proportion") + 
    xlab("Predicted survival probability") + 
    theme_classic() 


##### combine plots
source("/Users/work/IDrive-Sync/Projects/MIMAH/code/AH_code/grid_arrange_fun.R")
library(gtable)
p1 <- grid_arrange_shared_legend(p_auc_young, p_calibration_young)
p2 <- grid_arrange_shared_legend(p_auc_old, p_calibration_old)
