library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(forcats)
library(purrr)

### data 
path_data <- "/Users/work/IDrive-Sync/Projects/MIMAH/code/AH_code/updating models"
load(paste0(path_data, "/recalibrated_models_default.Rdata"))


# subset data
data <- test.data %>% 
    select(Subject, meld.surv.updated, lille.surv.updated, clif.surv.updated, D90_surv)

#############################################   
################ Correlation  ###############
############################################# 

#ggplot(data, aes(x = meld.surv.updated, y = lille.surv.updated)) +
#    geom_point() +
#    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
#    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) 

ggscatter(data, x = "meld.surv.updated", y = "lille.surv.updated",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", cor.coef.coord = c(0.10, 0.9)) + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
    labs(x = "MELD survival probability ", y = "Lille survival probability")
coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_classic()

ggscatter(data, x = "meld.surv.updated", y = "clif.surv.updated",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", cor.coef.coord = c(0.10, 0.9)) + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
    labs(x = "MELD survival probability", y = "CLIF-C ACLF survival probability ")
coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_classic()

ggscatter(data, x = "lille.surv.updated", y = "clif.surv.updated",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", cor.coef.coord = c(0.10, 0.9)) + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
    labs(x = "Lille survival probability ", y = "CLIF-C ACLF survival probability ")
coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_classic()

#############################################   
################ Correlation  ###############
#############################################
# logn format
data_long <- gather(data, model, probs, meld.surv.updated:clif.surv.updated, factor_key = TRUE)

data_long2 <- data_long %>% 
    group_by(Subject) %>% 
    mutate(dev = sd(probs), # sd 
           mad_dev = mad(probs), #median absolute deviation
           max_min_dev = max(probs) - min(probs)) # range 

# plot
data_long2 %>% 
    mutate(Subject = fct_reorder(as.factor(Subject), desc(max_min_dev))) %>%
    ggplot(., aes(x = max_min_dev, y = reorder(Subject, max_min_dev))) +
    geom_point(alpha = 0.5) +
    labs(x = "Probabilities range", y = "Subject") +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())


# overall survival rate
data_long2 %>% ungroup() %>% summarise(mean(D90_surv))

# fun 
cut_off_surv_rate <- function(data, thres, var_name){
    # calculates mean of `var_name` and sample size for chosen `thres` value
    
    df <- data %>% 
        filter(max_min_dev > thres) %>% 
        ungroup() %>%
        summarise(sum_var = mean({{ var_name }}, na.rm = T))
    
    df_n <- data %>% 
        filter(max_min_dev > thres) %>% 
        ungroup() %>%
        summarise(n = length(unique(Subject)))
    
    out <- data.frame(surv_rate = df$sum_var, thres = thres, n = df_n$n)
    
    return(out)
}


thresholds <- seq(0, 0.5, by = 0.01)
res <- map_df(thresholds, cut_off_surv_rate, data = data_long2, var_name = D90_surv)

# plot
ggplot(res, aes(y = surv_rate, x = thres, label = n)) +
    geom_point() +
    geom_line() + 
    geom_text(check_overlap = TRUE, nudge_y = 0.02) +
    labs(x = "Probabilities range", y = "Survival Rate") +
    theme_classic()

# fun 
sum_fun <- function(df1, df2, thres){
    # choses ids based the `thres` value 
    
    df_prelim <- df1 %>% filter(max_min_dev > thres) %>% arrange(Subject)
    
    df_final <- df2 %>% filter(Subject %in% df_prelim$Subject)
    
    df_final$threshold <- thres
    
    return(df_final)
}

res <- map_df(thresholds, sum_fun, df1 = data_long2, df2 = test.data)

vars <- c("Subject", "Bilirubin.mg.dl", "Bilirubin.day.7", "delta.bili", "INR", 
          "Creatinine.mg.dl", "Albumin", "WBC", "INR", "protime", "HE",
          "MAP", "Sodium", "CLIF.OF", 
          "kidney.score", "liver.score", "brain.score", "coag.score", "circ.score",
          "Gender.x", "Age.at.randomisation..calc..x")

sum_df <- res %>% 
    group_by(threshold) %>% 
    select(all_of(vars)) %>% 
    summarise_all(mean) %>%
    gather(., variable, value, Bilirubin.mg.dl:Age.at.randomisation..calc..x, factor_key=TRUE)

vars_to_keep <- c("CLIF.OF")

sum_df %>% filter(variable %in% vars_to_keep) %>%
    ggplot(., aes(x = threshold, y = value)) +
    geom_point() +
    geom_line() +
    facet_wrap(. ~ variable, scales = "free_y") +
    labs(x = "Probabilities range") +
    theme_bw()


vars_to_remove <- c("Gender.x", "Age.at.randomisation..calc..x", 
                    "liver.score", "kidney.score", "brain.score", "coag.score",
                    "circ.score", 
                    "protime", "delta.bili")

sum_df %>% filter(!(variable %in% vars_to_remove)) %>%
    ggplot(., aes(x = threshold, y = value)) +
    geom_point() +
    geom_line() +
    facet_wrap(. ~ variable, scales = "free_y") +
    labs(x = "Probabilities range") +
    theme_bw()


#################
#df <- test.data %>% 
#    mutate(Group = ifelse(Subject %in% data_long3$Subject, 1, 0)) %>% 
#    select(Group, vars) %>%
#    group_by(Group) 

# chisq.test
#https://data-flair.training/blogs/chi-square-test-in-r/#:~:text=Chi%2DSquare%20test%20in%20R%20is%20a%20statistical%20method%20which,Green%2C%20Yes%2FNo%20etc.
#test <- chisq.test(table(df$CLIF.OF, df$Group), simulate.p.value = TRUE)
#test

#stacked_1 <- df %>% ungroup(.) %>% filter(Group == 1) %>% select(-Group) %>% stack(.) %>% filter(ind != "Group")
#stacked_1$Group <- "1"
#stacked_2 <- df %>% ungroup(.) %>% filter(Group == 0) %>% select(-Group) %>% stack(.) %>% filter(ind != "Group")
#stacked_2$Group <- "0"

#stacked <- rbind(stacked_1, stacked_2)

# stacked %>% filter(ind != "Subject") %>%  
#     ggplot(., aes(x = ind, y = values, fill = factor(Group))) +
#     geom_boxplot() +
#     facet_wrap(ind ~ ., scales = "free") +
#     stat_pwc(method = "t.test", p.adjust.method = "bonferroni")
# 
# stacked %>% filter(ind != "Subject") %>%  
#     ggplot(., aes(x = ind, y = values, fill = factor(Group))) +
#     geom_boxplot() +
#     geom_pwc(aes(group = Group), tip.length = 0,
#     method = "t_test", p.adjust.method = "none", p.adjust.by = "group",
#     hide.ns = TRUE) #+
#     #facet_wrap(ind ~ ., scales = "free") 
#     
# comp_means <- compare_means(c(Bilirubin.mg.dl, Creatinine.mg.dl, HE, INR, MAP, 
#                 CLIF.OF) ~ Group, data = df, 
#                 method = "t.test",  p.adjust.method = "holm")
# 
# 
# adj_values <- p.adjust(comp_means$p.adj, method = 'holm')
# comp_means$adj_values = adj_values
# comp_means


# # calibration
# # funs
# path_funs <- "/Users/work/IDrive-Sync/Projects/MIMAH/code/funs"
# source(paste0(path_funs, "/calibration_fun.R"))
# 
# #from long to wide 
# data_wide <- spread(data_long3, model, probs)
# str(data_wide)
# 
# # MELD 
# cal_MELD.surv <- calibration(data_wide$meld.surv.updated, y = data_wide$D90_surv)
# cal_MELD.surv$Score <- "MELD"
# 
# # Lille
# cal_Lille <- calibration(data_wide$lille.surv.updated, y = data_wide$D90_surv)
# cal_Lille$Score <- "Lille"
# 
# # CLIF-C ACLF
# cal_CLIF <- calibration(data_wide$clif.surv.updated, y = data_wide$D90_surv)
# cal_CLIF$Score <- "CLIF-C ACLF"
# 
# # combine dfs
# df_cal <- rbind(cal_MELD.surv, cal_Lille, cal_CLIF)
# 
# # plot with ribbon 
# df_cal %>%
#     ggplot(., aes(x = pred, y = obs, col = Score)) +
#     geom_line(linewidth = 1)  + 
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#     geom_ribbon(aes(ymin = lower, ymax = upper, fill = Score, linetype = NA),  
#                 alpha = 0.3, show.legend = F) + 
#     #scale_fill_manual("", values = col) + 
#     #scale_color_manual(name = "Score", values = col) + 
#     facet_grid(. ~ Score) +
#     coord_equal() +
#     xlim(0, 1) + 
#     ylim(0, 1) + 
#     ylab("Observed survival proportion") + 
#     xlab("Predicted survival probability") + 
#     theme_classic() 
# 
# library(gmish)
# ici(data_wide$meld.surv.updated, data_wide$D90_surv)
# ici(data_wide$lille.surv.updated, data_wide$D90_surv)
# ici(data_wide$clif.surv.updated, data_wide$D90_surv)
# 
# mean(abs(cal_MELD.surv$pred - cal_MELD.surv$obs))
# 
# 
# mean((data_wide$meld.surv.updated - data_wide$D90_surv)^2)
# mean((data_wide$lille.surv.updated - data_wide$D90_surv)^2)
# mean((data_wide$clif.surv.updated - data_wide$D90_surv)^2)
# 
# #library(DescTools)
# #BrierScore(data_wide$D90_surv, data_wide$meld.surv.updated)
# #BrierScore(data_wide$D90_surv, data_wide$lille.surv.updated)
# #BrierScore(data_wide$D90_surv, data_wide$clif.surv.updated)
# 
# library(pROC)
# # MELD (survival function 1 and 2 do not matter here)
# roc_meld <- roc(data_wide$D90_surv, data_wide$meld.surv.updated)
# roc_meld
# # Lille
# roc_lille <- roc(data_wide$D90_surv, data_wide$lille.surv.updated)
# roc_lille
# # CLIF-C ACLF
# roc_clif <- roc(data_wide$D90_surv, data_wide$clif.surv.updated)
# roc_clif