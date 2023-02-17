
library("ggsci")  # for color palette
library(grid)
library(gridExtra)

path1 <- "/Users/work/IDrive-Sync/Projects/MIMAH/code/AH_code/original models/"
load(paste0(path1, "ROC_original.Rdata"))

roc_list_orig <- roc.list

path2 <- "/Users/work/IDrive-Sync/Projects/MIMAH/code/AH_code/updating models/"
load(paste0(path2, "ROC_updated.Rdata"))

roc_list_update <- roc.list

model <- "meld" # meld, lille, clif

if(model == "meld"){
    # MELD
    model_list <- list("meld_orig" = roc_list_orig$MELD,
                       "meld_update" = roc_list_update$MELD)
    
    P_value_boot <- roc.test(roc_list_orig$MELD, roc_list_update$MELD, method = "bootstrap") # bootstrap p_value
    
}else if(model == "lille"){
    #Lille
    model_list <- list("lille_orig" = roc_list_orig$Lille,
                       "lille_update" = roc_list_update$Lille)
   
    P_value_boot <- roc.test(roc_list_orig$Lille, roc_list_update$Lille, method = "bootstrap")
    
}else{
    # CLIF 
    model_list <- list("clif_orig" = roc_list_orig$`CLIF-C ACLF`,
                       "clif_update" = roc_list_update$`CLIF-C ACLF`)
    
    P_value_boot <- roc.test(roc_list_orig$`CLIF-C ACLF`, roc_list_update$`CLIF-C ACLF`, method = "bootstrap")
    
}

g.list <- ggroc(model_list)

# ROC plot of all models combined
g.list +
    geom_line(lwd = 1) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed") + 
    theme_classic()

# add confidence bands
ci.list <- lapply(model_list, ci.se, specificities = seq(0, 1, l = 25))

dat.ci.list <- lapply(ci.list, function(ciobj) 
    data.frame(x = as.numeric(rownames(ciobj)),
               lower = ciobj[, 1],
               upper = ciobj[, 3]))

df <- plyr::ldply(dat.ci.list, data.frame, .id = "name")

# individual plot
ggroc(model_list) + 
    ggtitle(model) +
    facet_grid(. ~ name) +
    theme_minimal() + 
    geom_ribbon(data = df, aes(x = x, ymin = lower, ymax = upper, fill = name), alpha = 0.3, inherit.aes = F) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed") + 
    #geom_text(x = 0.5, y = 0.5, label = paste0("p.value")) +
    annotate(geom = "text", x = 0.3, y = 0.5, label =  paste0("p.value = ", round(P_value_boot$p.value, 2)), fontsize = 12) +
    labs(x = "Specificity", y = "Sensitivity") + 
    scale_fill_jco() + 
    scale_color_jco() +
    coord_equal() +
    theme_classic() +
    theme(legend.position = "none") 


#### arranged plots
linewidth <- 1.3
text_size <- 4.5

roc1 <- ggroc(model_list, lwd = linewidth) + 
    ggtitle("CLIF-C ACLF") +
    #facet_grid(. ~ name) +
    #theme_minimal() + 
    geom_ribbon(data = df, aes(x = x, ymin = lower, ymax = upper, fill = name), alpha = 0.3, inherit.aes = F) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed") + 
    #geom_text(x = 0.5, y = 0.5, label = paste0("p.value")) +
    annotate(geom = "text", x = 0.3, y = 0.5, label =  paste0("p.value = ", round(P_value_boot$p.value, 2)), 
             size = text_size) +
    labs(x = " ", y = " ") + 
    scale_fill_jco() + 
    scale_color_jco() +
    coord_equal() +
    theme_classic2() +
    theme(legend.position = "none") 

roc2 <- ggroc(model_list, lwd = linewidth) + 
    ggtitle("Lille") +
    #facet_grid(. ~ name) +
    #theme_minimal() + 
    geom_ribbon(data = df, aes(x = x, ymin = lower, ymax = upper, fill = name), alpha = 0.3, inherit.aes = F) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed") + 
    #geom_text(x = 0.5, y = 0.5, label = paste0("p.value")) +
    annotate(geom = "text", x = 0.3, y = 0.5, label =  paste0("p.value = ", round(P_value_boot$p.value, 2)), 
             size = text_size) +
    labs(x = " ", y = " ") + 
    scale_fill_jco() + 
    scale_color_jco() +
    coord_equal() +
    theme_classic2() +
    theme(legend.position = "none") 

roc3 <- ggroc(model_list, lwd = linewidth) + 
    ggtitle("MELD") +
    #facet_grid(. ~ name) +
    #theme_minimal() + 
    geom_ribbon(data = df, aes(x = x, ymin = lower, ymax = upper, fill = name), alpha = 0.3, inherit.aes = F) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed") + 
    #geom_text(x = 0.5, y = 0.5, label = paste0("p.value")) +
    annotate(geom = "text", x = 0.3, y = 0.5, label =  paste0("p.value = ", round(P_value_boot$p.value, 2)), 
             size = text_size) +
    labs(x = " ", y = " ") + 
    scale_fill_jco() + 
    scale_color_jco() +
    coord_equal() +
    theme_classic2() +
    theme(legend.position = "none") 

grid.arrange(arrangeGrob(roc1), 
             arrangeGrob(roc2), 
             arrangeGrob(roc3),
             ncol = 3, left = "Sensitivity", bottom = "Specificity")



