# Main libraries
library(tidyverse)
library(mice)
library(reshape2)

# Library to use Python with R
library(reticulate)

# Multidimensionnal reduction method libary
library(FactoMineR)

# Survival analysis libraries
library(survival)
library(survminer)

# import specific function to compute metrics
py_run_string("import sys")
py_run_string("sys.path.append('../..')")

py_run_string("from lsg.security_metrics.security_metrics_class import Security_metrics ")
py_run_string("from lsg.dimension.projection import Projection")
py_run_string("from lsg.security_metrics.record_to_avatar_distance import record_to_avatar_distance")
py_run_string("from lsg.security_metrics.local_cloaking import local_cloaking")

file_data <- "datasets/"
data <- read.csv(paste0("../../datasets/AIDS/aids_original_data.csv"), sep = ";", na.strings = c("NA", "", NA))
data <- data[, -1]

# Avatar data with k=20
avatar <- read.csv(paste0("../../datasets/AIDS/aids_avatarized_base_k20_nf5.csv"), sep = ",", na.strings = c("NA", "", NA))
avatar$arms <- as.factor(avatar$arms)

colors <- read.csv("../../color.csv", stringsAsFactors = FALSE)
rownames(colors) <- colors$type

axis_title_size <- 19
axis_text_size <- 19
legend_text_size <- 19
legend_title_size <- 19
axis_text_size_small <- 16

# imput missing data to project the individuals in the multidimensionnal space
data_imp <- complete(mice(data, m = 1, maxit = 50, seed = 500, printFlag = FALSE), 1)
avatar_imp <- complete(mice(avatar, m = 1, maxit = 50, seed = 500, printFlag = FALSE), 1)

data$cd496 <- data_imp$cd496
avatar$cd496 <- avatar_imp$cd496

data_tot <- rbind(data, avatar)
data_tot_num <- data_tot
categorical <- c("hemo", "homo", "drugs", "karnof", "oprior", "z30", "zprior", "race", 
                 "gender", "str2", "strat", "symptom", "treat", "offtrt", "r", "cens", "arms")
data_tot[categorical] <- lapply(data_tot[categorical], factor)

pca <- FAMD(data_tot, ncp = 5, graph = FALSE, ind.sup = (nrow(data_tot)/2 + 1):nrow(data_tot))
res_ind <- as.data.frame(pca$ind$coord)
res_ind_sup <- as.data.frame(pca$ind.sup$coord)

res_ind_2D <- res_ind[, 1:2]
res_ind_2D["type"] <- "Original"
res_ind_sup["type"] <- "Avatar"
res_ind_2D <- rbind(res_ind_2D, res_ind_sup[, c("Dim.1", "Dim.2", "type")])

set.seed(42)
rows <- sample(nrow(res_ind_2D))
res_ind_2D <- res_ind_2D[rows,]

plotAa <- ggplot(res_ind_2D, aes(x = Dim.1, y = Dim.2, fill = factor(type))) + 
    # point
    geom_point(size = 3, shape = 21, alpha = 1) + 
    # fill by type 
    aes(fill = factor(type)) + 
    scale_fill_manual(values = c(colors["avatar", "color"], colors["original", "color"])) +
    # theme and details
    xlab(paste0("Dim. 1 (", formatC(pca$eig[1, 2], format = "f", digits = 2), "%)")) +
    ylab(paste0("Dim. 2 (", formatC(pca$eig[2, 2], format = "f", digits = 2), "%)"))+ 
    labs(fill = "") +
    theme_bw() +
    theme(legend.position = c(0.86, 0.15), 
          legend.text = element_text(size = legend_text_size, color = "black", family = "sans"), 
          axis.text = element_text(size = axis_text_size_small, color = "black", family = "sans"), 
          axis.title = element_text(size = axis_title_size, color = "black", family = "sans"), 
          legend.background = element_rect(fill = "white", linetype = "solid"))


data$arms <- as.factor(data$arms)
summary_cox <- summary(coxph(Surv(time = days/7, event = cens) ~ arms, data = data))
res_ori <- cbind(summary_cox$coefficients[, c("exp(coef)", "Pr(>|z|)")], summary_cox$conf.int[, c("lower .95", "upper .95")])

data_arms01 <- data[data$arms %in% c(0, 1),]
data_arms01$arms <- as.factor(as.character(data_arms01$arms))
summary_cox <- summary(coxph(Surv(time = days/7, event = cens) ~ arms, data = data_arms01))
res_original_01 <- cbind(summary_cox$coefficients, summary_cox$conf.int)

summary_cox <- summary(coxph(Surv(time = days/7, event = cens) ~ arms, data = avatar))
res_avat <- cbind(summary_cox$coefficients[, c("exp(coef)", "Pr(>|z|)")], summary_cox$conf.int[, c("lower .95", "upper .95")])

avatar_arms01 <- avatar[avatar$arms %in% c(0, 1),]
avatar_arms01$arms <- as.factor(as.character(avatar_arms01$arms))
summary_cox <- summary(coxph(Surv(time = days/7, event = cens) ~ arms, data = avatar_arms01))
res_avatar_01 <- cbind(summary_cox$coefficients, summary_cox$conf.int)

data_arms01$type <- "Original"
avatar_arms01$type <- "Avatar"
data_concat_01 <- rbind(data_arms01, avatar_arms01)
data_concat_01$type <- factor(data_concat_01$type, levels = c("Original", "Avatar"))

# Survival curve calculation
surv <- survfit(Surv(time = days/7, event = cens) ~ factor(arms) + type, data = data_concat_01)

# Display time: about 15 seconds

color <- c("#7155c8", "#c85573")

options(repr.plot.width = 10, repr.plot.height = 7)
ggsurv <- ggsurvplot(surv, data = data_concat_01, conf.int = FALSE, 
                     risk.table = FALSE, 
                     legend.title = "Treatments",
                     censor.cex = 100, 
                     linetype = c("type"), 
                     size = 0.7, xlim = c(0, 155), ylim = c(0.5, 1), break.x.by = 40, 
                     ggtheme = theme_minimal(), 
                     xlab = "Time (week)", 
                     ylab = "Proportion not reaching\nprimary end point",
                     censor.shape = "", 
                     legend.labs = c("Arms 0", "Arms 0", "Arms 1", "Arms 1"), 
                     palette = color)

plotAc <- ggsurv$plot <- ggsurv$plot +
    theme_minimal() +
    scale_linetype_discrete(name = c("Data source"), labels = c("Original", "Avatar")) +
    geom_label(aes(x = 15.5, y = 0.623, label = "Hazard ratio"), size = 6.5, family = "sans", fontface = "bold", label.size = NA) +
    geom_label(aes(x = 66, y = 0.574, label = paste0("Original: HR [CI: 95%] = ", formatC(res_original_01[2], format = "f", digits = 2),
                                                      " [", formatC(res_original_01[8], format = "f", digits = 2),
                                                      "-", formatC(res_original_01[9], format = "f", digits = 2),
                                                      "] ; p-value = ", formatC(res_original_01[5], format = "e", digits = 2), 
                                                      "\nAvatar: HR [CI: 95%] = ", formatC(res_avatar_01[2], format = "f", digits = 2),
                                                      " [", formatC(res_avatar_01[8], format = "f", digits = 2),
                                                      "-", formatC(res_avatar_01[9], format = "f", digits = 2),
                                                      "] ; p-value = ", formatC(res_avatar_01[5], format = "e", digits = 2), "  ")), 
               size = 6.2, family = "sans", label.size = NA) +

    theme(legend.position = c(0.14, 0.5), 
          legend.background = element_rect(fill = "white", size = 0.5, linetype = "blank"), 
          legend.text = element_text(size = legend_text_size,  color = "black", family = "sans"), 
          legend.title = element_text(size = legend_title_size,  color = "black", family = "sans", face = "bold"),
          axis.text = element_text(size = axis_text_size_small, color = "black", family = "sans"), 
          axis.title.y = element_text(vjust = 2), 
          axis.title = element_text(size = axis_title_size, color = "black", family = "sans"), 
          axis.line = element_line(colour = "black", size = 0.5, linetype = "solid", arrow = arrow(type = "closed", length = unit(5, "pt"))),
          text = element_text(size = 14))


data$type <- "Original"
avatar$type <- "Avatar"
data_concat <- rbind(data, avatar)
data_concat$type <- factor(data_concat$type, levels = c("Original", "Avatar"))

# Survival curve calculation
surv <- survfit(Surv(time = days/7, event = cens) ~ factor(arms) + type, data = data_concat)

color <- c("#7155c8", "#c85573", "#95aa5b", "#c37e47")

ggsurv <- ggsurvplot(surv, 
                     data = data_concat, 
                     conf.int = FALSE, 
                     risk.table = FALSE, 
                     legend.title = "Treatments", 
                     legend = c(0.2, 0.3), 
                     censor.cex = 100, 
                     linetype = c("type"), 
                     size = 0.8, 
                     xlim = c(0, 155), 
                     ylim = c(0.5, 1),
                     break.x.by = 40, 
                     ggtheme = theme_minimal(), 
                     xlab = "Time (week)", 
                     ylab = "Proportion not reaching\nprimary end point", 
                     censor.shape = "", 
                     legend.labs = c("Arms 0", "Arms 0", "Arms 1", "Arms 1", "Arms 2", "Arms 2", "Arms 3", "Arms 3"), 
                     palette = color)

ggsurv$plot <- ggsurv$plot + 
    theme_minimal() + 
    theme(legend.position = c(0.105, 0.5), 
          legend.background = element_rect(fill = "white", size = 0.5, linetype = "blank"),
          legend.text = element_text(size = legend_text_size, color = "black", family = "sans"), 
          axis.text = element_text(size = axis_text_size, color = "black",  family = "sans"), 
          axis.title = element_text(size = axis_title_size, color = "black", family = "sans", face = "bold"), 
          axis.title.y = element_text(vjust = 2), 
          legend.title = element_text(size = legend_title_size, face = "bold", family = "sans"), 
          axis.line = element_line(colour = "black", size = 0.5, linetype = "solid", arrow = arrow(type = "closed", length = unit(5, "pt")))) + 
    scale_linetype_discrete(name = "Data source")


rownames(res_ori) <- paste0(rownames(res_ori), "_original")
rownames(res_avat) <- paste0(rownames(res_avat), "_avatar")
res_concat <- as.data.frame(rbind(res_ori, res_avat))
res_concat["Hazard Ratio"] <- formatC(res_concat[, 1], format = "f", digits = 2)
res_concat["95% Confidence Interval"] <- paste0(formatC(res_concat[,3], format = "f", digits = 2), "-", formatC(res_concat[, 4], format = "f", digits = 2))
res_concat["p-value"] <- formatC(res_concat[,2], format = "e", digits = 2)
res_concat <- res_concat[c(1, 4, 2, 5, 3, 6), 5:7]
row.names(res_concat) <- c("Original arms 1", "Avatar arms 1", "Original arms 2", "Avatar arms 2", "Original arms 3", "Avatar arms 3")
res_concat

metrics <- py$Security_metrics()
metrics$fit(data, avatar, nf = 5L)
df_local_cloaking <- data.frame(metrics$local_cloaking["hit_counts"])

color_vector <- c(rep("#c14f4f", 2), rep(colors["avatar", "color"], 98))
intersect_median_density <- 44


plotBa <- ggplot(df_local_cloaking, aes(hit_counts.1)) + 
    # add histogram
    geom_histogram(bins = 100, fill = color_vector, color = "#1D1D1B") + 
    #  add density
    geom_density(aes(y = ..count..), lwd = 1.5, color = "#1D1D1B") + 
    # add median
    geom_segment(aes(x = median(hit_counts.1), y = 0, xend = median(hit_counts.1), yend = 45), color = "#1D1D1B", size = 1.5) +
    # theme, axis, leged
    theme_minimal() + 
    xlab("Local cloaking") + 
    ylab("Number of individuals") + 
    theme(axis.title = element_text(size = axis_title_size, family = "sans", color = "black"), 
          axis.text = element_text(size = axis_text_size_small, family = "sans", color = "black"), 
          axis.line = element_line(colour = "black", size = 0.5, linetype = "solid", arrow = arrow(type = "closed", length = unit(5, "pt")))) + 
    xlim(c(-1, 100))  +

    # add arrow and text 
    ## segment median 
    geom_segment(aes(x = median(hit_counts.1),
                     y = 0, xend = median(hit_counts.1),
                     yend = 45), 
                 color = '#1D1D1B', size = 2) +
    ## arrow median
    geom_segment(aes(x = median(hit_counts.1) + 25,
                     y = intersect_median_density + 12, 
                     xend = median(hit_counts.1)+1,
                     yend = 45+2), 
                 color = '#1D1D1B', size = 0.8, 
                 arrow = arrow(type='closed', length = unit(10,'pt'))) +
    ## text median
    geom_label(aes(x =  median(hit_counts.1) + 38, 
                   y =  intersect_median_density + 15,
                   label =paste0("Median = ", median(hit_counts.1))), 
               size = 8, 
               family = "sans", 
               label.size = NA) + 
    ## arrow hidden rate
    geom_segment(aes(x = median(hit_counts.1) + 25, 
                     y = intersect_median_density + 38, 
                     xend = 1, yend = sum(df_local_cloaking[,1] == 0) ), 
                 color = '#1D1D1B', size = 0.8, 
                 arrow = arrow(type='closed', 
                               length = unit(10,'pt'))) +
    ## text hidden rate
    geom_label(aes(x =  median(hit_counts.1) + 43.5, 
                   y = intersect_median_density + 41,
                   label = paste0("Hidden rate = ", round(metrics$hidden_rate, 0), "%")), 
               size = 8, 
               family = "sans", 
               label.size = NA)

avatar_tot <- read.csv('../../datasets/AIDS/aids_avatarized_25time_k20_nf_5.csv')
avatar_tot[categorical] <- lapply(avatar_tot[categorical], factor)
data[categorical] <- lapply(data[categorical], factor)
pr <- py$Projection()
coord_original <- pr$fit_transform(data[,-ncol(data)], nf = 5L)[[1]]

# Computation time: about 2 minutes
defaultW <- getOption("warn")
options(warn = -1)
df_local_cloaking_k2 <- data.frame(row.names = 1:nrow(data))
for (k in unique(avatar_tot$iter)) {
    avatar <- avatar_tot[avatar_tot$iter == k,]
    avatar_imp <- complete(mice(avatar, m = 1, maxit = 50, seed = 500, printFlag = FALSE), 1)
    avatar$cd496 <- avatar_imp$cd496
    
    coord_avatar <- pr$transform(avatar[, 1:(ncol(avatar) - 1)])
    distances <- py$record_to_avatar_distance(coord_original, coord_avatar)
    
    local_cloaking <- py$local_cloaking(coord_original, coord_avatar, distances)
    df_local_cloaking_k2[paste0("local_cloaking_", k, "_", 20)] <- local_cloaking$hit_counts[, 1]
}   
options(warn = defaultW)

df_local_cloaking_k2_LC0 <- as.data.frame(table(rowSums(df_local_cloaking_k2 == 0)))
df_local_cloaking_k2_LC0 <- tibble(df_local_cloaking_k2_LC0 %>% 
    mutate(Var1 = as.character(fct_collapse(Var1, "5" = as.character(5:max(as.numeric(df_local_cloaking_k2_LC0$Var1) - 1))))) %>%
    group_by(Var1) %>%
    mutate(Freq = sum(Freq)) %>%
    distinct())
df_local_cloaking_k2_LC0["Freq2"] <- paste0(formatC(df_local_cloaking_k2_LC0$Freq/sum(df_local_cloaking_k2_LC0$Freq)*100, format = "f", digits = 1), "%")
df_local_cloaking_k2_LC0[df_local_cloaking_k2_LC0$Var1 == 5, "Var1"] <- "\U2265 5"
df_local_cloaking_k2_LC0$Var1 <- factor(df_local_cloaking_k2_LC0$Var1, levels = c(as.numeric(0:4), "\U2265 5"))

plotBc <- ggplot(df_local_cloaking_k2_LC0, aes(x = as.factor(Var1), y = Freq)) + 
    geom_bar(stat = "identity", fill = colors["original", "color"], colour = "black", show.legend = FALSE) + 
    geom_text(aes(label = Freq2), stat = "identity", vjust = -0.5, size = 6) + 
    ylab("Number of individuals") + 
    xlab("Number of local cloaking at zero over 25 avatarizations") + 
    theme_minimal() + 
    theme(axis.title = element_text(size = axis_title_size, family = "sans", color = "black"), 
          axis.text = element_text(size = axis_text_size_small, family = "sans", color = "black"), 
          axis.line = element_line(colour = "black", size = 0.5, 
                                   linetype = "solid", 
                                   arrow = arrow(type = "closed", 
                                                 length = unit(5, "pt"))))