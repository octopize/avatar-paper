# Main libraries
library(reticulate)
library(tidyverse)
library(mice)
library(reshape2)

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
py_run_string("from lsg.security_metrics.avatars_are_first_hit import avatars_are_first_hit")
py_run_string("from lsg.security_metrics.hidden_rate import hidden_rate")

file_data <- "../../datasets/"
file_results <- "../../figure/"

data <- read.csv(paste0(file_data, "AIDS/aids_original_data.csv"), sep = ";", na.strings = c("NA", "", NA))
data <- data[,-1]

save <- FALSE
compute <- FALSE
seq_k <- c(4, 200, 400, 600, 800, 1000)

colors = read.csv("../../color.csv", stringsAsFactors = FALSE)
rownames(colors) <- colors$type

axis_title_size <- 19
axis_text_size <- 19
legend_text_size <- 19
legend_title_size <- 19
axis_text_size_small <- 16

data_imp <- complete(mice(data, m = 1, maxit = 50, seed = 500, printFlag = FALSE), 1)
data$cd496 <- data_imp$cd496

categorical <- c('hemo','homo', 'drugs', 'karnof', 'oprior', 'z30', 'zprior', 'race', 'gender', 'str2', 'strat', 'symptom', 'treat', 'offtrt', 'r','cens', 'arms')
data[categorical] <- lapply(data[categorical], factor)

write_pvalue <- function(val, dec) {
  if (is.na(val)) {
    return(NA)
  } else if (val < 0.001) {
    return("<0.001")
  } else {
    return(formatC(val, format = 'f', digits = dec))
  }
}

# Small k
k <- 4 ; iter <- 0
avatar_tot <- read.table(unz(paste0(file_data, "AIDS/aids_avatarized.zip"), 
                             paste0("aids_avatarized/k", k, "_aids_avatarized.csv")), sep = ",", header = TRUE)
avatar_small <- avatar_tot[avatar_tot$iter_k == iter,]

# Imputation
avatar_small_imp <- complete(mice(avatar_small, m = 1, maxit = 50, seed = 500, printFlag = FALSE), 1)
avatar_small$cd496 <- avatar_small_imp$cd496
avatar_small[categorical] <- lapply(avatar_small[categorical], factor)

# High k = N/2
k <- 1166 ; iter <- 1
avatar_tot <- read.table(unz(paste0(file_data, "AIDS/aids_avatarized.zip"), 
                             paste0("aids_avatarized/k", k, "_aids_avatarized.csv")), sep = ",", header = TRUE)
avatar_med <- avatar_tot[avatar_tot$iter_k == iter,]

# Imputation
avatar_med_imp <- complete(mice(avatar_med, m = 1, maxit = 50, seed = 500, printFlag = FALSE), 1)
avatar_med$cd496 <- avatar_med_imp$cd496
avatar_med[categorical] <- lapply(avatar_med[categorical], factor)

# binding data, avatar small and avatar high
data_small_med <- rbind(data[,1:26], avatar_small[,1:26], avatar_med[,1:26])
row.names(data_small_med) <- 1:nrow(data_small_med)

# perform projection
famd <- FAMD(data_small_med, ncp = 5, graph = FALSE, ind.sup = (nrow(data_small_med)/3 +1):nrow(data_small_med))

res_ind_2D <- as.data.frame(famd$ind.sup$coord)[1:nrow(data_small_med)/3, c('Dim.1', 'Dim.2')]
res_ind_sup <- as.data.frame(famd$ind.sup$coord)[(nrow(data_small_med)/3):nrow(famd$ind.sup$coord),]
res_ind_2D['type'] <- 'k = 4'
res_ind_sup['type'] <- 'k = 1166'
res_ind_2D <- rbind(res_ind_2D, res_ind_sup[, c('Dim.1', 'Dim.2', 'type')])
res_ind_2D$type <- factor(res_ind_2D$type, levels = c('k = 4', 'k = 1166'))

set.seed(42)
rows <- sample(nrow(res_ind_2D))
res_ind_2D <- res_ind_2D[rows, ]

plotCa <- ggplot(res_ind_2D, aes(x = Dim.1, y = Dim.2, fill = type)) +
    # point
    geom_point(size = 3, shape = 21, alpha = 1) +
    # fill by k value
    aes(fill = factor(type)) + 
    scale_fill_manual(values = c('k = 4' = colors["k_small", "color"], 'k = 1166' = colors["k_high", "color"])) +
    # theme and details
    xlab(paste0("Dim. 1 (", formatC(famd$eig[1,2], format = "f", digits = 2) ,"%)")) +
    ylab(paste0("Dim. 2 (", formatC(famd$eig[2,2], format = "f", digits = 2) ,"%)")) +
    labs(fill = "") +
    theme_bw() +
    theme(legend.position = c(0.9, 0.12),
          legend.text = element_text(size = legend_text_size, color = "black", family = "sans"),
          legend.key.size = unit(0.8, "cm"),
          axis.text = element_text(size = axis_text_size_small, color = "black", family = "sans"),
          axis.title = element_text(size = axis_title_size, color = "black", family = "sans"),
          legend.background = element_rect(fill = "white", linetype = "solid"))

data[] <- lapply(data, function(x){if(is.factor(x)) as.numeric(as.character(x)) else x})

n_arms <- as.vector(table(data$arms))
cens_arms <- paste0(as.vector(table(data$cens, data$arms)[2,]), " (", 
                    as.vector(round(prop.table(table(data$cens, data$arms), 2)[2, ]*100, 1)), ")")

summary_cox <- summary(coxph(Surv(time = days/7, event = cens)~factor(arms), data = data))

HR_arms <- c("ref", round(as.vector(summary_cox$coefficients[,2]), 2))
CI_low_arms <- c(NA, round(as.vector(summary_cox$conf.int[,3]), 2))
CI_high_arms <- c(NA, round(as.vector(summary_cox$conf.int[,4]), 2))
p_arms <- c(NA, sapply(as.vector(summary_cox$coefficients[,5]), write_pvalue, dec = 3))

pr <- py$Projection()
coord_original <- pr$fit_transform(data, nf = 5L)[[1]]

res_df <- read.csv("../../datasets/results_df/hazard_ratio.csv")

res_df["names"] <- paste(res_df$k, res_df$iter, sep = "_")

df_res_HR <- res_df[, c("k", "iter", "res_HR1", "res_HR2", "res_HR3")]
df_res_HR_CI_low <- res_df[, c("k", "iter", "res_HR_CI_low1", "res_HR_CI_low2", "res_HR_CI_low3" )]
df_res_HR_CI_high <- res_df[, c("k", "iter", "res_HR_CI_high1", "res_HR_CI_high2", "res_HR_CI_high3")]

# compute confidence intervals
mean_hr1 <- res_df %>%
  group_by(k) %>%
  summarise(mean = mean(res_HR1),
            std = sqrt(var(res_HR1)),
            lower = mean(res_HR1) - qnorm(.975)*std/sqrt(n()),
            upper = mean(res_HR1) + qnorm(.975)*std/sqrt(n()))
mean_hr1$inf <- mean_hr1$mean - mean_hr1$std
mean_hr1$sup <- mean_hr1$mean + mean_hr1$std

poly_df <- rbind(setNames(mean_hr1[,c("k", "inf")], c('x','y')),
                 setNames(mean_hr1[order(nrow(mean_hr1):1), c("k", "sup")], c('x','y')))

g1 <- ggplot(res_df) +
    # red line as a ploygone
    stat_smooth(aes(k, res_HR1), method = "lm", se=FALSE, 
                formula = y ~ poly(x, 8, raw=TRUE), colour="red") +
    # add upper and lower bounds
    stat_smooth(data = mean_hr1, aes(k, mean - std), col = "#696969", lty = 2, se = FALSE, 
                method = 'gam', formula = y ~ s(x, bs = "cs")) +
    stat_smooth(data = mean_hr1, aes(k, mean + std), col = "#696969", lty = 2, se = FALSE,
                method = 'gam', formula = y ~ s(x, bs = "cs")) 

gg1 <- ggplot_build(g1)

# get stat smooth (2,3) coordinate in order to draw the polygone between the smooth curve
df2 <- data.frame(x = gg1$data[[2]]$x, ymin = gg1$data[[2]]$y, ymax = gg1$data[[3]]$y) 
df3 <- data.frame(x = gg1$data[[2]]$x, ymin = 0.39, ymax = 0.63) 

# draw final figure 
# we must redo all the code lines to manage the layer order (point on the grey zone)

plotCc <- ggplot() +
    # original zone 
    geom_ribbon(data = df3, aes(x = x, ymin = ymin, ymax = ymax),
                fill = colors["original", "color"], alpha = 0.07) + 
       
    # avatar zone 
    geom_ribbon(data = df2, aes(x = x, ymin = ymin, ymax = ymax),
                fill = colors["avatar", "color"], alpha = 0.12) + 
    theme_minimal() +

    # add original mean
    geom_segment(aes(x = 0, y = 0.49, xend = 1175, yend = 0.49, colour = colors["original", "color"]), lwd = 1, lty = 1, alpha = 0.8) +

    # add mean point by k
    geom_point(data = mean_hr1, aes(x = k, y = mean), 
               colour = "white", alpha = 0.9, size = 2.5) +
    geom_point(data = mean_hr1, aes(x = k, y = mean), 
               colour = colors["avatar", "color"], alpha = 0.9, size = 2) +
    # avatar mean smooth
    stat_smooth(data = res_df, aes(k, res_HR1), method = "lm", se = FALSE, 
                formula = y ~ poly(x, 8, raw = TRUE), colour = "#35b696", size = 1.5) +
    stat_smooth(data = res_df, aes(k, res_HR1, colour = colors["avatar", "color"]), method = "lm", se = FALSE, 
                formula = y ~ poly(x, 8, raw = TRUE), size = 1) +
        
    # labs theme and axis
    ylab("Hazard ratio between\narms 0 and arms 1") +
    xlab("Number of nearest neighbors (k)") +
    ylim(0.35, 1) +
    #xlim(4, )
    scale_x_continuous(breaks = seq_k)  +
    coord_cartesian(xlim = c(0, 1075)) +
    scale_color_identity(guide = 'legend', name = NULL, breaks= c(colors["original", "color"], colors["avatar", "color"]), labels= c('Original', 'Avatar')) +
    theme(legend.position = c(0.13, 0.9),
          axis.title.y = element_text(angle = 90, vjust = 1, size = axis_title_size, family = "sans"),
          axis.text = element_text(size = axis_text_size_small, colour = "black"), 
          axis.title.x = element_text(size = axis_title_size, family = "sans"),
          legend.text = element_text(size = legend_text_size, color = "black", family = "sans"),
          legend.key.size = unit(1, 'cm'),
          axis.line = element_line(colour = "black", size = 0.5, linetype = "solid", 
                                   arrow = arrow(type = 'closed', length = unit(5,'pt'))))


df_local_cloaking_k100 <- data.frame(row.names = 1:nrow(data))
pr <- py$Projection()
hidden_rate_list <- c()


coord_original <- pr$fit_transform(data, nf = 5L)[[1]]
for (k in seq_k) {
    iter <- 0
    avatar_tot <- read.table(unz(paste0(file_data, "AIDS/aids_avatarized.zip"), 
                                 paste0("aids_avatarized/k", k, "_aids_avatarized.csv")), sep = ",", header = TRUE)
    avatar <- avatar_tot[avatar_tot$iter_k == iter,]
    avatar_imp <- complete(mice(avatar, m = 1, maxit = 50, seed = 500, printFlag = FALSE), 1)
    avatar$cd496 <- avatar_imp$cd496
    coord_avatar <- pr$transform(avatar[, 1:(ncol(avatar) - 3)])
    
    # local cloaking
    distances <- py$record_to_avatar_distance(coord_original, coord_avatar)
    local_cloaking <- py$local_cloaking(coord_original, coord_avatar, distances)
    df_local_cloaking_k100[paste0("local_cloaking_", k)] <- local_cloaking$hit_counts[, 1]
}

df_melt <- melt(df_local_cloaking_k100)

plotCe <- ggplot(df_melt, aes(x = variable, y = value)) + 
    geom_violin(trim = TRUE, alpha = 0.12, fill = colors['avatar', 'color']) +
    geom_boxplot(width = 0.4, fill = colors['avatar', 'color'], outlier.size = 0.1,alpha = 0.8) +
    scale_x_discrete(labels = as.character(seq_k)) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), breaks = c(1,5,10,25,50,100,300, 600, 1000, 2000)) +
    xlab("Number of nearest neighbors (k)") +
    ylab("Local cloaking") +
    theme_minimal() +  
    theme(legend.position = "none",
          axis.title.y = element_text(angle = 90, vjust = 1, size = axis_title_size),
          axis.text = element_text(size = axis_text_size_small, colour = "black"), 
          axis.title.x = element_text(size = axis_title_size),
          axis.line = element_line(colour = "black", size = 0.5, linetype = "solid", 
                                   arrow = arrow(type = 'closed', length = unit(5, 'pt'))))+
    coord_cartesian(ylim = c(-1, 2000))