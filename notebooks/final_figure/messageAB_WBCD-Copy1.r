library(librarian)

# Main libraries
librarian::shelf( mice, reshape2)

# Library to use Python with R
librarian::shelf(reticulate,
                # Multidimensional reduction method libary
                FactoMineR,
                # Display a corrplot
                corrplot, 
                 # data manipulation
                 tibble,
                 dplyr
                )

# import specific function to compute metrics 
py_run_string("import sys")
py_run_string("sys.path.append('../..')")
py_run_string("from lsg.security_metrics.security_metrics_class import Security_metrics ")
py_run_string("from lsg.security_metrics.record_to_avatar_distance import record_to_avatar_distance")
py_run_string("from lsg.security_metrics.local_cloaking import local_cloaking")

data <- read.csv("../../datasets/WBCD/breast_cancer_wisconsin.csv", sep = ",", na.strings = c("NA", "", NA))
data <- data[, -1]
avatar <- read.csv("../../datasets/WBCD/breast_cancer_wisconsin_avatarized_k20.csv", sep = ",", na.strings = c("NA", "", NA))

avatar$Class <- as_factor(avatar$Class)
data$Class <- as_factor(data$Class)

colors <- read.csv("../../color.csv", stringsAsFactors = FALSE)
rownames(colors) <- colors$type

axis_title_size <- 19
axis_text_size <- 19
legend_text_size <- 19
legend_title_size <- 19
axis_text_size_small <- 16

# Merge of datasets 
data_tot <- rbind(data, avatar)

# Perform FAMD with avatar as supplemental indivudals
famd <- FAMD(data_tot, ncp = 5, graph = FALSE, ind.sup = (nrow(data_tot)/2 +1):nrow(data_tot))
res_ind <- as.data.frame(famd$ind$coord)
res_ind_sup <- as.data.frame(famd$ind.sup$coord)
res_ind['type'] <- 'Original'
res_ind_sup['type'] <- 'Avatar'
res_ind_tot <- rbind(res_ind, res_ind_sup)

set.seed(43)
rows <- sample(nrow(res_ind_tot))
res_ind_tot <- res_ind_tot[rows,]

plotAb <- ggplot(res_ind_tot, aes(x = Dim.1, y = Dim.2, fill = type)) +
    # add points
    geom_point(size = 3, shape = 21, alpha = 1) +
    # fill according data source
    aes(fill = factor(type)) +
    scale_fill_manual(values = c(colors["avatar","color"], colors["original","color"])) +
    # add axis label with exlplained variability
    xlab(paste0("Dim. 1 (", round(famd$eig[1,2], 2) ,"%)")) +
    ylab(paste0("Dim. 2 (", round(famd$eig[2,2], 2) ,"%)")) +
    ylim(c(-2.6,5.3)) + 
    # theme and figure details
    theme_bw() +
    theme(legend.position = c(0.9, 0.12), 
          legend.title = element_blank(),
          legend.key.size = unit(0.8, "cm"),
          legend.text = element_text(size = legend_text_size, color = "black", family = "sans"),
          axis.text = element_text(size =  axis_text_size_small, color = "black", family = "sans"),
          axis.title = element_text(size = axis_title_size, color = "black", family = "sans"))


df_scores_70 <- read.csv("../../datasets/results_df/WBCD_f-score_variousSplit.csv")
df_res_70 <- read.csv("../../datasets/results_df/WBCD_perf_variousSplit.csv")

auc_original <- df_res_70[(df_res_70$perf == "auc") & (df_res_70$type == "Original"), "X0"]
auc_avatar <- df_res_70[(df_res_70$perf == "auc") & (df_res_70$type == "Avatar"), "X0"]
    
df_scores_70$type <- as.character(df_scores_70$type)
df_scores_70[df_scores_70$type == "Original", "type"] <- paste0("Original (AUC = ", round(mean(auc_original)*100, 2), ")")
df_scores_70[df_scores_70$type == "Avatar", "type"] <- paste0("Avatar (AUC = ", round(mean(auc_avatar)*100, 2), ")")
df_scores_70$order <- c(5,7,8,3,2,9,6,4,1)

df_scores_70 <- df_scores_70 %>% 
    group_by(type, feature) %>%
    mutate(mean = mean(F.score), sd = sd(F.score))

plotAd <- ggplot(data = df_scores_70, aes(x = reorder(feature, order), y = mean, fill = type)) +
    geom_bar(stat = "identity", position = position_dodge(), color = "black") +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0, lwd = 1, position=position_dodge(.9)) +
    coord_flip() + 
    theme_bw() +
    ylab("F-score") +
    xlab(NULL) +
    scale_fill_manual(name = NULL, values = c(colors["original","color"], colors["avatar","color"]),breaks= c(paste0("Original (AUC = ", round(mean(auc_original)*100, 2), ")"), paste0("Avatar (AUC = ", round(mean(auc_avatar)*100, 2), ")")))  +
    theme(legend.position = c(0.75, 0.1), 
          legend.title = element_blank(),
          legend.key.size = unit(0.8, "cm"),
          legend.text = element_text(size = legend_text_size, color = "black", family = "sans"),
          axis.text.x = element_text(size =  axis_text_size_small, color = "black", family = "sans"),
          axis.text.y = element_text(size =  axis_text_size, color = "black", family = "sans"),
          axis.title = element_text(size = axis_title_size, color = "black", family = "sans"),
          axis.ticks.y = element_blank())

metrics_wbcd <- py$Security_metrics()
metrics_wbcd$fit(data, avatar, nf = 2L)

df_local_cloaking_wbcd <- data.frame(metrics_wbcd$local_cloaking["hit_counts"])
color_vector <- c(rep("#c14f4f", 2), rep(colors['avatar', 'color'], 50))

intersect_median_density_wbcd <- 10.1

plotBb <- ggplot(df_local_cloaking_wbcd, aes(hit_counts.1)) +
    # add histogram and density
    geom_histogram(bins = 52, fill = color_vector, 
                   color = "#1D1D1B") +
    geom_density(aes(y = ..count..), lwd = 1, color = "#1D1D1B") +
    # theme and axis
    xlim(c(-1,100)) +
    theme_minimal() +
    xlab("Local cloaking") +
    ylab("Number of individuals") +
    theme(axis.title = element_text(size = axis_title_size, family = "sans", color = "black"),
          axis.text = element_text(size = axis_text_size_small, family = "sans", color = "black"),
          axis.line = element_line(colour = "black", 
                                   size = 0.8, 
                                   linetype = "solid", 
                                   arrow = arrow(type='closed', length = unit(5,'pt')))) +

    ## add arrow and text 
    # median segment
    geom_segment(aes(x = median(hit_counts.1),
                     y = 0, 
                     xend = median(hit_counts.1),
                     yend = intersect_median_density_wbcd), 
                 color = '#1D1D1B', size = 2) +

    # median arrow and text
    geom_segment(aes(x = median(hit_counts.1) + 25,
                     y = intersect_median_density_wbcd + 15.5, 
                     xend = median(hit_counts.1) + 0.5,
                     yend = intersect_median_density_wbcd), 
                 color = '#1D1D1B', size = 0.8, 
                 arrow = arrow(type='closed', length = unit(10, 'pt'))) +
    geom_label(aes(x =  median(hit_counts.1) + 37.5, 
                   y =  intersect_median_density_wbcd + 16.5,
                   label =paste0("Median = ", median(hit_counts.1))), 
               size = 8, 
               family = "sans", 
               label.size = NA) + 

    # hidden rate arrow and text 
    geom_segment(aes(x = median(hit_counts.1) + 25, 
                     y = intersect_median_density_wbcd + 23, 
                     xend = 0.5, 
                     yend = sum(df_local_cloaking_wbcd[,1] == 0) ), 
                 color = '#1D1D1B', size = 0.8, 
                 arrow = arrow(type='closed', length = unit(10,'pt'))) +
    geom_label(aes(x =  median(hit_counts.1) + 43, 
                   y = intersect_median_density_wbcd + 24,
                   label = paste0("Hidden rate = ", round(metrics_wbcd$hidden_rate, 0), "%")), 
               size = 8, 
               family = "sans", 
               label.size = NA) +
    
    scale_color_manual(name = "statistics", values = c(median = "red"))


avatar_tot <- read.csv('../../datasets/WBCD/wbcd_avatarized_25time_k20_nf_2.csv')
avatar_tot["Class"] <- lapply(avatar_tot["Class"], factor)
pr <- py$Projection()
coord_original <- pr$fit_transform(data, nf = 2L)[[1]]

df_local_cloaking_k2 <- data.frame(row.names = 1:nrow(data))
for (k in unique(avatar_tot$iter)) {
    avatar <- avatar_tot[avatar_tot$iter == k,]
    coord_avatar <- pr$transform(avatar[, 1:(ncol(avatar) - 1)])
    distances <- py$record_to_avatar_distance(coord_original, coord_avatar)
    local_cloaking <- py$local_cloaking(coord_original, coord_avatar, distances)
    df_local_cloaking_k2[paste0("local_cloaking_", k)] <- local_cloaking$hit_counts[, 1]
}

df_local_cloaking_k2_LC0_wbcd <- as.data.frame(table(rowSums(df_local_cloaking_k2 == 0)))
df_local_cloaking_k2_LC0_wbcd["Freq2"] <- paste0(formatC(df_local_cloaking_k2_LC0_wbcd$Freq/sum(df_local_cloaking_k2_LC0_wbcd$Freq)*100, format = "f", digits = 1), "%")

plotBd <- ggplot(df_local_cloaking_k2_LC0_wbcd, aes(x = Var1, y = Freq)) +
    # histogram
    geom_bar(fill = colors["original", "color"], colour = "black", show.legend = FALSE, stat = "identity") + 
    # adding percentage as text on each bar
    geom_text(aes(label = sprintf('%.1f%%', after_stat(df_local_cloaking_k2_LC0_wbcd$Freq / sum(df_local_cloaking_k2_LC0_wbcd$Freq) * 100))),
              vjust = -0.5, size = 6) +
    # axis and theme 
    ylab("Number of individuals") + 
    xlab("Number of local cloaking at zero over 25 avatarizations") +
    theme_minimal() +
    theme(axis.title = element_text(size = axis_title_size, family = "sans", color = "black"),
          axis.text = element_text(size = axis_text_size_small, family = "sans", color = "black"),
          axis.line = element_line(colour = "black", 
                                   size = 0.5, 
                                   linetype = "solid", 
                                   arrow = arrow(type = 'closed', length = unit(5,'pt')))) +
    scale_x_discrete(breaks = c("0", "1", "2", "3", "4", "\U2265 5"), limits = c("0", "1", "2", "3", "4", "\U2265 5"))