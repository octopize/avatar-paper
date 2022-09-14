library(librarian)

# Main libraries
librarian::shelf(mice, reshape2)

# Library to use Python with R
librarian::shelf(
  reticulate,
  # Multidimensional reduction method library
  FactoMineR,
  # Survival analysis
  survival,
  survminer,
  # Display a corrplot
  corrplot
)




# import specific function to compute metrics
py_run_string("import sys")
py_run_string("sys.path.append('..')")
py_run_string("from metrics.privacy_metrics import SecurityMetrics ")
py_run_string("from metrics.privacy_metrics.record_to_avatar_distance import record_to_avatar_distance")
py_run_string("from metrics.privacy_metrics.local_cloaking import get_local_cloaking")

data <- read.csv("../datasets/WBCD/breast_cancer_wisconsin.csv", sep = ",", na.strings = c("NA", "", NA))
data <- data[, -1]
avatar <- read.csv("../datasets/WBCD/breast_cancer_wisconsin_avatarized_k20.csv", sep = ",", na.strings = c("NA", "", NA))
synthpop <- read.csv("../datasets/WBCD/wbcd_synthpop_base.csv", sep = ",", na.strings = c("NA", "", NA))
ctgan <- read.csv("../datasets/WBCD/wbcd_CTGAN_base_2.csv", sep = ",", na.strings = c("NA", "", NA))

avatar$Class <- as.factor(avatar$Class)
data$Class <- as.factor(data$Class)
synthpop$Class <- as.factor(synthpop$Class)
ctgan$Class <- as.factor(ctgan$Class)

colors <- read.csv("../color.csv", stringsAsFactors = FALSE)
rownames(colors) <- colors$type

axis_title_size <- 20
axis_text_size <- 20
legend_text_size <- 20
legend_title_size <- 20

##### ORIGINAL - AVATAR comparison

# Merge of datasets
data_tot <- rbind(data, avatar)

#  Perform FAMD with avatar as supplemental individuals
famd <- FAMD(data_tot, ncp = 5, graph = FALSE, ind.sup = (nrow(data_tot) / 2 + 1):nrow(data_tot))
res_ind <- as.data.frame(famd$ind$coord)
res_ind_sup <- as.data.frame(famd$ind.sup$coord)
res_ind["type"] <- "Original"
res_ind_sup["type"] <- "Avatar"
res_ind_tot <- rbind(res_ind, res_ind_sup)

set.seed(43)
rows <- sample(nrow(res_ind_tot))
res_ind_tot <- res_ind_tot[rows, ]

options(repr.plot.width = 10, repr.plot.height = 7)
plotAb <- ggplot(res_ind_tot, aes(x = Dim.1, y = Dim.2, fill = type)) +
  # add points
  geom_point(size = 3, shape = 21, alpha = 1) +
  #  fill according data source
  aes(fill = factor(type)) +
  scale_fill_manual(values = c(colors["avatar", "color"], colors["original", "color"])) +
  # add axis label with explained variability
  xlab(paste0("Dim. 1 (", round(famd$eig[1, 2], 2), "%)")) +
  ylab(paste0("Dim. 2 (", round(famd$eig[2, 2], 2), "%)")) +
  ylim(c(-2.6, 5.3)) +
  # theme and figure details
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.12),
    legend.title = element_blank(),
    legend.key.size = unit(0.8, "cm"),
    legend.text = element_text(size = legend_text_size, color = "black", family = ""),
    axis.text = element_text(size = axis_text_size, color = "black", family = ""),
    axis.title = element_text(size = axis_title_size, color = "black", family = "")
  )


# ggsave(file="../figures/WBCD_pca_2D.svg", plot=plot, width=10, height=7, dpi = 320)

##### ORIGINAL - SYNTHPOP comparison



#  Perform FAMD with synthpop as supplemental individuals
get_aids_2D_projection = function(data, synthetic_name = "synthetic"){
  famd <- FAMD(data, ncp = 5, graph = FALSE, ind.sup = (nrow(data) / 2 + 1):nrow(data))
  res_ind <- as.data.frame(famd$ind$coord)
  res_ind_sup <- as.data.frame(famd$ind.sup$coord)
  res_ind["type"] <- "Original"
  res_ind_sup["type"] <- synthetic_name
  res_ind_tot <- rbind(res_ind, res_ind_sup)

  set.seed(43)
  rows <- sample(nrow(res_ind_tot))
  res_ind_tot <- res_ind_tot[rows, ]
  return(list("coord" = res_ind_tot, "model" = famd))
  }

get_plot_projection = function(projection, name, save = FALSE){
  options(repr.plot.width = 10, repr.plot.height = 7)
  projection_plot <- ggplot(res_ind_tot, aes(x = Dim.1, y = Dim.2, fill = type)) +
    # add points
    geom_point(size = 3, shape = 21, alpha = 1) +
    #  fill according data source
    aes(fill = factor(type)) +
    scale_fill_manual(values = c(colors["original", "color"], colors[name, "color"])) +
    # add axis label with explained variability
    xlab(paste0("Dim. 1 (", round(famd$eig[1, 2], 2), "%)")) +
    ylab(paste0("Dim. 2 (", round(famd$eig[2, 2], 2), "%)")) +
    ylim(c(-2.6, 5.3)) +
    # theme and figure details
    theme_bw() +
    theme(
      legend.position = c(0.9, 0.12),
      legend.title = element_blank(),
      legend.key.size = unit(0.8, "cm"),
      legend.text = element_text(size = legend_text_size, color = "black", family = ""),
      axis.text = element_text(size = axis_text_size, color = "black", family = ""),
      axis.title = element_text(size = axis_title_size, color = "black", family = "")
    )

    if (save){
      ggsave(file=paste0("../figures/WBCD_",name,"_pca_2D.svg"), plot=projection_plot, width=10, height=7, dpi = 320)
    }

    return(projection_plot)
}

# Merge of datasets
data_tot <- rbind(data, synthpop)
projection <- get_aids_2D_projection(data_tot, synthetic_name="Synthpop" )
projection_plot <- get_plot_projection(projections, name="synthpop")

##### ORIGINAL - CTGAN comparison

# Merge of datasets
data_tot <- rbind(data, ctgan)

#  Perform FAMD with ctgan as supplemental individuals
famd <- FAMD(data_tot, ncp = 5, graph = FALSE, ind.sup = (nrow(data_tot) / 2 + 1):nrow(data_tot))
res_ind <- as.data.frame(famd$ind$coord)
res_ind_sup <- as.data.frame(famd$ind.sup$coord)
res_ind["type"] <- "Original"
res_ind_sup["type"] <- "CTGAN"
res_ind_tot <- rbind(res_ind, res_ind_sup)

set.seed(43)
rows <- sample(nrow(res_ind_tot))
res_ind_tot <- res_ind_tot[rows, ]

options(repr.plot.width = 10, repr.plot.height = 7)
ctgan_plotAb <- ggplot(res_ind_tot, aes(x = Dim.1, y = Dim.2, fill = type)) +
  # add points
  geom_point(size = 3, shape = 21, alpha = 1) +
  #  fill according data source
  aes(fill = factor(type)) +
  scale_fill_manual(values = c(colors["ctgan", "color"], colors["original", "color"])) +
  # add axis label with explained variability
  xlab(paste0("Dim. 1 (", round(famd$eig[1, 2], 2), "%)")) +
  ylab(paste0("Dim. 2 (", round(famd$eig[2, 2], 2), "%)")) +
  ylim(c(-2.6, 5.3)) +
  # theme and figure details
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.12),
    legend.title = element_blank(),
    legend.key.size = unit(0.8, "cm"),
    legend.text = element_text(size = legend_text_size, color = "black", family = ""),
    axis.text = element_text(size = axis_text_size, color = "black", family = ""),
    axis.title = element_text(size = axis_title_size, color = "black", family = "")
  )


# ggsave(file="../figures/WBCD_ctgan_pca_2D.svg", plot=ctgan_plotAb, width=10, height=7, dpi = 320)


##### ORIGINAL - AVATAR comparison

## AUC and F-score results
df_scores_70 <- read.csv("../datasets/results_df/WBCD_f-score_variousSplit.csv")
df_res_70 <- read.csv("../datasets/results_df/WBCD_perf_variousSplit.csv")

auc_original <- df_res_70[(df_res_70$perf == "auc") & (df_res_70$type == "Original"), "X0"]
auc_avatar <- df_res_70[(df_res_70$perf == "auc") & (df_res_70$type == "Avatar"), "X0"]

df_scores_70$type <- as.character(df_scores_70$type)
df_scores_70[df_scores_70$type == "Original", "type"] <- paste0("Original (AUC = ", round(mean(auc_original) * 100, 2), ")")
df_scores_70[df_scores_70$type == "Avatar", "type"] <- paste0("Avatar (AUC = ", round(mean(auc_avatar) * 100, 2), ")")
df_scores_70$order <- c(5, 7, 8, 3, 2, 9, 6, 4, 1)

df_scores_70 <- df_scores_70 %>%
  group_by(type, feature) %>%
  mutate(mean = mean(F.score), sd = sd(F.score))

plotAd <- ggplot(data = df_scores_70, aes(x = reorder(feature, order), y = mean, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0, lwd = 1, position = position_dodge(.9)) +
  coord_flip() +
  theme_bw() +
  ylab("F-score") +
  xlab(NULL) +
  scale_fill_manual(name = NULL, values = c(colors["original", "color"], colors["avatar", "color"]), breaks = c(paste0("Original (AUC = ", round(mean(auc_original) * 100, 2), ")"), paste0("Avatar (AUC = ", round(mean(auc_avatar) * 100, 2), ")"))) +
  theme(
    legend.position = c(0.72, 0.1),
    legend.title = element_blank(),
    legend.key.size = unit(0.8, "cm"),
    legend.text = element_text(size = legend_text_size, color = "black", family = ""),
    axis.text.x = element_text(size = axis_text_size, color = "black", family = ""),
    axis.text.y = element_text(size = axis_text_size, color = "black", family = ""),
    axis.title = element_text(size = axis_title_size, color = "black", family = ""),
    axis.ticks.y = element_blank()
  )

metrics <- py$SecurityMetrics()
metrics$fit(data, avatar, nf = 2L)

avatar_local_cloaking_wbcd <- c()
hit_counts_wbcd <- r_to_py(metrics$local_cloaking$hit_counts)

for (i in seq_along(hit_counts_wbcd)) {
  avatar_local_cloaking_wbcd <- c(avatar_local_cloaking_wbcd, py_to_r(hit_counts_wbcd[i - 1])[1])
}

df_local_cloaking <- data.frame(avatar = avatar_local_cloaking_wbcd)
color_vector <- c(rep("#c14f4f", 2), rep(colors["avatar", "color"], 50))

intersect_median_density <- 10.1
options(repr.plot.width = 10, repr.plot.height = 7)


axis_title_size <- 21
axis_text_size <- 22
legend_text_size <- 19
legend_title_size <- 19

plotBb <- ggplot(df_local_cloaking, aes(avatar)) +
  #  add histogram and density
  geom_histogram(
    bins = 52, fill = color_vector,
    color = "#1D1D1B"
  ) +
  geom_density(aes(y = ..count..), lwd = 1, color = "#1D1D1B") +
  # theme and axis
  xlim(c(-1, 100)) +
  theme_minimal() +
  xlab("Local cloaking") +
  ylab("Number of individuals") +
  theme(
    axis.title = element_text(size = axis_title_size, family = ""),
    axis.text = element_text(size = axis_text_size, family = ""),
    axis.line = element_line(
      colour = "black",
      size = 0.8,
      linetype = "solid",
      arrow = arrow(type = "closed", length = unit(5, "pt"))
    )
  ) +

  ## add arrow and text
  # median segment
  geom_segment(aes(
    x = median(avatar),
    y = 0,
    xend = median(avatar),
    yend = intersect_median_density
  ),
  color = "#1D1D1B", size = 2
  ) +

  # median arrow and text
  geom_segment(aes(
    x = median(avatar) + 25,
    y = intersect_median_density + 15.5,
    xend = median(avatar) + 0.5,
    yend = intersect_median_density
  ),
  color = "#1D1D1B", size = 0.8,
  arrow = arrow(type = "closed", length = unit(10, "pt"))
  ) +
  geom_label(aes(
    x = median(avatar) + 37.5,
    y = intersect_median_density + 15,
    label = paste0("Median = ", median(avatar))
  ),
  size = 9,
  family = "",
  label.size = NA
  ) +

  # hidden rate arrow and text
  geom_segment(aes(
    x = median(avatar) + 25,
    y = intersect_median_density + 30,
    xend = 0.5,
    yend = sum(df_local_cloaking[, 1] == 0)
  ),
  color = "#1D1D1B", size = 0.8,
  arrow = arrow(type = "closed", length = unit(10, "pt"))
  ) +
  geom_label(aes(
    x = median(avatar) + 45,
    y = intersect_median_density + 30,
    label = paste0("Hidden rate = ", round(metrics$hidden_rate, 0), " %")
  ),
  size = 9,
  family = "",
  label.size = NA
  ) +
  scale_color_manual(name = "statistics", values = c(median = "red"))




# ggsave(file="../figures/WBCD_local_cloaking.svg", plot=plot, width=10, height=7, dpi = 320)

# Number of individuals with a local cloaking under 5
prop_wbcd <- sum(df_local_cloaking$avatar <= 5) / dim(df_local_cloaking)[1] * 100

avatar_tot <- read.csv("../datasets/WBCD/wbcd_avatarized_25time_k20_nf_2.csv")
avatar_tot["Class"] <- lapply(avatar_tot["Class"], factor)

#  Projection
coord_original <- FAMD(data, ncp = 2, graph = FALSE)
coord_original <- as.data.frame(coord_original$ind$coord)
coord_original <- reticulate::r_to_py(as.data.frame(coord_original))


df_local_cloaking_k2 <- data.frame(row.names = 1:nrow(data))
for (k in unique(avatar_tot$iter)) {
  avatar <- avatar_tot[avatar_tot$iter == k, ]
  # metrics
  metrics <- py$SecurityMetrics()
  metrics$fit(data, avatar[, -ncol(avatar)], nf = 5L)

  avatar_local_cloaking_tmp <- c()
  hit_counts_wbcd <- r_to_py(metrics$local_cloaking$hit_counts)

  for (i in seq_along(hit_counts_wbcd)) {
    avatar_local_cloaking_tmp <- c(avatar_local_cloaking_tmp, py_to_r(hit_counts_wbcd[i - 1])[1])
  }
  df_local_cloaking_k2[paste0("local_cloaking_", k, "_", 20)] <- data.frame(avatar_local_cloaking_tmp)
}

df_local_cloaking_k2_LC0_total <- as.data.frame(table(rowSums(df_local_cloaking_k2 == 0)))
df_local_cloaking_k2_LC0_total["Freq2"] <- paste0(
  formatC(df_local_cloaking_k2_LC0_total$Freq / sum(
    df_local_cloaking_k2_LC0_total$Freq
  ) * 100, format = "f", digits = 1), "%"
)

cum_freq <- df_local_cloaking_k2_LC0_total[, ]
levels(cum_freq$Var1) <- c(levels(cum_freq$Var1), "\U2265 10")
cum_freq$Var1[cum_freq$Var1 == 10] <- "\U2265 10"


# number of zero local cloaking over 25 avatarizations

options(repr.plot.width = 10, repr.plot.height = 7)
plotBd <- ggplot(cum_freq, aes(x = Var1, y = Freq)) +
  #  histogram
  geom_bar(fill = colors["original", "color"], colour = "black", show.legend = FALSE, stat = "identity") +
  #  adding percentage as text on each bar
  geom_text(aes(label = sprintf(
    "%.1f%%", after_stat(cum_freq$Freq / sum(cum_freq$Freq) * 100)
  )),
  vjust = -0.5, size = 6
  ) +
  # axis and theme
  ylab("Number of individuals") +
  xlab("Number of local cloaking at zero over 25 avatarizations") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = axis_title_size, family = ""),
    axis.text = element_text(size = axis_text_size, family = ""),
    axis.line = element_line(
      colour = "black",
      size = 0.5,
      linetype = "solid",
      arrow = arrow(type = "closed", length = unit(5, "pt"))
    )
  ) +
  scale_x_discrete(
    breaks = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "\U2265 10"),
    limits = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "\U2265 10")
  )



# ggsave(file="../figures/WBCD_zeroCloaking_percentage.svg", plot=plot, width=10, height=7, dpi = 320)


print("execution is done")
