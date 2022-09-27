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
  corrplot,
  stringr,
  tibble,
  dplyr
)




# import specific function to compute metrics
py_run_string("import sys")
py_run_string("sys.path.append('..')")
py_run_string("from metrics.privacy_metrics import SecurityMetrics ")
py_run_string("from metrics.privacy_metrics.record_to_avatar_distance import record_to_avatar_distance")
py_run_string("from metrics.privacy_metrics.local_cloaking import get_local_cloaking")


axis_title_size <- 20
axis_text_size <- 20
legend_text_size <- 20
legend_title_size <- 20
save <- FALSE


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



#  Perform FAMD with synthetic as supplemental individuals

get_aids_2D_projection <- function(data, synthetic_name = "synthetic") {
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

get_plot_projection <- function(projection, name, save = FALSE) {
  options(repr.plot.width = 10, repr.plot.height = 7)

  if (name == "ctgan") {
    category <- "CT-GAN"
  }
  if (name == "avatar") {
    category <- "Avatar"
  }
  if (name == "synthpop") {
    category <- "Synthpop"
  }
  cols <- c(colors["original", "color"], colors[name, "color"])
  names(cols) <- c("Original", category)

  
  projection_plot <- ggplot(projection$coord, aes(x = Dim.1, y = Dim.2, fill = type)) +
    # add points
    geom_point(size = 3, shape = 21, alpha = 1) +
    #  fill according data source
    aes(fill = factor(type)) +
    scale_fill_manual(values = cols) +
    # add axis label with explained variability
    xlab(paste0("Dim. 1 (", round(projection$model$eig[1, 2], 2), "%)")) +
    ylab(paste0("Dim. 2 (", round(projection$model$eig[2, 2], 2), "%)")) +
    ylim(c(-2.6, 5.3)) +
    # theme and figure details
    theme_bw() +
    theme(
      legend.position = c(0.9, 0.12),
      legend.title = element_blank(),
      legend.key.size = unit(0.8, "cm"),
      legend.text = element_text(size = legend_text_size, color = "black", family = "sans"),
      axis.text = element_text(size = axis_text_size, color = "black", family = "sans"),
      axis.title = element_text(size = axis_title_size, color = "black", family = "sans")
    )

  if (save) {
    ggsave(file = paste0("../figures/WBCD_", name, "_pca_2D.svg"), plot = projection_plot, width = 10, height = 7, dpi = 320)
  }

  return(projection_plot)
}
##### ORIGINAL - AVATAR comparison
data_tot <- rbind(data, avatar)
avatar_projection <- get_aids_2D_projection(data_tot, synthetic_name = "Avatar")
plotAb <- get_plot_projection(avatar_projection, name = "avatar")


##### ORIGINAL - SYNTHPOP comparison
data_tot <- rbind(data, synthpop)
synthpop_projection <- get_aids_2D_projection(data_tot, synthetic_name = "Synthpop")
projection_plot_synthpop <- get_plot_projection(synthpop_projection, name = "synthpop")

##### ORIGINAL - CTGAN comparison
data_tot <- rbind(data, ctgan)
ctgan_projection <- get_aids_2D_projection(data_tot, synthetic_name = "CT-GAN")
projection_plot_ctgan <- get_plot_projection(ctgan_projection, name = "ctgan")





## AUC and F-score analysis
reshape_predictions_and_f1score <- function(predictions, f1_score, synthetic_name) {
  auc_original <- predictions[(predictions$perf == "auc") & (predictions$type == "Original"), "X0"]
  auc_synthetic <- predictions[(predictions$perf == "auc") & (predictions$type == synthetic_name), "X0"]

  f1_score$type <- as.character(f1_score$type)
  f1_score[f1_score$type == "Original", "type"] <- paste0("Original (AUC = ", round(mean(auc_original) * 100, 2), ")")
  f1_score[f1_score$type == synthetic_name, "type"] <- paste0(synthetic_name, " (AUC = ", round(mean(auc_synthetic) * 100, 2), ")")
  f1_score$order <- c(5, 7, 8, 3, 2, 9, 6, 4, 1) # original descending order

  return(list("f1_score" = f1_score, "auc_original" = auc_original, "auc_synthetic" = auc_synthetic))
}

get_f1score_plot <- function(scores_and_auc, synthetic_color) {
  stats_f1_scores <- scores_and_auc$f1_score %>%
    group_by(type, feature) %>%
    mutate(mean = mean(F.score), sd = sd(F.score))

  if (synthetic_color == "ctgan") {
    category <- paste0("CT-GAN (AUC = ", round(mean(scores_and_auc$auc_synthetic) * 100, 2), ")")
  }
  if (synthetic_color == "avatar") {
    category <- paste0("Avatar (AUC = ", round(mean(scores_and_auc$auc_synthetic) * 100, 2), ")")
  }
  if (synthetic_color == "synthpop") {
    category <- paste0("Synthpop (AUC = ", round(mean(scores_and_auc$auc_synthetic) * 100, 2), ")")
  }
  cols <- c(colors["original", "color"], colors[synthetic_color, "color"])
  original_category <- paste0("Original (AUC = ", round(mean(scores_and_auc$auc_original) * 100, 2), ")")
  names(cols) <- c(original_category, category)

  plot <- ggplot(data = stats_f1_scores, aes(x = reorder(feature, order), y = mean, fill = type)) +
    geom_bar(stat = "identity", position = position_dodge(), color = "black") +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0, lwd = 1, position = position_dodge(.9)) +
    coord_flip() +
    theme_bw() +
    ylab("F-score") +
    xlab(NULL) +
    scale_fill_manual(
      name = NULL,
      values = cols,
      # breaks = c(
      #   paste0("Original (AUC = ", round(mean(scores_and_auc$auc_original) * 100, 2), ")"),
      #   paste0("Avatar (AUC = ", round(mean(scores_and_auc$auc_synthetic) * 100, 2), ")"))
    ) +
    theme(
      legend.position = c(0.72, 0.1),
      legend.title = element_blank(),
      legend.key.size = unit(0.8, "cm"),
      legend.text = element_text(size = legend_text_size, color = "black", family = "sans"),
      axis.text.x = element_text(size = axis_text_size, color = "black", family = "sans"),
      axis.text.y = element_text(size = axis_text_size, color = "black", family = "sans"),
      axis.title = element_text(size = axis_title_size, color = "black", family = "sans"),
      axis.ticks.y = element_blank()
    )
  return(plot)
}

f1_scores_comparison <- read.csv("../datasets/results_df/comparative_f1_scores.csv")
predictions_comparison <- read.csv("../datasets/results_df/comparative_predictions.csv")
##### ORIGINAL - AVATAR comparison

f1_scores_avatar <- subset(f1_scores_comparison, type == "Original" | type == "Avatar")
predictions_avatar <- subset(predictions_comparison, type == "Original" | type == "Avatar")
scores_and_auc_avatar <- reshape_predictions_and_f1score(predictions_avatar, f1_scores_avatar, synthetic_name = "Avatar")
plotAd <- get_f1score_plot(scores_and_auc_avatar, synthetic_color = "avatar")


##### ORIGINAL - Synthpop comparison

f1_scores_synthpop <- subset(f1_scores_comparison, type == "Original" | type == "Synthpop")
predictions_synthpop <- subset(predictions_comparison, type == "Original" | type == "Synthpop")
scores_and_auc_synthpop <- reshape_predictions_and_f1score(predictions_synthpop, f1_scores_synthpop, synthetic_name = "Synthpop")
f1score_plot_synthpop <- get_f1score_plot(scores_and_auc_synthpop, synthetic_color = "synthpop")

##### ORIGINAL - ctgan comparison

f1_scores_ctgan <- subset(f1_scores_comparison, type == "Original" | type == "CT-GAN")
predictions_ctgan <- subset(predictions_comparison, type == "Original" | type == "CT-GAN")
scores_and_auc_ctgan <- reshape_predictions_and_f1score(predictions_ctgan, f1_scores_ctgan, synthetic_name = "CT-GAN")
f1score_plot_ctgan <- get_f1score_plot(scores_and_auc_ctgan, synthetic_color = "ctgan")

## Comparative utility plot 
vec <- c("Bare Nuclei", "Uniformity of Cell Shape", "Uniformity of Cell Size",
         "Bland Chromatin", "Clump Thickness", "Normal Nucleoli", "Marginal Adhesion",
         "Single Epithelial Cell Size", "Mitoses") 


f1_original <- scores_and_auc_avatar$f1_score %>% 
  filter(str_detect(type, "^Original"))        

# Reorder data frame
f1_original <- left_join(data.frame(feature = vec),f1_original,by = "feature")

f1_avatar <- scores_and_auc_avatar$f1_score %>%
  filter(str_detect(type, "^Avatar")) 

# Reorder data frame
f1_avatar <- left_join(data.frame(feature = vec),f1_avatar,by = "feature")

f1_synthpop <- scores_and_auc_synthpop$f1_score %>%
  filter(str_detect(type, "^Synthpop"))        
# Reorder data frame
f1_synthpop <- left_join(data.frame(feature = vec),f1_synthpop,by = "feature")

f1_ctgan <- scores_and_auc_ctgan$f1_score %>%
  filter(str_detect(type, "^CT-GAN"))        
# Reorder data frame
f1_ctgan <- left_join(data.frame(feature = vec),f1_ctgan,by = "feature")

levels <- c(f1_original$type[1], f1_avatar$type[1], f1_synthpop$type[1], f1_ctgan$type[1] )
levels <-gsub('\\(AUC', '\nAUC', levels )
levels <-gsub(')', '', levels )

F1_scores <- rbind(f1_original, f1_avatar, f1_synthpop, f1_ctgan)
F1_scores$feature <- factor(F1_scores$feature, levels=vec )
F1_scores$type <-gsub('\\(AUC', '\nAUC', F1_scores$type )
F1_scores$type <-gsub(')', '', F1_scores$type )
F1_scores$type = factor(F1_scores$type, levels=rev(levels)) # order levels

wrapping_size <- 20
cols = c(colors['original', 'color'],
         colors['avatar', 'color'],
         colors['synthpop', 'color'],
         colors['ctgan', 'color'])
color_names <-gsub('\\(AUC', '\nAUC', c(levels[1],levels[2],levels[3],levels[4] ))
color_names <- gsub(')', '', color_names)
names(cols) <- color_names


wbcd_comparative_utility <- ggplot(F1_scores, aes(x = F.score, y = type , fill = type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = cols )+
  theme_bw() +
  ylab(NULL) +
  xlab("F scores") +
  facet_grid(feature~.,switch = 'y') +
  theme(text = element_text(size = 16), 
        strip.text.y.left = element_text(angle = 0),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.text = element_text(size = legend_text_size -2, color = "black", family = "sans"),
        axis.text = element_text(size = axis_text_size, color = "black", family = "sans"),
        axis.title = element_text(size = axis_title_size, color = "black", family = "sans"),
        legend.background = element_rect(fill = "white", linetype = "solid"),
        legend.title=element_blank(),
        strip.text.y =  element_text(size = 11, color = "black", family = "sans"),
        legend.spacing.y = unit(0.8, 'cm'),
       ) +
  guides(fill = guide_legend(byrow = TRUE))

if (save) {
  ggsave(file = "../figures/wbcd_comparative_plot.svg", plot = wbcd_comparative_utility, width = 10, height = 7, dpi = 290)
}

dcr_nndr_results_wbcd <- read.csv('../datasets/results_df/WBCD_DCR_NNDR_comparison_results.csv')
wbcd_comparative_privacy <- get_table_plot(dcr_nndr_results_wbcd, save = TRUE, title = 'WBCD DCR-NNDR Results')


## Specific Privacy metrics computation
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
    axis.title = element_text(size = axis_title_size, family = "sans"),
    axis.text = element_text(size = axis_text_size, family = "sans"),
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
  family = "sans",
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
  family = "sans",
  label.size = NA
  ) +
  scale_color_manual(name = "statistics", values = c(median = "red"))



if (save){
  ggsave(file="../figures/WBCD_local_cloaking.svg", plot=plot, width=10, height=7, dpi = 320)
}

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
    axis.title = element_text(size = axis_title_size, family = "sans"),
    axis.text = element_text(size = axis_text_size, family = "sans"),
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


if (save) {
 ggsave(file="../figures/WBCD_zeroCloaking_percentage.svg", plot=plot, width=10, height=7, dpi = 320)
}


print("execution is done")
