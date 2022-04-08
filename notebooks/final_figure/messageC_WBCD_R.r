
library(librarian)

# Main libraries
librarian::shelf(mice, reshape2) # tidyverse,

# Library to use Python with R
librarian::shelf(
  reticulate,
  # Multidimensional reduction method libary
  FactoMineR,
  # visualization
  ggplot2,
  # data table manipulation
  tibble,
  dplyr,
)


# import specific function to compute metrics
py_run_string("import sys")
py_run_string("sys.path.append('../..')")

py_run_string("from lsg.security_metrics import SecurityMetrics")
py_run_string("from lsg.security_metrics.record_to_avatar_distance import record_to_avatar_distance")
py_run_string("from lsg.security_metrics.local_cloaking import get_local_cloaking")
py_run_string("from lsg.security_metrics.avatars_are_k_hit import avatars_are_k_hit")

save <- FALSE

file_data <- "../../datasets/"
file_results <- "../../figure/"

data <- read.csv(paste0(file_data, "/WBCD/breast_cancer_wisconsin.csv"), sep = ",", na.strings = c("NA", "", NA))
data <- data[, -1]

seq_k <- c(4, 50, 100, 150, 200, 300)

colors <- read.csv("../../color.csv", stringsAsFactors = FALSE)
rownames(colors) <- colors$type

axis_title_size <- 16
axis_text_size <- 15
legend_text_size <- 16
legend_title_size <- 16

# Small k
k <- 4
iter <- 5

tryCatch(
  {
    avatar_tot <- read.table(paste0(file_data, "tmp/wbcd_avatarized/k", k, "_wbcd_avatarized.csv"),
      sep = ",", header = TRUE
    )
  },
  error = function(cnd) {
    print("The dataset was not extracted")
    untar(paste0(file_data, "WBCD/wbcd_avatarized.tar.bz2"),
      files = paste0("./wbcd_avatarized/k", k, "_wbcd_avatarized.csv"),
      exdir = paste0(file_data, "tmp/")
    )
    avatar_tot <- read.table(paste0(file_data, "tmp/wbcd_avatarized/k", k, "_wbcd_avatarized.csv"),
      sep = ",", header = TRUE
    )
  }
)

avatar_small <- avatar_tot[avatar_tot$iter_k == iter, ]

# high k
k <- 342
iter <- 5
tryCatch(
  {
    avatar_tot <- read.table(paste0(file_data, "tmp/wbcd_avatarized/k", k, "_wbcd_avatarized.csv"),
      sep = ",", header = TRUE
    )
  },
  error = function(cnd) {
    print("The dataset was not extracted")
    untar(paste0(file_data, "WBCD/wbcd_avatarized.tar.bz2"),
      files = paste0("./wbcd_avatarized/k", k, "_wbcd_avatarized.csv"),
      exdir = paste0(file_data, "tmp/")
    )
    avatar_tot <- read.table(paste0(file_data, "tmp/wbcd_avatarized/k", k, "_wbcd_avatarized.csv"),
      sep = ",", header = TRUE
    )
  }
)
avatar_high <- avatar_tot[avatar_tot$iter_k == iter, ]

# binding original, small k avatar and high k avatar datasets.
data_small_high <- rbind(data[, 1:10], avatar_small[, 1:10], avatar_high[, 1:10])
data_small_high$Class <- as.factor(data_small_high$Class)

row.names(data_small_high) <- 1:nrow(data_small_high)

# Perform FAMD in the priginal data space
FAMD <- FAMD(data_small_high, ncp = 5, graph = FALSE, ind.sup = (nrow(data_small_high) / 3 + 1):nrow(data_small_high))
res_ind <- as.data.frame(FAMD$ind.sup$coord)[1:683, ]
res_ind_sup <- as.data.frame(FAMD$ind.sup$coord)[684:1366, ]

res_ind_2D <- res_ind[, 1:2]
res_ind_2D["type"] <- "k = 4"
res_ind_sup["type"] <- "k = 342"
res_ind_2D <- rbind(res_ind_2D, res_ind_sup[, c("Dim.1", "Dim.2", "type")])

set.seed(42)
rows <- sample(nrow(res_ind_2D))
res_ind_2D <- res_ind_2D[rows, ]

options(repr.plot.width = 10, repr.plot.height = 7)

plotCb <- ggplot(res_ind_2D, aes(x = Dim.1, y = Dim.2, fill = type)) +
  geom_point(size = 3, shape = 21, alpha = 1) +
  aes(fill = factor(type)) +
  scale_fill_manual(values = c("k = 4" = colors["k_small", "color"], "k = 342" = colors["k_high", "color"])) +
  xlab(paste0("Dim. 1 (", round(FAMD$eig[1, 2], 2), "%)")) +
  ylab(paste0("Dim. 2 (", round(FAMD$eig[2, 2], 2), "%)")) +
  labs(fill = "") +
  theme_bw() +
  ylim(c(-2.6, 5.3)) +
  theme(
    legend.position = c(0.9, 0.12),
    legend.title = element_blank(),
    legend.key.size = unit(0.8, "cm"),
    legend.text = element_text(size = legend_text_size, color = "black", family = ""),
    axis.text = element_text(size = axis_text_size, color = "black", family = ""),
    axis.title = element_text(size = axis_title_size, color = "black", family = "")
  )

# ggsave(file="../../figure/WBCD_pca2D_Ksmall_high.svg", plot = plot, width = 10, height = 7, dpi = 320)

df <- read.csv("../../datasets/results_df/original_data_perf_wbcd.csv")
auc_original <- df[df$test == "auc", "X0"]

result_k <- read.csv("../../datasets/results_df/accuracy_k_wbcd.csv")




# subset to plot only avatar train - avatar test
# (and remove avatar train, original test)
result_k <- result_k[result_k["test_type"] == "avatar", ]
result_k <- result_k[result_k["k"] < 343, ]

result_k <- result_k %>%
  group_by(test_type, k) %>%
  mutate(
    mean = mean(auc), sd = sd(auc), ci_low = mean(auc) - 1.96 * sd(auc) / sqrt(n()), ci_high = mean(auc) + 1.96 * sd(auc) / sqrt(n()),
    percentile_97_5 = quantile(x = auc, probs = 0.975), percentile_2_5 = quantile(x = auc, probs = 0.025)
  )

result_k$color <- colors["avatar", "color"]

df2 <- data.frame(
  x = c(min(result_k$k), max(result_k$k)),
  sd_low = mean(auc_original) - sd(auc_original), sd_high = mean(auc_original) + sd(auc_original),
  ci_low = mean(auc_original) - 1.96 * sd(auc_original) / sqrt(length(auc_original)), ci_high = mean(auc_original) + 1.96 * sd(auc_original) / sqrt(length(auc_original)),
  percentile_97_5 = as.numeric(as.character(quantile(x = auc_original, probs = 0.975))), percentile_2_5 = as.numeric(as.character(quantile(x = auc_original, probs = 0.025)))
)

plotCd <- ggplot() +
  geom_ribbon(data = df2, aes(x = x, ymin = sd_low, ymax = sd_high), fill = colors["original", "color"], alpha = 0.07) +
  geom_ribbon(data = result_k, aes(x = k, ymin = mean - sd, ymax = mean + sd), fill = colors["avatar", "color"], alpha = 0.07) +
  geom_segment(aes(x = min(result_k$k), xend = max(result_k$k), y = mean(auc_original), yend = mean(auc_original), color = colors["original", "color"]), size = 1) +
  geom_line(data = result_k, aes(x = k, y = mean, color = color), size = 1) +
  geom_point(data = result_k, aes(x = k, y = auc), size = 3, color = "white", fill = colors["avatar", "color"], pch = 21) +
  # labs theme and axis
  ylab("Area under the ROC curve") +
  xlab("Number of nearest neighbors (k)") +
  # scale_x_continuous(breaks = seq_k)  +
  scale_y_continuous(breaks = seq(0.9, 1, 0.02), limits = c(0.896, 1)) +
  scale_color_identity(guide = "legend", name = NULL, breaks = c(colors["original", "color"], colors["avatar", "color"]), labels = c("Original", "Avatar")) +
  theme_minimal() +
  theme(
    legend.position = c(0.12, 0.15),
    axis.title.y = element_text(angle = 90, vjust = 2, size = axis_title_size, family = ""),
    axis.text = element_text(size = axis_text_size, colour = "black", family = ""),
    axis.title.x = element_text(size = axis_title_size, family = ""),
    legend.text = element_text(size = legend_text_size, family = ""),
    legend.key.size = unit(1, "cm"),
    axis.line = element_line(
      colour = "black", size = 0.5, linetype = "solid",
      arrow = arrow(type = "closed", length = unit(5, "pt"))
    )
  )


df_local_cloaking_k100 <- data.frame(row.names = 1:nrow(data))
data$Class <- as.factor(data$Class)


#  Original projection
coord_original <- FAMD(data, ncp = 2, graph = FALSE)
coord_original <- as.data.frame(coord_original$ind$coord)

for (k in seq_k) {
  iter <- 1
  tryCatch(
    {
      avatar_tot <- read.table(paste0(file_data, "tmp/wbcd_avatarized/k", k, "_wbcd_avatarized.csv"),
        sep = ",", header = TRUE
      )
    },
    error = function(cnd) {
      print("The dataset was not extracted")
      untar(paste0(file_data, "WBCD/wbcd_avatarized.tar.bz2"),
        files = paste0("./wbcd_avatarized/k", k, "_wbcd_avatarized.csv"),
        exdir = paste0(file_data, "tmp/")
      )
      avatar_tot <- read.table(paste0(file_data, "tmp/wbcd_avatarized/k", k, "_wbcd_avatarized.csv"),
        sep = ",", header = TRUE
      )
    }
  )


  avatar <- avatar_tot[avatar_tot$iter_k == iter, ]
  avatar <- subset(avatar, select = -c(k, id_ind, iter_k))

  metrics_tmp <- py$SecurityMetrics()
  metrics_tmp$fit(data, avatar, nf = 5L)

  hit_counts <- r_to_py(metrics_tmp$local_cloaking$hit_counts)
  avatar_local_cloaking <- c()

  for (i in seq_along(hit_counts)) {
    avatar_local_cloaking <- c(avatar_local_cloaking, py_to_r(hit_counts[i - 1])[1])
  }

  df_local_cloaking_k100[paste0("local_cloaking_", k)] <- avatar_local_cloaking
}

df_local_cloaking_k100["Class"] <- data$Class

df_melt <- melt(df_local_cloaking_k100, id.vars = "Class")


options(repr.plot.width = 10, repr.plot.height = 7)

plotCf <- ggplot(df_melt, aes(x = variable, y = value)) +
  geom_violin(trim = TRUE, alpha = 0.12, fill = colors["avatar", "color"]) +
  geom_boxplot(width = 0.4, fill = colors["avatar", "color"], outlier.size = 0.1, alpha = 0.8) +
  scale_x_discrete(labels = as.character(seq_k)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), breaks = c(1, 5, 10, 25, 50, 100, 300, 600)) +
  xlab("Number of nearest neighbors (k)") +
  ylab("Local Cloaking") +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.title.y = element_text(angle = 90, vjust = 1, size = axis_title_size),
    axis.text = element_text(size = axis_text_size),
    axis.title.x = element_text(size = axis_title_size),
    axis.line = element_line(
      colour = "black", size = 0.5, linetype = "solid",
      arrow = arrow(type = "closed", length = unit(5, "pt"))
    )
  ) +
  coord_cartesian(ylim = c(0, 700))



# ggsave(file = "../../figure/WBCD_localCloaking_violinplot.svg", plot = plot, width = 10, height = 7)

print("execution is done")
