library(librarian)

# Main libraries
librarian::shelf(mice, reshape2) # tidyverse,

# Library to use Python with R
librarian::shelf(
  reticulate,
  # Multidimensional reduction method libary
  FactoMineR,
  # Survival analysis
  survival,
  survminer,
  # data table manipulation
  tibble,
  dplyr,
  # svg visualizationD
  svglite
)




# import specific function to compute metrics
py_run_string("import sys")
py_run_string("sys.path.append('../..')")

py_run_string("from lsg.security_metrics import SecurityMetrics")
py_run_string("from lsg.security_metrics.record_to_avatar_distance import record_to_avatar_distance")
py_run_string("from lsg.security_metrics.local_cloaking import get_local_cloaking")
py_run_string("from lsg.security_metrics.avatars_are_k_hit import avatars_are_k_hit")
py_run_string("from lsg.security_metrics.hidden_rate import hidden_rate")

file_data <- "../../datasets/"
file_results <- "../../figure/"

data <- read.csv(paste0(file_data, "AIDS/aids_original_data.csv"), sep = ";", na.strings = c("NA", "", NA))
data <- data[, -1]

save <- FALSE
seq_k <- c(4, 200, 400, 600, 800, 1000)

colors <- read.csv("../../color.csv", stringsAsFactors = FALSE)
rownames(colors) <- colors$type

axis_title_size <- 16
axis_text_size <- 15
legend_text_size <- 16
legend_title_size <- 16

data_imp <- complete(mice(data, m = 1, maxit = 50, seed = 500, printFlag = FALSE), 1)
data$cd496 <- data_imp$cd496

categorical <- c("hemo", "homo", "drugs", "karnof", "oprior", "z30", "zprior", "race", "gender", "str2", "strat", "symptom", "treat", "offtrt", "r", "cens", "arms")
data[categorical] <- lapply(data[categorical], factor)

write_pvalue <- function(val, dec) {
  if (is.na(val)) {
    return(NA)
  } else if (val < 0.001) {
    return("<0.001")
  } else {
    return(formatC(val, format = "f", digits = dec))
  }
}

# Small k
k <- 4
iter <- 0

tryCatch(
  {
    avatar_tot <- read.table(paste0(file_data, "tmp/aids_avatarized/k", k, "_aids_avatarized.csv"),
      sep = ",", header = TRUE
    )
  },
  error = function(cnd) {
    print("The dataset was not extracted")
    untar(paste0(file_data, "AIDS/aids_avatarized.tar.bz2"),
      files = paste0("./aids_avatarized/k", k, "_aids_avatarized.csv"),
      exdir = paste0(file_data, "tmp/")
    )
    avatar_tot <- read.table(paste0(file_data, "tmp/aids_avatarized/k", k, "_aids_avatarized.csv"),
      sep = ",", header = TRUE
    )
  }
)
avatar_small <- avatar_tot[avatar_tot$iter_k == iter, ]

# Imputation
avatar_small_imp <- complete(mice(avatar_small, m = 1, maxit = 50, seed = 500, printFlag = FALSE), 1)
avatar_small$cd496 <- avatar_small_imp$cd496
avatar_small[categorical] <- lapply(avatar_small[categorical], factor)

# High k = N/2
k <- 1166
iter <- 1
tryCatch(
  {
    avatar_tot <- read.table(paste0(file_data, "tmp/aids_avatarized/k", k, "_aids_avatarized.csv"),
      sep = ",", header = TRUE
    )
  },
  error = function(cnd) {
    print("The dataset was not extracted")
    untar(paste0(file_data, "AIDS/aids_avatarized.tar.bz2"),
      files = paste0("./aids_avatarized/k", k, "_aids_avatarized.csv"),
      exdir = paste0(file_data, "tmp/")
    )
    avatar_tot <- read.table(paste0(file_data, "tmp/aids_avatarized/k", k, "_aids_avatarized.csv"),
      sep = ",", header = TRUE
    )
  }
)
avatar_med <- avatar_tot[avatar_tot$iter_k == iter, ]

# Imputation
avatar_med_imp <- complete(mice(avatar_med, m = 1, maxit = 50, seed = 500, printFlag = FALSE), 1)
avatar_med$cd496 <- avatar_med_imp$cd496
avatar_med[categorical] <- lapply(avatar_med[categorical], factor)

#  binding data, avatar small and avatar high
data_small_med <- rbind(data[, 1:26], avatar_small[, 1:26], avatar_med[, 1:26])
row.names(data_small_med) <- 1:nrow(data_small_med)

# perform projection
famd <- FAMD(data_small_med, ncp = 5, graph = FALSE, ind.sup = (nrow(data_small_med) / 3 + 1):nrow(data_small_med))

res_ind_2D <- as.data.frame(famd$ind.sup$coord)[1:(nrow(data_small_med) / 3), c("Dim.1", "Dim.2")]
res_ind_sup <- as.data.frame(famd$ind.sup$coord)[(nrow(data_small_med) / 3 + 1):nrow(famd$ind.sup$coord), ]
res_ind_2D["type"] <- "k= 4"
res_ind_sup["type"] <- "k= 1166"
res_ind_2D <- rbind(res_ind_2D, res_ind_sup[, c("Dim.1", "Dim.2", "type")])

set.seed(42)
# Sample rows for the visualization
rows <- sample(nrow(res_ind_2D))
res_ind_2D <- res_ind_2D[rows, ]

options(repr.plot.width = 10, repr.plot.height = 7)

plotCa <- ggplot(res_ind_2D, aes(x = Dim.1, y = Dim.2, fill = type)) +
  # point
  geom_point(size = 3, shape = 21, alpha = 1) +
  #  fill by k value
  aes(fill = factor(type)) +
  scale_fill_manual(values = c("k= 4" = colors["k_small", "color"], "k= 1166" = colors["k_high", "color"])) +
  #  theme and details
  xlab(paste0("Dim. 1 (", formatC(famd$eig[1, 2], format = "f", digits = 2), "%)")) +
  ylab(paste0("Dim. 2 (", formatC(famd$eig[2, 2], format = "f", digits = 2), "%)")) +
  labs(fill = "") +
  theme_bw() +
  theme(
    legend.position = c(0.86, 0.15),
    legend.text = element_text(size = legend_text_size, color = "black", family = ""),
    legend.key.size = unit(0.8, "cm"),
    axis.text = element_text(size = axis_text_size, color = "black", family = ""),
    axis.title = element_text(size = axis_title_size, color = "black", family = ""),
    legend.background = element_rect(fill = "white", linetype = "solid")
  )


# ggsave(file="../../figure/aids_famd2D_k_small_high.svg", plot=plot, width=10, height=7, dpi = 320)

data[] <- lapply(data, function(x) {
  if (is.factor(x)) as.numeric(as.character(x)) else x
})

n_arms <- as.vector(table(data$arms))
cens_arms <- paste0(
  as.vector(table(data$cens, data$arms)[2, ]), " (",
  as.vector(round(prop.table(table(data$cens, data$arms), 2)[2, ] * 100, 1)), ")"
)

summary_cox <- summary(coxph(Surv(time = days / 7, event = cens) ~ factor(arms), data = data))

HR_arms <- c("ref", round(as.vector(summary_cox$coefficients[, 2]), 2))
CI_low_arms <- c(NA, round(as.vector(summary_cox$conf.int[, 3]), 2))
CI_high_arms <- c(NA, round(as.vector(summary_cox$conf.int[, 4]), 2))
p_arms <- c(NA, sapply(as.vector(summary_cox$coefficients[, 5]), write_pvalue, dec = 3))

table2_data <- t(data.frame(
  "N" = n_arms, "Event" = cens_arms, "Hazard_ratio" = HR_arms,
  "Confidence_interval_95" = c(NA, paste0(CI_low_arms[-1], "-", CI_high_arms[-1])),
  "p-value" = p_arms, row.names = c("0", "1", "2", "3")
))
table2_data <- tibble(as.data.frame(table2_data))
rownames(table2_data) <- c("N", "Event", "Hazard ratio", "95% confidence interval", "p-value")
table2_data

data[categorical] <- lapply(data[categorical], factor)
coord_original <- FAMD(data, ncp = 5, graph = FALSE)
coord_original <- as.data.frame(coord_original$ind$coord)


compute <- FALSE
# results are saved in a df
if (compute) {
  ## parallele computing (not requiered if you import results datasets)
  library(doParallel)
  library(parallel)
  library(foreach)

  parralele_function <- function(k) {
    tryCatch(
      {
        avatar_tot <- read.table(paste0(file_data, "tmp/aids_avatarized/k", k, "_aids_avatarized.csv"),
          sep = ",", header = TRUE
        )
      },
      error = function(cnd) {
        print("The dataset was not extracted")
        untar(paste0(file_data, "AIDS/aids_avatarized.tar.bz2"),
          files = paste0("./aids_avatarized/k", k, "_aids_avatarized.csv"),
          exdir = paste0(file_data, "tmp/")
        )
        avatar_tot <- read.table(paste0(file_data, "tmp/aids_avatarized/k", k, "_aids_avatarized.csv"),
          sep = ",", header = TRUE
        )
      }
    )

    for (iter in 0:9) {
      avatar <- avatar_tot[avatar_tot$iter_k == iter, 1:26]

      ## imputation
      avatar_imp <- tidyr::complete(mice::mice(avatar, m = 1, maxit = 50, seed = 500, printFlag = FALSE, print = FALSE), 1)
      avatar$cd496 <- avatar_imp$cd496

      ## survival analysis
      summary_cox <- summary(survival::coxph(survival::Surv(time = days / 7, event = cens) ~ factor(arms), data = avatar))
      res_HR <- c(round(as.vector(summary_cox$coefficients[, 2]), 2))
      res_HR_CI_low <- c(round(as.vector(summary_cox$conf.int[, 3]), 2))
      res_HR_CI_high <- c(round(as.vector(summary_cox$conf.int[, 4]), 2))
      res_pval <- sapply(as.vector(summary_cox$coefficients[, 5]), write_pvalue, dec = 3)

      res_row <- c(k, iter, res_HR, res_HR_CI_low, res_HR_CI_high, res_pval)
      res_df[nrow(res_df) + 1, ] <- res_row
    }
    return(res_df)
  }

  # Initialization
  res_df <- data.frame(matrix(ncol = 14, nrow = 0))

  colnames(res_df) <- c(
    "k", "iter", "res_HR1", "res_HR2", "res_HR3",
    "res_HR_CI_low1", "res_HR_CI_low2", "res_HR_CI_low3",
    "res_HR_CI_high1", "res_HR_CI_high2", "res_HR_CI_high3",
    "res_pval1", "res_pval2", "res_pval3"
  )

  ## data loading
  p <- seq(1, 1999, 2)

  P <- 8 # define number of core
  cl <- makeCluster(P) # start parallele computing

  clusterExport(cl, c("write_pvalue", "res_df"))

  system.time(
    res_par <- clusterApply(cl, p, fun = parralele_function) # eval parralele_function accross k
  )

  stopCluster(cl)

  ### Save results
  class(res_par)
  res_df <- bind_rows(res_par)
  tail(res_df)
}
# Saved results
res_df <- read.csv("../../datasets/results_df/hazard_ratio.csv")

res_df["names"] <- paste(res_df$k, res_df$iter, sep = "_")

df_res_HR <- res_df[, c("k", "iter", "res_HR1", "res_HR2", "res_HR3")]
df_res_HR_CI_low <- res_df[, c("k", "iter", "res_HR_CI_low1", "res_HR_CI_low2", "res_HR_CI_low3")]
df_res_HR_CI_high <- res_df[, c("k", "iter", "res_HR_CI_high1", "res_HR_CI_high2", "res_HR_CI_high3")]

#  compute confidence intervals
mean_hr1 <- res_df %>%
  group_by(k) %>%
  summarise(
    mean = mean(res_HR1),
    std = sqrt(var(res_HR1)),
    lower = mean(res_HR1) - qnorm(.975) * std / sqrt(n()),
    upper = mean(res_HR1) + qnorm(.975) * std / sqrt(n())
  )
mean_hr1$inf <- mean_hr1$mean - mean_hr1$std
mean_hr1$sup <- mean_hr1$mean + mean_hr1$std

poly_df <- rbind(
  setNames(mean_hr1[, c("k", "inf")], c("x", "y")),
  setNames(mean_hr1[order(nrow(mean_hr1):1), c("k", "sup")], c("x", "y"))
)

options(repr.plot.width = 10, repr.plot.height = 7)

g1 <- ggplot(res_df) +
  # red line as a ploygone
  stat_smooth(aes(k, res_HR1),
    method = "lm", se = FALSE,
    formula = y ~ poly(x, 8, raw = TRUE), colour = "red"
  ) +
  # add upper and lower bounds
  stat_smooth(
    data = mean_hr1, aes(k, mean - std), col = "#696969", lty = 2, se = FALSE,
    method = "gam", formula = y ~ s(x, bs = "cs")
  ) +
  stat_smooth(
    data = mean_hr1, aes(k, mean + std), col = "#696969", lty = 2, se = FALSE,
    method = "gam", formula = y ~ s(x, bs = "cs")
  )

gg1 <- ggplot_build(g1)

# get stat smooth (2,3) coordinate in order to draw the polygone between the smooth curve
df2 <- data.frame(x = gg1$data[[2]]$x, ymin = gg1$data[[2]]$y, ymax = gg1$data[[3]]$y)
df3 <- data.frame(x = gg1$data[[2]]$x, ymin = 0.39, ymax = 0.63)

# draw final figure
#  we must redo all the code lines to manage the layer order (point on the grey zone)

plotCc <- ggplot() +
  # original zone
  geom_ribbon(
    data = df3, aes(x = x, ymin = ymin, ymax = ymax),
    fill = colors["original", "color"], alpha = 0.07
  ) +
  # add original mean
  geom_segment(aes(x = 0, y = 0.49, xend = 1175, yend = 0.49),
    col = colors["original", "color"], lwd = 1, lty = 1, alpha = 0.8
  ) +

  # avatar zone
  geom_ribbon(
    data = df2, aes(x = x, ymin = ymin, ymax = ymax),
    fill = colors["avatar", "color"], alpha = 0.12
  ) +
  theme_minimal() +

  # add mean point by k
  geom_point(
    data = mean_hr1, aes(x = k, y = mean),
    colour = "white", alpha = 0.9, size = 2.5
  ) +
  geom_point(
    data = mean_hr1, aes(x = k, y = mean),
    colour = colors["avatar", "color"], alpha = 0.9, size = 2
  ) +
  # avatar mean smooth
  stat_smooth(
    data = res_df, aes(k, res_HR1), method = "lm", se = FALSE,
    formula = y ~ poly(x, 8, raw = TRUE), colour = "#35b696", size = 1.5
  ) +
  stat_smooth(
    data = res_df, aes(k, res_HR1), method = "lm", se = FALSE,
    formula = y ~ poly(x, 8, raw = TRUE), colour = colors["avatar", "color"], size = 1
  ) +

  # labs theme and axis
  ylab("Hazard ratio between\narms 0 and arms 1") +
  xlab("Number of nearest neighbors (k)") +
  ylim(0.35, 1) +
  # xlim(4, )
  scale_x_continuous(breaks = seq_k) +
  coord_cartesian(xlim = c(0, 1075)) +
  theme(
    axis.title.y = element_text(angle = 90, vjust = 1, size = axis_title_size, family = ""),
    axis.text = element_text(size = axis_text_size),
    axis.title.x = element_text(size = axis_title_size, family = ""),
    legend.text = element_text(size = legend_text_size, family = ""),
    legend.key.size = unit(1, "cm"),
    axis.line = element_line(
      colour = "black", size = 0.5, linetype = "solid",
      arrow = arrow(type = "closed", length = unit(5, "pt"))
    )
  )

# ggsave(file = "../../figure/aids_hazardratio_variousk.svg", plot = plot, width = 10, height = 7)

# If you never run this cell, it will download the necessary packages,
# It could take several minutes to download them

df_local_cloaking_k100 <- data.frame(row.names = 1:nrow(data))
hidden_rate_list <- c()


# coord_original <- pr$fit_transform(data, nf = 5L)[[1]]
coord_original <- FAMD(data, ncp = 5, graph = FALSE)
coord_original <- as.data.frame(coord_original$ind$coord)
for (k in seq_k) {
  iter <- 0

  # Extract dataset if necessary
  tryCatch(
    {
      avatar_tot <- read.table(paste0(file_data, "tmp/aids_avatarized/k", k, "_aids_avatarized.csv"),
        sep = ",", header = TRUE
      )
    },
    error = function(cnd) {
      print("The dataset was not extracted")
      untar(paste0(file_data, "AIDS/aids_avatarized.tar.bz2"),
        files = paste0("./aids_avatarized/k", k, "_aids_avatarized.csv"),
        exdir = paste0(file_data, "tmp/")
      )
      avatar_tot <- read.table(paste0(file_data, "tmp/aids_avatarized/k", k, "_aids_avatarized.csv"),
        sep = ",", header = TRUE
      )
    }
  )
  #  Subset and imput the avatar dataset
  avatar <- avatar_tot[avatar_tot$iter_k == iter, ]
  avatar_imp <- complete(mice(avatar, m = 1, maxit = 50, seed = 500, printFlag = FALSE), 1)
  avatar$cd496 <- avatar_imp$cd496

  metrics_tmp <- py$SecurityMetrics()
  metrics_tmp$fit(data, avatar[, 1:(ncol(avatar) - 3)], nf = 5L)

  hit_counts <- r_to_py(metrics_tmp$local_cloaking$hit_counts)
  avatar_local_cloaking <- c()

  for (i in seq_along(hit_counts)) {
    avatar_local_cloaking <- c(avatar_local_cloaking, py_to_r(hit_counts[i - 1])[1])
  }

  df_local_cloaking_k100[paste0("local_cloaking_", k)] <- avatar_local_cloaking
}

df_melt <- melt(df_local_cloaking_k100)


options(repr.plot.width = 10, repr.plot.height = 7)
plotCe <- ggplot(df_melt, aes(x = variable, y = value)) +
  geom_violin(trim = TRUE, alpha = 0.12, fill = colors["avatar", "color"]) +
  geom_boxplot(width = 0.4, fill = colors["avatar", "color"], outlier.size = 0.1, alpha = 0.8) +
  scale_x_discrete(labels = as.character(seq_k)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), breaks = c(1, 5, 10, 25, 50, 100, 300, 600, 1000, 2000)) +
  xlab("Number of nearest neighbors (k)") +
  ylab("Local cloaking") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.y = element_text(angle = 90, vjust = 1, size = axis_title_size),
    axis.text = element_text(size = axis_text_size),
    axis.title.x = element_text(size = axis_title_size),
    axis.line = element_line(
      colour = "black", size = 0.5, linetype = "solid",
      arrow = arrow(type = "closed", length = unit(5, "pt"))
    )
  ) +
  coord_cartesian(ylim = c(-1, 2000))

plotCe

# ggsave(file="../../figure/aids_localCloaking_violinplot.svg", plot = plot, width = 10, height = 7)


print("execution is done")
