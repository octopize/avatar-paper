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
  # Save corrplot as svg
  svglite,
  # Display a Sankey diagram
  networkD3,
  htmlwidgets,
  forcats,
  stringr,
  # data table management
  tibble,
  dplyr
)

# import specific function to compute metrics
py_run_string("import sys")
py_run_string("sys.path.append('..')")

py_run_string("from metrics.privacy_metrics import SecurityMetrics")
py_run_string("from metrics.privacy_metrics.record_to_avatar_distance import record_to_avatar_distance")
py_run_string("from metrics.privacy_metrics.local_cloaking import get_local_cloaking")

save=FALSE

## utils function
categorical <- c(
  "hemo", "homo", "drugs", "karnof", "oprior", "z30", "zprior", "race",
  "gender", "str2", "strat", "symptom", "treat", "offtrt", "r", "cens", "arms"
)

apply_dtypes <- function(data) {
  data[categorical] <- lapply(data[categorical], factor)
  return(data)
}

file_data <- "datasets/"
data <- read.csv(paste0("../datasets/AIDS/aids_original_data.csv"), sep = ";", na.strings = c("NA", "", NA))
data <- data[, -1]

# Avatar data with k=20
avatar <- read.csv(paste0("../datasets/AIDS/aids_avatarized_base_k20_nf5.csv"), sep = ",", na.strings = c("NA", "", NA))
avatar$arms <- as.factor(avatar$arms)

# Synthpop data
synthpop <- read.csv(paste0("../datasets/AIDS/aids_synthpop_base.csv"), sep = ",", na.strings = c("NA", "", NA))
synthpop$arms <- as.factor(synthpop$arms)

# CTGAN data
ctgan <- read.csv(paste0("../datasets/AIDS/aids_CTGAN_base_2.csv"), sep = ",", na.strings = c("NA", "", NA))
ctgan$arms <- as.factor(ctgan$arms)

colors <- read.csv("../color.csv", stringsAsFactors = FALSE)
rownames(colors) <- colors$type

axis_title_size <- 19
axis_text_size <- 19
legend_text_size <- 19
legend_title_size <- 19

# input missing data to project the individuals in the multidimensional space
data_imp <- complete(mice(data, m = 1, maxit = 50, seed = 500, printFlag = FALSE), 1)
avatar_imp <- complete(mice(avatar, m = 1, maxit = 50, seed = 500, printFlag = FALSE), 1)
synthpop_imp <- complete(mice(synthpop, m = 1, maxit = 50, seed = 500, printFlag = FALSE), 1)
ctgan_imp <- complete(mice(ctgan, m = 1, maxit = 50, seed = 500, printFlag = FALSE), 1)

data$cd496 <- data_imp$cd496
avatar$cd496 <- avatar_imp$cd496
synthpop$cd496 <- synthpop_imp$cd496
ctgan$cd496 <- ctgan_imp$cd496

##### ORIGINAL - SYNTHETIC projection analysis

get_aids_2D_projection <- function(data, synthetic_name = "synthetic") {
  data <- apply_dtypes(data)

  famd <- FAMD(data, ncp = 5, graph = FALSE, ind.sup = (nrow(data) / 2 + 1):nrow(data))
  res_ind <- as.data.frame(famd$ind$coord)
  res_ind_sup <- as.data.frame(famd$ind.sup$coord)

  res_ind_2D <- res_ind[, 1:2]
  res_ind_2D["type"] <- "Original"
  res_ind_sup["type"] <- synthetic_name
  res_ind_2D <- rbind(res_ind_2D, res_ind_sup[, c("Dim.1", "Dim.2", "type")])

  set.seed(42)
  rows <- sample(nrow(res_ind_2D))
  res_ind_2D <- res_ind_2D[rows, ]
  return(list("coord" = res_ind_2D, "model" = famd))
}


get_plot_projection <- function(projection, name, save = FALSE) {
  options(repr.plot.width = 10, repr.plot.height = 7)
  plot <- ggplot(projection$coord, aes(x = Dim.1, y = Dim.2, fill = factor(type))) +
    # point
    geom_point(size = 3, shape = 21, alpha = 1) +
    # fill by type
    aes(fill = factor(type)) +
    scale_fill_manual(values = c(colors[name, "color"], colors["original", "color"])) +
    # theme and details
    xlab(paste0("Dim. 1 (", formatC(projection$model$eig[1, 2], format = "f", digits = 2), "%)")) +
    ylab(paste0("Dim. 2 (", formatC(projection$model$eig[2, 2], format = "f", digits = 2), "%)")) +
    labs(fill = "") +
    theme_bw() +
    theme(
      legend.position = c(0.86, 0.15),
      legend.text = element_text(size = legend_text_size, color = "black", family = "sans"),
      axis.text = element_text(size = axis_text_size, color = "black", family = "sans"),
      axis.title = element_text(size = axis_title_size, color = "black", family = "sans"),
      legend.background = element_rect(fill = "white", linetype = "solid")
    )

  if (save) {
    ggsave(file = paste0("../figures/aids_pca2D_", name, ".svg"), plot = plot, width = 10, height = 7, dpi = 320)
  }
  return(plot)
}

# avatar
data_avatar <- rbind(data, avatar)
avatar_projection <- get_aids_2D_projection(data_avatar, synthetic_name = "Avatar")
plotAa <- get_plot_projection(avatar_projection, name = "avatar", save = save)


# synthpop
data_synthpop <- rbind(data, synthpop)
synthpop_projection <- get_aids_2D_projection(data_synthpop, synthetic_name = "Synthpop")
plotAa_synthpop <- get_plot_projection(synthpop_projection, name = "synthpop", save = save)

# ctgan
data_ctgan <- rbind(data, ctgan)
ctgan_projection <- get_aids_2D_projection(data_ctgan, synthetic_name = "CT-GAN")
plotAa_ctgan <- get_plot_projection(ctgan_projection, name = "ctgan", save = save)


##### ORIGINAL - SYNTHETIC survival curve analysis

get_survival_results <- function(data, synthetic, names = c("Original", "Synthetic")) {
  data$arms <- as.factor(data$arms)

  data_arms01 <- data[data$arms %in% c(0, 1), ]
  data_arms01$arms <- as.factor(as.character(data_arms01$arms))

  summary_cox <- summary(coxph(Surv(time = days / 7, event = cens) ~ arms, data = data_arms01))
  res_original_01 <- cbind(summary_cox$coefficients, summary_cox$conf.int)

  summary_cox <- summary(coxph(Surv(time = days / 7, event = cens) ~ arms, data = synthetic))
  res_avat <- cbind(summary_cox$coefficients[, c("exp(coef)", "Pr(>|z|)")], summary_cox$conf.int[, c("lower .95", "upper .95")])

  synthetic_arms01 <- synthetic[synthetic$arms %in% c(0, 1), ]
  synthetic_arms01$arms <- as.factor(as.character(synthetic_arms01$arms))
  summary_cox <- summary(coxph(Surv(time = days / 7, event = cens) ~ arms, data = synthetic_arms01))
  res_synthetic_01 <- cbind(summary_cox$coefficients, summary_cox$conf.int)

  data_arms01$type <- names[1]
  synthetic_arms01$type <- names[2]
  data_concat <- rbind(data_arms01, synthetic_arms01)
  data_concat$type <- factor(data_concat$type, levels = names)

  # Survival curve computation
  survival_curve <- survfit(Surv(time = days / 7, event = cens) ~ factor(arms) + type, data = data_concat)

  return(list("data" = data_concat, "curve" = survival_curve, "hazard_ratio_synthetic" = res_synthetic_01, "hazard_ratio_original" = res_original_01))
}


get_survival_plot <- function(survival_curve, survival_data, hazard_ratio_synthetic, hazard_ratio_original, names = c("Original", "Synthetic"), save=FALSE) {
  color <- c("#7155c8", "#c85573")

  options(repr.plot.width = 10, repr.plot.height = 7)
  ggsurv <- ggsurvplot(survival_curve,
    data = survival_data, conf.int = FALSE,
    risk.table = FALSE,
    legend.title = "Treatments",
    censor.cex = 100,
    linetype = c("type"),
    size = 0.7, xlim = c(0, 155), ylim = c(0.5, 1), break.x.by = 40,
    ggtheme = theme_minimal(),
    xlab = "Time (week)",
    ylab = "Proportion not reaching\nprimary end point",
    censor.shape = "",
    legend.labs = c("Arm 0", "Arm 0", "Arm 1", "Arm 1"),
    palette = color
  )

  survival_plot <- ggsurv$plot <- ggsurv$plot +
    theme_minimal() +
    scale_linetype_discrete(name = c("Data source"), labels = names) +
    geom_label(aes(x = 11, y = 0.625, label = "Hazard ratio"), size = 6, family = "sans", fontface = "bold", label.size = NA) +
    geom_label(aes(x = 59.5, y = 0.58, label = paste0(
      "Original: HR [CI: 95%] = ", formatC(hazard_ratio_original[2], format = "f", digits = 2),
      " [", formatC(hazard_ratio_original[8], format = "f", digits = 2),
      "-", formatC(hazard_ratio_original[9], format = "f", digits = 2),
      "] ; p-value = ", formatC(hazard_ratio_original[5], format = "e", digits = 2),
      paste0("\n", names[1], ": HR [CI: 95%] = "), formatC(hazard_ratio_synthetic[2], format = "f", digits = 2),
      " [", formatC(hazard_ratio_synthetic[8], format = "f", digits = 2),
      "-", formatC(hazard_ratio_synthetic[9], format = "f", digits = 2),
      "] ; p-value = ", formatC(hazard_ratio_synthetic[5], format = "e", digits = 2), "  "
    )),
    size = 6.2, family = "sans", label.size = NA
    ) +
    theme(
      legend.position = c(0.105, 0.5),
      legend.background = element_rect(fill = "white", size = 0.5, linetype = "blank"),
      legend.text = element_text(size = legend_text_size, color = "black", family = "sans"),
      legend.title = element_text(size = legend_title_size, color = "black", family = "sans", face = "bold"),
      axis.text = element_text(size = axis_text_size, color = "black", family = "sans"),
      axis.title.y = element_text(vjust = 2),
      axis.title = element_text(size = axis_title_size, color = "black", family = "sans"),
      axis.line = element_line(colour = "black", size = 0.5, linetype = "solid", arrow = arrow(type = "closed", length = unit(5, "pt"))),
      text = element_text(size = 14)
    )
  if (save){
    ggsave(file = paste0("../figures/tmp_aids_", names[2], "_survival.svg"), plot = survival_plot, width = 10, height = 7, dpi = 290)
  }

  return(survival_plot)
}


### Apply survival analysis

# Avatar
survival_avatar <- get_survival_results(data = data, synthetic = avatar, names = c("Original", "Avatar"))
# Display time: about 15 seconds
survival_avatar_data <- survival_avatar$data
survival_avatar_curve <- survival_avatar$curve
hazard_ratio_avatar <- survival_avatar$hazard_ratio_synthetic
hazard_ratio_original <- survival_avatar$hazard_ratio_original

plotAc <- get_survival_plot(
  survival_curve = survival_avatar_curve,
  survival_data = survival_avatar_data,
  hazard_ratio_synthetic = hazard_ratio_avatar,
  hazard_ratio_original = hazard_ratio_original,
  names = c("Original", "Avatar")
)

# Synthpop
survival_synthpop <- get_survival_results(data = data, synthetic = synthpop, names = c("Original", "Synthpop"))
# Display time: about 15 seconds
survival_data <- survival_synthpop$data
survival_curve <- survival_synthpop$curve
hazard_ratio_synthpop <- survival_synthpop$hazard_ratio_synthetic
hazard_ratio_original <- survival_synthpop$hazard_ratio_original

survival_plot_synthpop <- get_survival_plot(
  survival_curve = survival_curve,
  survival_data = survival_data,
  hazard_ratio_synthetic = hazard_ratio_synthpop, 
  hazard_ratio_original = hazard_ratio_original
)

# ctgan
survival_ctgan <- get_survival_results(data = data, synthetic = ctgan, names = c("Original", "CT-GAN"))
survival_data <- survival_ctgan$data
survival_curve <- survival_ctgan$curve
hazard_ratio_ctgan <- survival_ctgan$hazard_ratio_synthetic
hazard_ratio_original <- survival_ctgan$hazard_ratio_original

survival_plot_ctgan <- get_survival_plot(
  survival_curve = survival_curve,
  survival_data = survival_data,
  hazard_ratio_synthetic = hazard_ratio_ctgan, 
  hazard_ratio_original = hazard_ratio_original
)

## Comparative hazard ratio analysis
Hazard_Ratio <- c(hazard_ratio_original[2], hazard_ratio_avatar[2], hazard_ratio_synthpop[2], hazard_ratio_ctgan[2])
ci_high <- c(hazard_ratio_original[9], hazard_ratio_avatar[9], hazard_ratio_synthpop[9], hazard_ratio_ctgan[9])
ci_low <- c(hazard_ratio_original[8], hazard_ratio_avatar[8], hazard_ratio_synthpop[8], hazard_ratio_ctgan[8])
p_values <- c(hazard_ratio_original[5], hazard_ratio_avatar[5], hazard_ratio_synthpop[5], hazard_ratio_ctgan[5])

p_values_text <- paste0('p-value: ',  format.pval(p_values, digits = 3))
type <- c('Original', 'Avatar', 'Synthpop', 'CT-GAN')

aids_results <- data.frame(Hazard_Ratio, ci_high, ci_low, type, p_values_text)
aids_results$type_f = factor(aids_results$type, levels=c('CT-GAN', 'Synthpop', 'Avatar','Original'))

aids_comparative_plot <- ggplot(aids_results, aes(x = Hazard_Ratio, y = type_f, colour = type)) +
  geom_errorbarh(aes(xmax = ci_high, xmin = ci_low), size = 0.6, height
                 = 0.15) +
  geom_text(
      aes(label = p_values_text, x = as.numeric(Hazard_Ratio)), 
      hjust = 0.5,
      vjust = -2,
      color = "black", 
      family = "sans", 
      size = 4.5,
           ) +
  geom_point(size = 2) +
  scale_colour_manual(values = c('Original'= colors['original', 'color'],
                                 'Avatar'=colors['avatar', 'color'],
                                 'Synthpop'=colors['synthpop', 'color'] ,
                                 'CT-GAN'=colors['ctgan', 'color']),
                      breaks=c('Original', 'Avatar', 'Synthpop', 'CT-GAN')) +
  theme_bw() +
  xlim(c(0.1,0.8)) +
  ylab(NULL) +
  xlab("Hazard Ratio (95% CI)") +
  theme(
      text = element_text(size = 16),
      legend.position="none",
      legend.text = element_text(size = legend_text_size, color = "black", family = "sans"),
      axis.text = element_text(size = axis_text_size, color = "black", family = "sans"),
      axis.title = element_text(size = axis_title_size, color = "black", family = "sans"),
      legend.background = element_rect(fill = "white", linetype = "solid")
  )

if (save){
  ggsave(file = "../figures/aids_comparative_plot.svg", plot = aids_comparative_plot, width = 10, height = 7, dpi = 290)
}


## Supplementary graph : Arms 1-2-3-4 for avatar comparison
data_typed <- data.frame(data)
avatar_typed <- data.frame(avatar)
data_typed$type <- "Original"
avatar_typed$type <- "Avatar"
data_concat <- rbind(data_typed, avatar_typed)
data_concat$type <- factor(data_concat$type, levels = c("Original", "Avatar"))

# Survival curve calculation
surv <- survfit(Surv(time = days / 7, event = cens) ~ factor(arms) + type, data = data_concat)

color <- c("#7155c8", "#c85573", "#95aa5b", "#c37e47")

options(repr.plot.width = 12, repr.plot.height = 8)
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
  legend.labs = c("Arm 0", "Arm 0", "Arm 1", "Arm 1", "Arm 2", "Arm 2", "Arm 3", "Arm 3"),
  palette = color
)

ggsurv$plot <- ggsurv$plot +
  theme_minimal() +
  theme(
    legend.position = c(0.105, 0.5),
    legend.background = element_rect(fill = "white", size = 0.5, linetype = "blank"),
    legend.text = element_text(size = legend_text_size, color = "black", family = "sans"),
    axis.text = element_text(size = axis_text_size, color = "black", family = "sans"),
    axis.title = element_text(size = axis_title_size, color = "black", family = "sans", face = "bold"),
    axis.title.y = element_text(vjust = 2),
    legend.title = element_text(size = legend_title_size, face = "bold", family = "sans"),
    axis.line = element_line(colour = "black", size = 0.5, linetype = "solid", arrow = arrow(type = "closed", length = unit(5, "pt")))
  ) +
  scale_linetype_discrete(name = "Data source")

ggsurv


# Hazard ratio table  (only avatar dataset)
data$arms <- as.factor(as.character(data$arms))
avatar$arms <- as.factor(as.character(avatar$arms))
summary_cox_original <- summary(coxph(Surv(time = days / 7, event = cens) ~ arms, data = data))
res_original <- cbind(summary_cox_original$coefficients[, c("exp(coef)", "Pr(>|z|)")], summary_cox_original$conf.int[, c("lower .95", "upper .95")])
summary_cox_avatar <- summary(coxph(Surv(time = days / 7, event = cens) ~ arms, data = avatar))
res_avatar <- cbind(summary_cox_avatar$coefficients[, c("exp(coef)", "Pr(>|z|)")], summary_cox_avatar$conf.int[, c("lower .95", "upper .95")])


rownames(res_original) <- paste0(rownames(res_original), "_original")
rownames(res_avatar) <- paste0(rownames(res_avatar), "_avatar")
res_concat <- as.data.frame(rbind(res_original, res_avatar))
res_concat["Hazard Ratio"] <- formatC(res_concat[, 1], format = "f", digits = 2)
res_concat["95% Confidence Interval"] <- paste0(formatC(res_concat[, 3], format = "f", digits = 2), "-", formatC(res_concat[, 4], format = "f", digits = 2))
res_concat["p-value"] <- formatC(res_concat[, 2], format = "e", digits = 2)
res_concat <- res_concat[c(1, 4, 2, 5, 3, 6), 5:7]
row.names(res_concat) <- c("Original arm 1", "Avatar arm 1", "Original arm 2", "Avatar arm 2", "Original arm 3", "Avatar arm 3")


### PRIVACY METRICS for avatar data

metrics_aids <- py$SecurityMetrics()
data <- apply_dtypes(data)
avatar <- apply_dtypes(avatar)
metrics_aids$fit(data, avatar, nf = 5L)

hit_counts <- r_to_py(metrics_aids$local_cloaking$hit_counts)
avatar_local_cloaking <- c()
original_local_cloaking <- c()

for (i in seq_along(hit_counts)) {
  avatar_local_cloaking <- c(avatar_local_cloaking, py_to_r(hit_counts[i - 1])[1])
  original_local_cloaking <- c(original_local_cloaking, py_to_r(hit_counts[i - 1])[2])
}
df_ref_local_cloaking_aids <- data.frame(avatar = avatar_local_cloaking, original = original_local_cloaking)
color_vector <- c(rep("#c14f4f", 2), rep(colors["avatar", "color"], 98))

axis_title_size <- 21
axis_text_size <- 22
legend_text_size <- 19
legend_title_size <- 19
intersect_median_density_aids <- 52

## plot privacy metrics
options(repr.plot.width = 10, repr.plot.height = 7)
plotBa <- ggplot(df_ref_local_cloaking_aids, aes(avatar)) +
  # add histogram
  geom_histogram(bins = 100, fill = color_vector, color = "#1D1D1B") +
  #  add density
  geom_density(aes(y = ..count..), lwd = 1.5, color = "#1D1D1B") +
  # add median
  geom_segment(aes(x = median(avatar), y = 0, xend = median(avatar), yend = intersect_median_density_aids), color = "#1D1D1B", size = 2) +
  # theme, axis, legend
  theme_minimal() +
  xlab("Local cloaking") +
  ylab("Number of individuals") +
  theme(
    axis.title = element_text(size = axis_title_size, family = ""),
    axis.text = element_text(size = axis_title_size, family = ""),
    axis.line = element_line(colour = "black", size = 0.5, linetype = "solid", arrow = arrow(type = "closed", length = unit(5, "pt")))
  ) +
  xlim(c(-1, 100)) +

  # add arrow and text
  ## arrow median
  geom_segment(aes(
    x = median(avatar) + 25,
    y = intersect_median_density_aids + 50,
    xend = median(avatar) + 1,
    yend = intersect_median_density_aids
  ),
  color = "#1D1D1B", size = 0.8,
  arrow = arrow(type = "closed", length = unit(10, "pt"))
  ) +
  ## text median
  geom_label(aes(
    x = median(avatar) + 39.,
    y = intersect_median_density_aids + 50,
    label = paste0("Median = ", median(avatar))
  ),
  size = 9,
  family = "sans",
  label.size = NA
  ) +
  ## arrow hidden rate
  geom_segment(aes(
    x = median(avatar) + 25,
    y = intersect_median_density_aids + 88,
    xend = 1, yend = sum(df_ref_local_cloaking_aids[, 1] == 0)
  ),
  color = "#1D1D1B", size = 0.8,
  arrow = arrow(
    type = "closed",
    length = unit(10, "pt")
  )
  ) +
  ## text hidden rate
  geom_label(aes(
    x = median(avatar) + 45,
    y = intersect_median_density_aids + 88,
    label = paste0("Hidden rate = ", round(metrics_aids$hidden_rate, 0), " %")
  ),
  size = 9,
  label.size = NA
  ) +
  xlim(c(-1, 100))

if (save){
  ggsave(file="../figures/aids_local_cloaking.svg", plot=plot, width=10, height=7, dpi = 290)
}
 
# Proportion of individuals with a local cloaking under 5
prop_aids <- sum(df_ref_local_cloaking_aids$avatar <= 5) / dim(df_ref_local_cloaking_aids)[1] * 100

avatar_tot <- read.csv("../datasets/AIDS/aids_avatarized_25time_k20_nf_5.csv")
avatar_tot[categorical] <- lapply(avatar_tot[categorical], factor)
data[categorical] <- lapply(data[categorical], factor)
# Projection
coord_original <- FAMD(data[, -ncol(data)], ncp = 5, graph = FALSE)
coord_original <- as.data.frame(coord_original$ind$coord)
coord_original <- reticulate::r_to_py(as.data.frame(coord_original))


# Computation time: about 2 minutes

df_local_cloaking_aids <- data.frame(row.names = 1:nrow(data))
for (k in unique(avatar_tot$iter)) {
  avatar <- avatar_tot[avatar_tot$iter == k, ]
  avatar_imp <- complete(mice(avatar, m = 1, maxit = 50, seed = 500, printFlag = FALSE), 1)
  avatar$cd496 <- avatar_imp$cd496

  metrics <- py$SecurityMetrics()
  metrics$fit(data[, 1:(ncol(avatar) - 1)], avatar[, 1:(ncol(avatar) - 1)], nf = 5L)
  avatar_local_cloaking_tmp <- c()
  hit_counts_tmp <- r_to_py(metrics$local_cloaking$hit_counts)

  for (i in seq_along(hit_counts_tmp)) {
    avatar_local_cloaking_tmp <- c(avatar_local_cloaking_tmp, py_to_r(hit_counts_tmp[i - 1])[1])
  }
  df_local_cloaking_aids[paste0("local_cloaking_", k, "_", 20)] <- data.frame(avatar_local_cloaking_tmp)
}

df_local_cloaking_aids_LC0_total <- as.data.frame(table(rowSums(df_local_cloaking_aids == 0)))
df_local_cloaking_aids_LC0_cum <- tibble(df_local_cloaking_aids_LC0_total %>%
  mutate(
    Var1 = as.character(fct_collapse(
      Var1,
      "10" = as.character(10:max(as.numeric(df_local_cloaking_aids_LC0_total$Var1)))
    ))
  ) %>%
  group_by(Var1) %>%
  mutate(Freq = sum(Freq)) %>%
  distinct())

df_local_cloaking_aids_LC0_cum["Freq2"] <- paste0(
  formatC(df_local_cloaking_aids_LC0_cum$Freq / sum(df_local_cloaking_aids_LC0_cum$Freq) * 100, format = "f", digits = 1), "%"
)

df_local_cloaking_aids_LC0_cum[df_local_cloaking_aids_LC0_cum$Var1 == 10, "Var1"] <- "\U2265 10"
df_local_cloaking_aids_LC0_cum$Var1 <- factor(df_local_cloaking_aids_LC0_cum$Var1, levels = c(as.numeric(0:8), "\U2265 10"))

options(repr.plot.width = 10, repr.plot.height = 7)
plotBc <- ggplot(df_local_cloaking_aids_LC0_cum, aes(x = as.factor(Var1), y = Freq)) +
  geom_bar(stat = "identity", fill = colors["original", "color"], colour = "black", show.legend = FALSE) +
  geom_text(aes(label = Freq2), stat = "identity", vjust = -0.5, size = 6) +
  ylab("Number of individuals") +
  xlab("Number of local cloaking at zero over 25 avatarizations") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = axis_title_size, family = "sans"),
    axis.text = element_text(size = axis_text_size, family = "sans"),
    axis.line = element_line(
      colour = "black", size = 0.5,
      linetype = "solid",
      arrow = arrow(
        type = "closed",
        length = unit(5, "pt")
      )
    )
  ) +
  scale_x_discrete(
    breaks = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "\U2265 10"),
    limits = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "\U2265 10")
  )

if (save){
  ggsave(file="../figures/aids_zeroCloaking_percentage.svg", plot=plot, width=10, height=7, dpi = 320)
}

df_localCloaking_plot <- df_local_cloaking_aids
df_localCloaking_plot[df_localCloaking_plot >= 5] <- 5
df_localCloaking_plot["names"] <- as.numeric(rownames(df_localCloaking_plot))

freq_local_cloaking <- as.data.frame(table(melt(df_localCloaking_plot, id.vars = "names")$value))
freq_local_cloaking["Freq"] <- formatC(freq_local_cloaking[, "Freq"] / sum(freq_local_cloaking[, "Freq"]) * 100, format = "f", digits = 2)
freq_local_cloaking["Var1"] <- as.character(freq_local_cloaking[, "Var1"])
freq_local_cloaking[freq_local_cloaking["Var1"] == "5", "Var1"] <- ">= 5"
colnames(freq_local_cloaking) <- c("Local cloaking", "Probability (%)")
freq_local_cloaking

df_loc_6 <- df_local_cloaking_aids
df_loc_6[df_loc_6 >= 5] <- ">=5"
for (i in 1:ncol(df_loc_6)) {
  df_loc_6[, i] <- factor(df_loc_6[, i], levels = c(as.character(0:4), ">=5"))
}

table_recap <- data.frame()
for (i in 1:24) {
  if (length(table_recap) == 0) {
    table_recap <- table(df_loc_6[, i], df_loc_6[, i + 1])
  } else {
    table_recap <- table_recap + table(df_loc_6[, i], df_loc_6[, i + 1])
  }
}
table_recap <- as.matrix(round(table_recap / rowSums(table_recap), 3))
row.names(table_recap) <- paste0(c(as.character(0:4), "\U2265 5"), "  ")
colnames(table_recap) <- paste0(c(as.character(0:4), "\U2265 5"), "\n")


if (save) {
  svglite("../figures/aids_Probabilities.svg", width = 8, height = 8)
}
options(repr.plot.width = 10, repr.plot.height = 10)
corrplot <- corrplot(table_recap,
  method = "color", is.corr = FALSE, addCoef.col = "black", addCoefasPercent = TRUE, number.cex = 1.2, mar = c(7, 7, 7, 7), cl.pos = "n", tl.col = "black", tl.cex = 1.5, tl.srt = 0,
  col = colorRampPalette(c("white", "#f5f5f5", colors["avatar", "color"]))(50)
)
mtext(text = "Local cloaking first iteration", side = 2, line = 0, las = 0, cex = 1.9, adj = 0.18, outer = FALSE)
mtext(text = "Local cloaking second iteration", side = 3, line = 0, las = 0, cex = 1.8, adj = 0.30, outer = FALSE)

if (save) {
  dev.off()
}


df_table <- data.frame()
for (i in 1:4) {
  if (nrow(df_table) == 0) {
    df_table <- as.data.frame(table(df_loc_6[, i], df_loc_6[, i + 1]))
    df_table$avatar_dep <- paste0("A", i)
    df_table$avatar_fin <- paste0("A", i + 1)
  } else {
    temp <- as.data.frame(table(df_loc_6[, i], df_loc_6[, i + 1]))
    temp$avatar_dep <- paste0("A", i)
    temp$avatar_fin <- paste0("A", i + 1)
    df_table <- rbind(df_table, temp)
  }
}

df_table["deb"] <- str_replace(paste0(df_table$avatar_dep, " ", df_table$Var1), ">=", "\U2265")
df_table["fin"] <- str_replace(paste0(df_table$avatar_fin, " ", df_table$Var2), ">=", "\U2265")

df_table_graph <- data.frame(source = as.character(df_table$deb), target = factor(paste0(df_table$fin)), value = df_table$Freq, deb = df_table$Var1, fin = df_table$Var2)
df_table_graph <- df_table_graph[order(str_replace(df_table_graph$source, "\U2265", "")), ]

nodes <- data.frame(name = c(as.character(df_table_graph$source), as.character(df_table_graph$target)) %>% unique())
nodes$group <- sapply(nodes$name, function(x) {
  str_split(x, " ")[[1]][2]
})
df_table_graph$IDsource <- match(df_table_graph$source, nodes$name) - 1
df_table_graph$IDtarget <- match(df_table_graph$target, nodes$name) - 1

ColourScal <- "d3.scaleOrdinal() .range(["
nb <- 6
col <- colorRampPalette(c(colors["original", "color"], colors["avatar", "color"]))(nb)
for (i in 1:nb) {
  if (i != nb) {
    ColourScal <- paste0(ColourScal, "\"", col[i], "\"", ",")
  } else {
    ColourScal <- paste0(ColourScal, "\"", col[i], "\"", "])")
  }
}

options(repr.plot.width = 40, repr.plot.height = 7)

# Make the Network
sknetwork <- sankeyNetwork(
  Links = df_table_graph, Nodes = nodes,
  Source = "IDsource", Target = "IDtarget",
  Value = "value", NodeID = "name", NodeGroup = "group",
  sinksRight = FALSE, colourScale = ColourScal, nodeWidth = 15, fontSize = 13, nodePadding = 20, iteration = 0
)

if (save){
  saveWidget(sknetwork, file="../figures/sankeyBasic_count_A1-A5.html", selfcontained = FALSE)
}
 

print("execution is done")
