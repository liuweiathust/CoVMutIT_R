suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggpmisc))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(bignum))

# println <- function(line) {cat(line); cat("\n")}
# 
# mutation_freq_table <- readr::read_rds("data/predict_example_data.rds")
# 
# colnames(mutation_freq_table) <- c("mutation", "position", "freq_prev", "freq_next")
# 
# #%% subsample 1% of the example data to accelerate the prediction process
# mutation_freq_table <- mutation_freq_table[sample(nrow(mutation_freq_table), ceiling(nrow(mutation_freq_table) / 10)), ]
# 
# mutation_freq_table %<>% filter(freq_prev != 0)

#%%
BN_target_score_fun <- function(freq_table, F) {
  mut_freq_BN <- freq_table %>%
    filter(freq_prev != 0) %>%
    mutate(shape1 = (1 - F) / F * freq_prev, shape2 = (1 - F) / F * (1 - freq_next)) %>%
    mutate(p = dbeta(freq_next, shape1, shape2)) %>%
    filter(!is.infinite(p)) %>%
    filter(p != 0)
  score <- mut_freq_BN %>% pull(p) %>% bignum::as_bigfloat() %>%  prod()
  return (score)
}

#%%
# Find the index of maximum  number from a list of bigfloat
get_max_bigfloat_index <- function(bigfloat_list) {
  max_score <- bignum::as_bigfloat(0)
  max_index <- -1
  for (n in seq_along(bigfloat_list)) {
    if (bignum::as_bigfloat(bigfloat_list[n]) > max_score) {
      max_score <- bignum::as_bigfloat(bigfloat_list[n])
      max_index <- n
    }
  }
  return (max_index)
}

#%%
generate_shrink_factor_for_bignum_score <- function(score) {
  shrink_factor <- ifelse(score >= 1e300, score / 1e300, score)
  return (shrink_factor)
}

#%%
# 
target_func_label <- latex2exp::TeX('$target = \\prod(dbeta(P_{b},   \\frac{1 - F}{F} \\cdot P_{a},   \\frac{1 - F}{F} \\cdot (1 - P_{a}) ))$')

covmutit_predict <- function(mutation_freq_table, session = NULL) {
  
  N <- nrow(mutation_freq_table)
  
  current_step = 0

  # 1st iteration ---------------------------------------------------------------------------------------------------
  x1 <- seq(1e-5, 1, length=1000)
  y1 <- rep(0, length(x1))
  z1 <- rep(0, length(x1))
  
  for (i in seq_along(x1)) {
    F <- x1[i]
    y1[i] <- BN_target_score_fun(freq_table=mutation_freq_table, F=F)
    current_step <- current_step + 1
    shinyWidgets::updateProgressBar(session = session, id = "predict__progress_bar", value = current_step, total = 4000, title = "Predicting (Iteration: 1) ...")
  }
  
  F1 <- x1[get_max_bigfloat_index(y1)]
  max_score1 <- BN_target_score_fun(freq_table=mutation_freq_table, F=F1)
  
  shrink_factor1 <- generate_shrink_factor_for_bignum_score(max_score1)
  
  for (i in seq_along(x1)) {
    z1[i] <- as.numeric(bignum::as_bigfloat(y1[i]) / bignum::as_bigfloat(shrink_factor1))
    current_step <- current_step + 1
    shinyWidgets::updateProgressBar(session = session, id = "predict__progress_bar", value = current_step, total = 4000, title = "Predicting (Iteration: 1) ...")
  }
  
  
  # 2nd iteration ---------------------------------------------------------------------------------------------------
  x2 <- seq(max(F1-0.0005, 1e-6), min(F1+0.0005, 1), length=1000)
  y2 <- rep(0, length(x2))
  z2 <- rep(0, length(x2))
  
  for (i in seq_along(x2)) {
    F = x2[i]
    y2[i] <- BN_target_score_fun(freq_table=mutation_freq_table, F=F)
    current_step <- current_step + 1
    shinyWidgets::updateProgressBar(session = session, id = "predict__progress_bar", value = current_step, total = 4000, title = "Predicting (Iteration: 2) ...")
  }
  
  F2 <- x2[get_max_bigfloat_index(y2)]
  
  max_score2 <- BN_target_score_fun(freq_table=mutation_freq_table, F=F2)
  
  shrink_factor2 <- generate_shrink_factor_for_bignum_score(max_score2)
  
  for (i in seq_along(x2)) {
    z2[i] <- as.numeric(bignum::as_bigfloat(y2[i]) / bignum::as_bigfloat(shrink_factor2))
    current_step <- current_step + 1
    shinyWidgets::updateProgressBar(session = session, id = "predict__progress_bar", value = current_step, total = 4000, title = "Predicting (Iteration: 2) ...")
  }
  
  F <- F2
  
  p1 <- tibble(x = x1, y = z1) %>%
    ggplot(aes(x = x, y = y)) +
    geom_line() +
    scale_y_continuous(breaks = NULL) +
    ylab(target_func_label) +
    xlab("F") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  p2 <- tibble(x = x2, y = z2) %>%
    ggplot(aes(x = x, y = y)) +
    geom_line() +
    scale_y_continuous(breaks = NULL) +
    scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
    geom_vline(xintercept = F, size = 0.3, color = "gray", linetype = "dashed") +
    ylab("") +
    xlab("") +
    annotation_custom(
      grid::textGrob(paste("F", "=", signif(F, 3), sep = " "), gp = grid::gpar(col = "steelblue", fontsize = 9)),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, vjust=0.5, size=8),
    ) 
  
  F_estimate_plot <- ggdraw(p1) + draw_plot(p2 , 0.47, 0.47, 0.5, 0.5) 
  
  pvalue_table <- mutation_freq_table %>%
    filter(freq_prev != 0) %>%
    mutate(shape1=(1-F)/F*freq_prev, shape2=(1-F)/F*(1-freq_prev)) %>%
    mutate(pvalue=1-pbeta(freq_next, shape1, shape2)) %>%
    arrange(pvalue)
  
  freq_scatter_plot <- ggplot() +
    geom_abline(slope=1, intercept = 0, color="#CCCCCC", linetype=2, size=1) +
    geom_point(aes(x=freq_prev, y=freq_next), data=pvalue_table, size=2, alpha=0.8, color="#999999") +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    xlab("mutation frequency (prev-month)") +
    ylab("mutation frequency (next-month)") +
    coord_fixed() 
  
  max_value <- 20
  
  bn_manhattan_data <- pvalue_table %>% mutate(log10pvalue=ifelse(pvalue == 0, max_value, -log10(pvalue)))
  
  pvalue_manhattan_plot <- ggplot() +
    geom_vline(xintercept=21563, color="#CCCCCC", size=1.0, alpha=0.8, linetype=2) +
    geom_vline(xintercept=25384, color="#CCCCCC", size=1.0, alpha=0.8, linetype=2) +
    geom_point(aes(x=position, y=log10pvalue), size=2, color="#999999", alpha=0.85, data=bn_manhattan_data) +
    scale_y_continuous(limits = c(0, max_value)) +
    ylab(expression("-log"[10]~"(p-value)")) +
    xlab("SARS-CoV-2 genome position") +
    theme(
      legend.position = "none"
    )
  
  freq_change_density_plot <- mutation_freq_table %>% 
    mutate(freq_diff = freq_next - freq_prev) %>% 
    ggplot(aes(x=freq_diff)) +
    geom_density(color = DEFAULT_COLOR_PAL[1]) +
    scale_x_continuous(limits = c(-1, 1)) +
    xlab("freq_next - freq_prev") +
    ylab("Density")
  
  freq_density_plot <- mutation_freq_table %>% gather(group, freq, c(freq_prev, freq_next)) %>% 
    ggplot(aes(x=freq, color=group)) +
    geom_density(alpha = 0.4) +
    scale_fill_aaas() +
    scale_color_aaas() +
    facet_wrap(group ~ ., nrow = 2)

  return(list(
    F_estimate_plot = F_estimate_plot, 
    table = pvalue_table, 
    scatter_plot = freq_scatter_plot, 
    manhattan_plot = pvalue_manhattan_plot,
    freq_density_plot = freq_density_plot,
    freq_change_density_plot = freq_change_density_plot
  ))
}



