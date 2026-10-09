cat("\nEstimation complete.\n")
cat("Outputs saved to:\n", output_directory, "\n")
# Jason West
# Submitted for consideration to Journal of Energy Markets
#
# Corrected estimation script
# ============================================================

rm(list = ls())
graphics.off()

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

required_packages <- c(
  "data.table",
  "lubridate",
  "VineCopula",
  "strucchange",
  "ggplot2",
  "patchwork",
  "scales"
)

new_packages <- required_packages[
  !required_packages %in% installed.packages()[, "Package"]
]

if (length(new_packages) > 0) {
  install.packages(new_packages)
}

invisible(lapply(required_packages, library, character.only = TRUE))

set.seed(12345)

# ------------------------------------------------------------
# 2. User settings
# ------------------------------------------------------------
# ============================================================
# Regime-Aware Bivariate Copula Modelling 
# of Wind-Price Risk
# ============================================================

rm(list = ls())
graphics.off()

# ---- libraries ----
library(data.table)
library(lubridate)
library(fitdistrplus)
library(VineCopula)
library(copula)
library(ggplot2)
library(ggExtra)
library(viridis)

# ---- set dir ----
setwd('C:/Users/jwest.INTERNAL/OneDrive - Bureau of Meteorology/Documents/Data/Wind/Copulas')
data_directory <- 'C:/Users/jwest.INTERNAL/OneDrive - Bureau of Meteorology/Documents/Data/Wind/Copulas'

price_generation_file <- "pricegen.csv"
wind_file             <- "winds.csv"

# Rolling-window settings
rolling_years <- 2
update_months <- 1

# At least this proportion of rolling observations must lie
# within each Bai-Perron regime.
minimum_regime_fraction <- 0.15

# Copula families:
# 1  = Gaussian
# 2  = Student-t
# 3  = Clayton
# 4  = Gumbel
# 5  = Frank
# 6  = Joe
# 7  = BB1
# 8  = BB6
# 9  = BB7
# 10 = BB8
base_family_set <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

# Fixed empirical tail thresholds
tail_levels <- c(0.10, 0.05, 0.01)

# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

clip_uniform <- function(x, epsilon = 1e-6) {
  pmin(pmax(x, epsilon), 1 - epsilon)
}

make_pseudo_observations <- function(x) {
  clip_uniform(rank(x, ties.method = "average") / (length(x) + 1))
}

family_name_safe <- function(family) {
  tryCatch(
    BiCopName(family),
    error = function(e) paste("Family", family)
  )
}

parameter_count <- function(family) {
  if (family %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20,
                    27, 28, 29, 30, 37, 38, 39, 40)) {
    return(2L)
  }
  
  if (family == 0) {
    return(0L)
  }
  
  1L
}

extract_fit_statistics <- function(fit, model_name, n_obs) {
  
  k <- parameter_count(fit$family)
  
  log_likelihood <- if (!is.null(fit$logLik)) {
    as.numeric(fit$logLik)
  } else {
    sum(
      log(
        BiCopPDF(
          fit$u1,
          fit$u2,
          family = fit$family,
          par = fit$par,
          par2 = fit$par2
        )
      )
    )
  }
  
  aic_value <- if (!is.null(fit$AIC)) {
    as.numeric(fit$AIC)
  } else {
    -2 * log_likelihood + 2 * k
  }
  
  bic_value <- if (!is.null(fit$BIC)) {
    as.numeric(fit$BIC)
  } else {
    -2 * log_likelihood + log(n_obs) * k
  }
  
  tau_value <- BiCopPar2Tau(
    family = fit$family,
    par = fit$par,
    par2 = fit$par2
  )
  
  data.table(
    model       = model_name,
    family      = fit$family,
    family_name = family_name_safe(fit$family),
    par         = fit$par,
    par2        = fit$par2,
    kendall_tau = tau_value,
    logLik      = log_likelihood,
    AIC         = aic_value,
    BIC         = bic_value,
    n           = n_obs
  )
}

fit_candidate_models <- function(u_wind, u_price) {
  
  stopifnot(
    length(u_wind) == length(u_price),
    length(u_wind) > 20
  )
  
  u_wind  <- clip_uniform(u_wind)
  u_price <- clip_uniform(u_price)
  
  gaussian_fit <- BiCopEst(
    u1 = u_wind,
    u2 = u_price,
    family = 1,
    method = "mle"
  )
  
  student_fit <- BiCopEst(
    u1 = u_wind,
    u2 = u_price,
    family = 2,
    method = "mle"
  )
  
  bb8_fit <- BiCopSelect(
    u1 = u_wind,
    u2 = u_price,
    familyset = 10,
    rotations = TRUE,
    selectioncrit = "AIC",
    indeptest = FALSE
  )
  
  unrestricted_fit <- BiCopSelect(
    u1 = u_wind,
    u2 = u_price,
    familyset = base_family_set,
    rotations = TRUE,
    selectioncrit = "AIC",
    indeptest = FALSE
  )
  
  fit_list <- list(
    Gaussian = gaussian_fit,
    Student_t = student_fit,
    BB8 = bb8_fit,
    Unrestricted = unrestricted_fit
  )
  
  fit_table <- rbindlist(
    list(
      extract_fit_statistics(
        gaussian_fit, "Gaussian", length(u_wind)
      ),
      extract_fit_statistics(
        student_fit, "Student-t", length(u_wind)
      ),
      extract_fit_statistics(
        bb8_fit, "BB8", length(u_wind)
      ),
      extract_fit_statistics(
        unrestricted_fit, "AIC-selected", length(u_wind)
      )
    ),
    fill = TRUE
  )
  
  list(
    fits = fit_list,
    table = fit_table
  )
}

copula_tail_probabilities <- function(fit, q) {
  
  # Cannibalisation:
  # P(U_wind > 1-q, U_price < q)
  cannibalisation <- q - BiCopCDF(
    u1 = 1 - q,
    u2 = q,
    family = fit$family,
    par = fit$par,
    par2 = fit$par2
  )
  
  # Scarcity:
  # P(U_wind < q, U_price > 1-q)
  scarcity <- q - BiCopCDF(
    u1 = q,
    u2 = 1 - q,
    family = fit$family,
    par = fit$par,
    par2 = fit$par2
  )
  
  data.table(
    tail_level = q,
    cannibalisation_probability = max(0, cannibalisation),
    scarcity_probability = max(0, scarcity)
  )
}

tail_table_from_fits <- function(fit_list, regime_label) {
  
  rbindlist(
    lapply(names(fit_list), function(model_name) {
      
      fit <- fit_list[[model_name]]
      
      result <- rbindlist(
        lapply(
          tail_levels,
          function(q) copula_tail_probabilities(fit, q)
        )
      )
      
      result[, `:=`(
        regime = regime_label,
        model = model_name,
        family = fit$family,
        family_name = family_name_safe(fit$family)
      )]
      
      result
    }),
    fill = TRUE
  )
}

# ------------------------------------------------------------
# 4. Load, merge and clean data
# ------------------------------------------------------------

gen  <- fread(price_generation_file)
wind <- fread(wind_file)

gen[, date_time := mdy_hm(date_time, tz = "Australia/Melbourne")]
wind[, date_time := dmy_hm(date_time, tz = "Australia/Melbourne")]

dat <- merge(
  gen,
  wind,
  by = "date_time",
  all = FALSE
)

setorder(dat, date_time)

required_variables <- c(
  "date_time",
  "Regions_VIC_Price",
  "PORTWF",
  "Loc1"
)

missing_variables <- setdiff(required_variables, names(dat))

if (length(missing_variables) > 0) {
  stop(
    "Missing required variables: ",
    paste(missing_variables, collapse = ", ")
  )
}

dat <- dat[
  is.finite(Regions_VIC_Price) &
    is.finite(PORTWF) &
    is.finite(Loc1)
]

setnames(
  dat,
  old = c("Regions_VIC_Price", "PORTWF", "Loc1"),
  new = c("price", "generation", "wind")
)

dat <- unique(dat, by = "date_time")

stopifnot(
  nrow(dat) > 0,
  !anyDuplicated(dat$date_time)
)

cat("\nSample information\n")
cat("------------------\n")
cat("Observations:", nrow(dat), "\n")
cat("Start:", format(min(dat$date_time)), "\n")
cat("End:", format(max(dat$date_time)), "\n")

# ------------------------------------------------------------
# 5. Diagnostic correlations
# ------------------------------------------------------------

diagnostic_correlations <- data.table(
  comparison = c(
    "Price-generation",
    "Generation-wind",
    "Price-wind"
  ),
  pearson = c(
    cor(dat$price, dat$generation, method = "pearson"),
    cor(dat$generation, dat$wind, method = "pearson"),
    cor(dat$price, dat$wind, method = "pearson")
  ),
  kendall = c(
    cor(dat$price, dat$generation, method = "kendall"),
    cor(dat$generation, dat$wind, method = "kendall"),
    cor(dat$price, dat$wind, method = "kendall")
  )
)

print(diagnostic_correlations)

# ------------------------------------------------------------
# 6. Full-sample empirical pseudo-observations
# ------------------------------------------------------------

dat[, u_wind := make_pseudo_observations(wind)]
dat[, u_price := make_pseudo_observations(price)]

stopifnot(
  all(dat$u_wind > 0 & dat$u_wind < 1),
  all(dat$u_price > 0 & dat$u_price < 1)
)

# ------------------------------------------------------------
# 7. Full-sample copula estimates
# ------------------------------------------------------------

full_sample_models <- fit_candidate_models(
  dat$u_wind,
  dat$u_price
)

full_sample_fit_table <- copy(full_sample_models$table)

full_sample_fit_table[, `:=`(
  regime = "Full sample",
  start_date = min(dat$date_time),
  end_date = max(dat$date_time)
)]

setcolorder(
  full_sample_fit_table,
  c(
    "regime", "start_date", "end_date", "model",
    "family", "family_name", "par", "par2",
    "kendall_tau", "logLik", "AIC", "BIC", "n"
  )
)

print(full_sample_fit_table)

full_sample_tail_table <- tail_table_from_fits(
  full_sample_models$fits,
  regime_label = "Full sample"
)

# Relative probability compared with Gaussian
full_sample_tail_table[
  ,
  gaussian_cannibalisation :=
    cannibalisation_probability[model == "Gaussian"],
  by = tail_level
]

full_sample_tail_table[
  ,
  relative_to_gaussian :=
    cannibalisation_probability / gaussian_cannibalisation
]

print(full_sample_tail_table)

# ------------------------------------------------------------
# 8. Monthly endpoints for calendar-based rolling estimation
# ------------------------------------------------------------
dat <- dat[!is.na(dat$date_time), ]

first_possible_end <- min(dat$date_time) %m+% years(rolling_years)
last_possible_end  <- max(dat$date_time)

rolling_end_dates <- seq(
  from = floor_date(first_possible_end, unit = "month"),
  to = floor_date(last_possible_end, unit = "month"),
  by = paste(update_months, "months")
)

rolling_end_dates <- rolling_end_dates[
  rolling_end_dates >= first_possible_end &
    rolling_end_dates <= last_possible_end
]

if (length(rolling_end_dates) < 10) {
  stop("Insufficient monthly endpoints for rolling estimation.")
}

# ------------------------------------------------------------
# 9. Two-year rolling copula estimation
# ------------------------------------------------------------

rolling_results <- vector(
  mode = "list",
  length = length(rolling_end_dates)
)

for (j in seq_along(rolling_end_dates)) {
  
  window_end <- rolling_end_dates[j]
  window_start <- window_end %m-% years(rolling_years)
  
  window_data <- dat[
    date_time > window_start &
      date_time <= window_end
  ]
  
  if (nrow(window_data) < 1000) {
    warning(
      "Skipping rolling window ending ",
      format(window_end),
      ": insufficient observations."
    )
    next
  }
  
  # Recompute ranks within each rolling window.
  u_wind_window <- make_pseudo_observations(window_data$wind)
  u_price_window <- make_pseudo_observations(window_data$price)
  
  unrestricted_fit <- tryCatch(
    BiCopSelect(
      u1 = u_wind_window,
      u2 = u_price_window,
      familyset = base_family_set,
      rotations = TRUE,
      selectioncrit = "AIC",
      indeptest = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(unrestricted_fit)) {
    warning(
      "Copula estimation failed for window ending ",
      format(window_end)
    )
    next
  }
  
  gaussian_fit <- BiCopEst(
    u1 = u_wind_window,
    u2 = u_price_window,
    family = 1,
    method = "mle"
  )
  
  student_fit <- BiCopEst(
    u1 = u_wind_window,
    u2 = u_price_window,
    family = 2,
    method = "mle"
  )
  
  unrestricted_stats <- extract_fit_statistics(
    unrestricted_fit,
    "AIC-selected",
    nrow(window_data)
  )
  
  gaussian_stats <- extract_fit_statistics(
    gaussian_fit,
    "Gaussian",
    nrow(window_data)
  )
  
  student_stats <- extract_fit_statistics(
    student_fit,
    "Student-t",
    nrow(window_data)
  )
  
  rolling_results[[j]] <- data.table(
    window_start = window_start,
    window_end = window_end,
    n = nrow(window_data),
    
    selected_family = unrestricted_stats$family,
    selected_family_name = unrestricted_stats$family_name,
    selected_par = unrestricted_stats$par,
    selected_par2 = unrestricted_stats$par2,
    selected_tau = unrestricted_stats$kendall_tau,
    selected_logLik = unrestricted_stats$logLik,
    selected_AIC = unrestricted_stats$AIC,
    selected_BIC = unrestricted_stats$BIC,
    
    gaussian_tau = gaussian_stats$kendall_tau,
    gaussian_logLik = gaussian_stats$logLik,
    gaussian_AIC = gaussian_stats$AIC,
    gaussian_BIC = gaussian_stats$BIC,
    
    student_tau = student_stats$kendall_tau,
    student_logLik = student_stats$logLik,
    student_AIC = student_stats$AIC,
    student_BIC = student_stats$BIC
  )
}

rolling_results <- Filter(Negate(is.null), rolling_results)
rolling_table <- rbindlist(rolling_results, fill = TRUE)

if (nrow(rolling_table) < 10) {
  stop("Too few successful rolling estimates for break analysis.")
}

rolling_table[, rolling_index := .I]

print(rolling_table)

# ------------------------------------------------------------
# 10. Bai-Perron structural-break analysis
# ------------------------------------------------------------
library(strucchange)

break_model_full <- breakpoints(
  selected_tau ~ 1,
  data = rolling_table,
  h = minimum_regime_fraction
)

# BIC-selected number of breaks
break_model_selected <- breakpoints(break_model_full)

cat("\nBai-Perron break analysis\n")
cat("-------------------------\n")
print(summary(break_model_full))
print(break_model_selected)

break_indices <- break_model_selected$breakpoints
break_indices <- break_indices[!is.na(break_indices)]

if (length(break_indices) == 0) {
  
  warning(
    "BIC selected no structural breaks. ",
    "The subsequent analysis will contain one regime."
  )
  
  break_dates <- as.POSIXct(
    character(0),
    tz = "Australia/Melbourne"
  )
  
  rolling_table[, regime := 1L]
  
} else {
  
  break_dates <- rolling_table$window_end[break_indices]
  
  rolling_table[, regime := findInterval(
    rolling_index,
    vec = break_indices
  ) + 1L]
}

cat("\nEstimated break dates\n")
cat("---------------------\n")
print(break_dates)

# Confidence intervals for break locations
break_confidence_intervals <- tryCatch(
  confint(break_model_selected),
  error = function(e) NULL
)

if (!is.null(break_confidence_intervals)) {
  print(break_confidence_intervals)
}

# ------------------------------------------------------------
# 11. Assign regimes to original observations
# ------------------------------------------------------------

dat[, regime := 1L]

if (length(break_dates) > 0) {
  for (b in seq_along(break_dates)) {
    dat[date_time > break_dates[b], regime := b + 1L]
  }
}

regime_date_table <- dat[
  ,
  .(
    start_date = min(date_time),
    end_date = max(date_time),
    observations = .N
  ),
  by = regime
]

print(regime_date_table)

# ------------------------------------------------------------
# 12. Re-estimate competing copulas within each regime
# ------------------------------------------------------------

regime_model_objects <- list()
regime_fit_tables <- list()
regime_tail_tables <- list()

for (r in sort(unique(dat$regime))) {
  
  regime_data <- dat[regime == r]
  
  if (nrow(regime_data) < 100) {
    warning("Regime ", r, " has insufficient observations.")
    next
  }
  
  # Regime-specific empirical pseudo-observations
  u_wind_regime <- make_pseudo_observations(regime_data$wind)
  u_price_regime <- make_pseudo_observations(regime_data$price)
  
  regime_models <- fit_candidate_models(
    u_wind_regime,
    u_price_regime
  )
  
  regime_label <- paste("Regime", r)
  
  regime_fit_table <- copy(regime_models$table)
  
  regime_fit_table[, `:=`(
    regime = regime_label,
    regime_number = r,
    start_date = min(regime_data$date_time),
    end_date = max(regime_data$date_time)
  )]
  
  regime_tail_table <- tail_table_from_fits(
    regime_models$fits,
    regime_label = regime_label
  )
  
  regime_tail_table[, regime_number := r]
  
  regime_tail_table[
    ,
    gaussian_cannibalisation :=
      cannibalisation_probability[model == "Gaussian"],
    by = tail_level
  ]
  
  regime_tail_table[
    ,
    relative_to_gaussian :=
      cannibalisation_probability / gaussian_cannibalisation
  ]
  
  regime_model_objects[[regime_label]] <- regime_models$fits
  regime_fit_tables[[regime_label]] <- regime_fit_table
  regime_tail_tables[[regime_label]] <- regime_tail_table
}

regime_fit_table <- rbindlist(
  regime_fit_tables,
  fill = TRUE
)

regime_tail_table <- rbindlist(
  regime_tail_tables,
  fill = TRUE
)

setorder(
  regime_fit_table,
  regime_number,
  AIC
)

setorder(
  regime_tail_table,
  regime_number,
  tail_level,
  model
)

print(regime_fit_table)
print(regime_tail_table)

# ------------------------------------------------------------
# 13. Observed empirical tail frequencies by regime
# ------------------------------------------------------------

# Use fixed full-sample empirical thresholds for every regime.
wind_thresholds <- quantile(
  dat$wind,
  probs = 1 - tail_levels,
  na.rm = TRUE,
  names = FALSE
)

lower_price_thresholds <- quantile(
  dat$price,
  probs = tail_levels,
  na.rm = TRUE,
  names = FALSE
)

upper_price_thresholds <- quantile(
  dat$price,
  probs = 1 - tail_levels,
  na.rm = TRUE,
  names = FALSE
)

low_wind_thresholds <- quantile(
  dat$wind,
  probs = tail_levels,
  na.rm = TRUE,
  names = FALSE
)

observed_tail_table <- rbindlist(
  lapply(sort(unique(dat$regime)), function(r) {
    
    regime_data <- dat[regime == r]
    
    rbindlist(
      lapply(seq_along(tail_levels), function(k) {
        
        q <- tail_levels[k]
        
        data.table(
          regime_number = r,
          regime = paste("Regime", r),
          tail_level = q,
          
          high_wind_threshold = wind_thresholds[k],
          low_price_threshold = lower_price_thresholds[k],
          low_wind_threshold = low_wind_thresholds[k],
          high_price_threshold = upper_price_thresholds[k],
          
          observed_cannibalisation = mean(
            regime_data$wind > wind_thresholds[k] &
              regime_data$price < lower_price_thresholds[k]
          ),
          
          observed_scarcity = mean(
            regime_data$wind < low_wind_thresholds[k] &
              regime_data$price > upper_price_thresholds[k]
          ),
          
          observations = nrow(regime_data)
        )
      })
    )
  })
)

print(observed_tail_table)

# ------------------------------------------------------------
# 14. Model-implied versus observed tail probabilities
# ------------------------------------------------------------

comparison_tail_table <- merge(
  regime_tail_table,
  observed_tail_table[
    ,
    .(
      regime_number,
      tail_level,
      observed_cannibalisation,
      observed_scarcity
    )
  ],
  by = c("regime_number", "tail_level"),
  all.x = TRUE
)

comparison_tail_table[
  ,
  cannibalisation_error :=
    cannibalisation_probability -
    observed_cannibalisation
]

comparison_tail_table[
  ,
  scarcity_error :=
    scarcity_probability -
    observed_scarcity
]

print(comparison_tail_table)

# ------------------------------------------------------------
# 15. Rolling dependence figure
# ------------------------------------------------------------

break_plot_data <- data.table(
  break_date = break_dates
)

p_rolling <- ggplot(
  rolling_table,
  aes(x = window_end, y = selected_tau)
) +
  geom_line(
    colour = "firebrick",
    linewidth = 0.8
  ) +
  geom_point(
    aes(colour = selected_family_name),
    size = 1.6,
    alpha = 0.8
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_vline(
    data = break_plot_data,
    aes(xintercept = break_date),
    linetype = "dashed",
    colour = "black",
    linewidth = 0.55
  ) +
  labs(
    x = NULL,
    y = "Kendall's tau",
    colour = "AIC-selected family",
    title = "Rolling wind-price dependence and structural breaks",
    subtitle = paste0(
      rolling_years,
      "-year calendar windows updated monthly; ",
      "vertical lines show BIC-selected Bai-Perron breaks"
    )
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

print(p_rolling)

# ------------------------------------------------------------
# 16. Information-criterion comparison by regime
# ------------------------------------------------------------

p_aic <- ggplot(
  regime_fit_table,
  aes(
    x = model,
    y = AIC,
    fill = model
  )
) +
  geom_col(width = 0.7) +
  facet_wrap(
    ~ regime,
    scales = "free_y"
  ) +
  labs(
    x = NULL,
    y = "AIC",
    title = "Copula model comparison by dependence regime"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    plot.title = element_text(face = "bold")
  )

print(p_aic)

# ------------------------------------------------------------
# 17. Regime-specific cannibalisation probabilities
# ------------------------------------------------------------
library(scales)

p_tail <- ggplot(
  regime_tail_table,
  aes(
    x = factor(tail_level),
    y = cannibalisation_probability,
    colour = model,
    group = model
  )
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_wrap(~ regime) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.001)
  ) +
  labs(
    x = "Lower-price tail probability",
    y = "P(high wind, low price)",
    colour = "Copula model",
    title = "Regime-specific cannibalisation probabilities",
    subtitle = "Thresholds are fixed at corresponding empirical quantiles"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

print(p_tail)

# ------------------------------------------------------------
# 18. Selected-family frequencies across rolling windows
# ------------------------------------------------------------

rolling_family_frequency <- rolling_table[
  ,
  .N,
  by = .(
    selected_family,
    selected_family_name
  )
][
  order(-N)
]

rolling_family_frequency[
  ,
  share := N / sum(N)
]

print(rolling_family_frequency)

# ------------------------------------------------------------
# 19. Sensitivity to rolling-window length
# ------------------------------------------------------------

estimate_rolling_tau <- function(
    data,
    years_in_window,
    update_months = 1) {
  
  first_end <- min(data$date_time) %m+% years(years_in_window)
  
  end_dates <- seq(
    from = floor_date(first_end, "month"),
    to = floor_date(max(data$date_time), "month"),
    by = paste(update_months, "months")
  )
  
  end_dates <- end_dates[
    end_dates >= first_end &
      end_dates <= max(data$date_time)
  ]
  
  result <- vector("list", length(end_dates))
  
  for (i in seq_along(end_dates)) {
    
    end_date <- end_dates[i]
    start_date <- end_date %m-% years(years_in_window)
    
    d <- data[
      date_time > start_date &
        date_time <= end_date
    ]
    
    if (nrow(d) < 1000) {
      next
    }
    
    u1 <- make_pseudo_observations(d$wind)
    u2 <- make_pseudo_observations(d$price)
    
    fit <- tryCatch(
      BiCopSelect(
        u1,
        u2,
        familyset = base_family_set,
        rotations = TRUE,
        selectioncrit = "AIC",
        indeptest = FALSE
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit)) {
      next
    }
    
    result[[i]] <- data.table(
      window_years = years_in_window,
      window_start = start_date,
      window_end = end_date,
      n = nrow(d),
      family = fit$family,
      family_name = family_name_safe(fit$family),
      tau = BiCopPar2Tau(
        fit$family,
        fit$par,
        fit$par2
      )
    )
  }
  
  rbindlist(
    Filter(Negate(is.null), result),
    fill = TRUE
  )
}

window_sensitivity <- rbindlist(
  lapply(
    c(1, 2, 3),
    function(y) {
      estimate_rolling_tau(
        data = dat,
        years_in_window = y,
        update_months = update_months
      )
    }
  ),
  fill = TRUE
)

p_sensitivity <- ggplot(
  window_sensitivity,
  aes(
    x = window_end,
    y = tau,
    colour = factor(window_years)
  )
) +
  geom_line(linewidth = 0.75) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey40"
  ) +
  labs(
    x = NULL,
    y = "Kendall's tau",
    colour = "Window length\n(years)",
    title = "Sensitivity of rolling dependence to window length"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

print(p_sensitivity)

# ------------------------------------------------------------
# 20. Save reproducible outputs
# ------------------------------------------------------------

output_directory <- file.path(
  data_directory,
  "corrected_copula_results"
)

if (!dir.exists(output_directory)) {
  dir.create(
    output_directory,
    recursive = TRUE
  )
}

fwrite(
  diagnostic_correlations,
  file.path(output_directory, "diagnostic_correlations.csv")
)

fwrite(
  full_sample_fit_table,
  file.path(output_directory, "full_sample_copula_fits.csv")
)

fwrite(
  full_sample_tail_table,
  file.path(output_directory, "full_sample_tail_probabilities.csv")
)

fwrite(
  rolling_table,
  file.path(output_directory, "rolling_copula_estimates.csv")
)

fwrite(
  regime_date_table,
  file.path(output_directory, "regime_dates.csv")
)

fwrite(
  regime_fit_table,
  file.path(output_directory, "regime_copula_fits.csv")
)

fwrite(
  regime_tail_table,
  file.path(output_directory, "regime_tail_probabilities.csv")
)

fwrite(
  observed_tail_table,
  file.path(output_directory, "observed_regime_tail_frequencies.csv")
)

fwrite(
  comparison_tail_table,
  file.path(output_directory, "model_observed_tail_comparison.csv")
)

fwrite(
  rolling_family_frequency,
  file.path(output_directory, "rolling_family_frequencies.csv")
)

fwrite(
  window_sensitivity,
  file.path(output_directory, "rolling_window_sensitivity.csv")
)

if (length(break_dates) > 0) {
  fwrite(
    data.table(
      break_number = seq_along(break_dates),
      break_date = break_dates,
      rolling_index = break_indices
    ),
    file.path(output_directory, "bai_perron_break_dates.csv")
  )
}

ggsave(
  filename = file.path(
    output_directory,
    "rolling_dependence_breaks.png"
  ),
  plot = p_rolling,
  width = 9,
  height = 5.5,
  dpi = 400
)

ggsave(
  filename = file.path(
    output_directory,
    "regime_model_aic.png"
  ),
  plot = p_aic,
  width = 9,
  height = 5.5,
  dpi = 400
)

ggsave(
  filename = file.path(
    output_directory,
    "regime_tail_probabilities.png"
  ),
  plot = p_tail,
  width = 9,
  height = 5.5,
  dpi = 400
)

ggsave(
  filename = file.path(
    output_directory,
    "rolling_window_sensitivity.png"
  ),
  plot = p_sensitivity,
  width = 9,
  height = 5.5,
  dpi = 400
)

saveRDS(
  full_sample_models,
  file.path(
    output_directory,
    "full_sample_model_objects.rds"
  )
)

saveRDS(
  regime_model_objects,
  file.path(
    output_directory,
    "regime_model_objects.rds"
  )
)

saveRDS(
  break_model_full,
  file.path(
    output_directory,
    "bai_perron_full_model.rds"
  )
)

saveRDS(
  break_model_selected,
  file.path(
    output_directory,
    "bai_perron_selected_model.rds"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    output_directory,
    "session_info.txt"
  )
)


###  COMPARISON TAIL #####
### ADDENDUM ###

# ============================================================
# Correct observed probabilities using regime-specific quantiles
# ============================================================

observed_rank_tail_table <- rbindlist(
  lapply(sort(unique(dat$regime)), function(r) {
    
    regime_data <- copy(dat[regime == r])
    
    regime_data[, u_wind_regime :=
                  make_pseudo_observations(wind)]
    
    regime_data[, u_price_regime :=
                  make_pseudo_observations(price)]
    
    rbindlist(
      lapply(tail_levels, function(q) {
        
        data.table(
          regime_number = r,
          regime = paste("Regime", r),
          tail_level = q,
          
          observed_cannibalisation = mean(
            regime_data$u_wind_regime > 1 - q &
              regime_data$u_price_regime < q
          ),
          
          observed_scarcity = mean(
            regime_data$u_wind_regime < q &
              regime_data$u_price_regime > 1 - q
          ),
          
          observations = nrow(regime_data)
        )
      })
    )
  })
)

comparison_rank_tail_table <- merge(
  regime_tail_table,
  observed_rank_tail_table,
  by = c("regime_number", "regime", "tail_level"),
  all.x = TRUE
)

comparison_rank_tail_table[
  ,
  cannibalisation_error :=
    cannibalisation_probability -
    observed_cannibalisation
]

comparison_rank_tail_table[
  ,
  scarcity_error :=
    scarcity_probability -
    observed_scarcity
]

print(observed_rank_tail_table)
print(comparison_rank_tail_table)

sample_start_date <- min(dat$date_time, na.rm = TRUE)
sample_end_date   <- max(dat$date_time, na.rm = TRUE)

full_sample_fit_table[, `:=`(
  regime = "Full sample",
  start_date = sample_start_date,
  end_date = sample_end_date
)]

full_sample_fit_table

### ANCILLARY DETAILS
names(gen)
names(wind)

str(gen[, .(date_time, Regions_VIC_Price, PORTWF)])
str(wind[, .(date_time, Loc1)])

range(gen$date_time, na.rm = TRUE)
range(wind$date_time, na.rm = TRUE)

nrow(gen)
nrow(wind)
nrow(dat)

data.table(
  source_price_rows = nrow(gen),
  source_wind_rows  = nrow(wind),
  matched_rows      = nrow(merge(gen, wind, by = "date_time")),
  final_rows        = nrow(dat)
)

colSums(is.na(dat[, .(
  price,
  generation,
  wind
)]))

sum(duplicated(gen$date_time))
sum(duplicated(wind$date_time))
sum(duplicated(dat$date_time))

attr(gen$date_time, "tzone")
attr(wind$date_time, "tzone")
attr(dat$date_time, "tzone")

table(hour(dat$date_time))
summary(diff(dat$date_time))

break_confidence_intervals <- tryCatch(
  confint(break_model_selected),
  error = function(e) e
)

print(break_confidence_intervals)
class(break_confidence_intervals)

# ============================================================
# Moving-block bootstrap uncertainty for copula tail probabilities
# ============================================================

library(data.table)
library(VineCopula)

set.seed(20261007)

# ------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------

B <- 99                 # use 199 for an initial test run
block_length <- 168      # one week of hourly observations
tail_levels <- c(0.10, 0.05, 0.01)

family_set <- c(
  1,  # Gaussian
  2,  # Student-t
  3,  # Clayton
  4,  # Gumbel
  5,  # Frank
  6,  # Joe
  7,  # BB1
  8,  # BB6
  9,  # BB7
  10  # BB8
)

output_directory <- file.path(
  data_directory,
  "corrected_copula_results"
)

if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

# ------------------------------------------------------------
# 2. Check required data
# ------------------------------------------------------------

required_columns <- c(
  "date_time",
  "wind",
  "price",
  "regime"
)

missing_columns <- setdiff(required_columns, names(dat))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

setorder(dat, date_time)

stopifnot(
  all(is.finite(dat$wind)),
  all(is.finite(dat$price)),
  !anyDuplicated(dat$date_time)
)

cat("Final observations:", nrow(dat), "\n")
cat("Sample start:", format(min(dat$date_time)), "\n")
cat("Sample end:", format(max(dat$date_time)), "\n")

# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

clip_uniform <- function(x, epsilon = 1e-6) {
  pmin(pmax(x, epsilon), 1 - epsilon)
}

make_pseudo_observations <- function(x) {
  clip_uniform(
    rank(x, ties.method = "average") /
      (length(x) + 1)
  )
}

family_name_safe <- function(family) {
  tryCatch(
    BiCopName(family),
    error = function(e) paste0("Family_", family)
  )
}

# Moving-block bootstrap indices.
# Blocks are sampled jointly for wind and price.
moving_block_indices <- function(n, block_length) {
  
  if (block_length >= n) {
    stop("Block length must be smaller than sample size.")
  }
  
  number_of_blocks <- ceiling(n / block_length)
  
  possible_starts <- seq_len(
    n - block_length + 1
  )
  
  sampled_starts <- sample(
    possible_starts,
    size = number_of_blocks,
    replace = TRUE
  )
  
  indices <- unlist(
    lapply(
      sampled_starts,
      function(s) {
        s:(s + block_length - 1)
      }
    ),
    use.names = FALSE
  )
  
  indices[seq_len(n)]
}

fit_bootstrap_models <- function(wind, price) {
  
  u_wind <- make_pseudo_observations(wind)
  u_price <- make_pseudo_observations(price)
  
  gaussian_fit <- BiCopEst(
    u1 = u_wind,
    u2 = u_price,
    family = 1,
    method = "mle"
  )
  
  student_fit <- BiCopEst(
    u1 = u_wind,
    u2 = u_price,
    family = 2,
    method = "mle"
  )
  
  selected_fit <- BiCopSelect(
    u1 = u_wind,
    u2 = u_price,
    familyset = family_set,
    rotations = TRUE,
    selectioncrit = "AIC",
    indeptest = FALSE
  )
  
  list(
    Gaussian = gaussian_fit,
    Student_t = student_fit,
    Selected = selected_fit
  )
}

copula_tail_probabilities <- function(fit, q) {
  
  cannibalisation <- q - BiCopCDF(
    u1 = 1 - q,
    u2 = q,
    family = fit$family,
    par = fit$par,
    par2 = fit$par2
  )
  
  scarcity <- q - BiCopCDF(
    u1 = q,
    u2 = 1 - q,
    family = fit$family,
    par = fit$par,
    par2 = fit$par2
  )
  
  data.table(
    tail_level = q,
    cannibalisation_probability =
      max(0, min(q, cannibalisation)),
    scarcity_probability =
      max(0, min(q, scarcity))
  )
}

extract_bootstrap_probabilities <- function(
    fit_list,
    replication,
    sample_label,
    regime_number = NA_integer_) {
  
  result <- rbindlist(
    lapply(
      names(fit_list),
      function(model_name) {
        
        fit <- fit_list[[model_name]]
        
        probability_table <- rbindlist(
          lapply(
            tail_levels,
            function(q) {
              copula_tail_probabilities(fit, q)
            }
          )
        )
        
        probability_table[, `:=`(
          replication = replication,
          sample = sample_label,
          regime_number = regime_number,
          model = model_name,
          family = fit$family,
          family_name =
            family_name_safe(fit$family),
          par = fit$par,
          par2 = fit$par2
        )]
        
        probability_table
      }
    ),
    fill = TRUE
  )
  
  # Selected-minus-Gaussian differences
  wide <- dcast(
    result,
    replication +
      sample +
      regime_number +
      tail_level ~ model,
    value.var = "cannibalisation_probability"
  )
  
  if (all(c("Selected", "Gaussian") %in% names(wide))) {
    wide[, selected_minus_gaussian :=
           Selected - Gaussian]
    
    wide[, selected_to_gaussian :=
           Selected / Gaussian]
  }
  
  list(
    probabilities = result,
    differences = wide
  )
}

summarise_bootstrap <- function(
    x,
    grouping_variables,
    value_variable) {
  
  x[
    is.finite(get(value_variable)),
    .(
      bootstrap_mean =
        mean(get(value_variable)),
      bootstrap_sd =
        sd(get(value_variable)),
      lower_95 =
        quantile(
          get(value_variable),
          probs = 0.025,
          na.rm = TRUE,
          names = FALSE
        ),
      median =
        quantile(
          get(value_variable),
          probs = 0.500,
          na.rm = TRUE,
          names = FALSE
        ),
      upper_95 =
        quantile(
          get(value_variable),
          probs = 0.975,
          na.rm = TRUE,
          names = FALSE
        ),
      successful_replications = .N
    ),
    by = grouping_variables
  ]
}

# ------------------------------------------------------------
# 4. Original-sample estimates
# ------------------------------------------------------------

original_full_fits <- fit_bootstrap_models(
  wind = dat$wind,
  price = dat$price
)

original_full <- extract_bootstrap_probabilities(
  fit_list = original_full_fits,
  replication = 0L,
  sample_label = "Full sample"
)

original_full_probabilities <-
  original_full$probabilities

original_full_differences <-
  original_full$differences

# ------------------------------------------------------------
# 5. Full-sample moving-block bootstrap
# ------------------------------------------------------------

full_probability_results <- vector(
  mode = "list",
  length = B
)

full_difference_results <- vector(
  mode = "list",
  length = B
)

full_error_log <- vector(
  mode = "character",
  length = B
)

for (b in seq_len(B)) {
  
  if (b %% 25 == 0) {
    cat(
      "Full-sample bootstrap:",
      b,
      "of",
      B,
      "\n"
    )
  }
  
  bootstrap_result <- tryCatch({
    
    index <- moving_block_indices(
      n = nrow(dat),
      block_length = block_length
    )
    
    bootstrap_data <- dat[index]
    
    fit_list <- fit_bootstrap_models(
      wind = bootstrap_data$wind,
      price = bootstrap_data$price
    )
    
    extract_bootstrap_probabilities(
      fit_list = fit_list,
      replication = b,
      sample_label = "Full sample"
    )
    
  }, error = function(e) {
    full_error_log[b] <<- conditionMessage(e)
    NULL
  })
  
  if (!is.null(bootstrap_result)) {
    full_probability_results[[b]] <-
      bootstrap_result$probabilities
    
    full_difference_results[[b]] <-
      bootstrap_result$differences
  }
}

full_bootstrap_probabilities <- rbindlist(
  Filter(
    Negate(is.null),
    full_probability_results
  ),
  fill = TRUE
)

full_bootstrap_differences <- rbindlist(
  Filter(
    Negate(is.null),
    full_difference_results
  ),
  fill = TRUE
)

# ------------------------------------------------------------
# 6. Full-sample uncertainty summaries
# ------------------------------------------------------------

full_probability_summary <- summarise_bootstrap(
  x = full_bootstrap_probabilities,
  grouping_variables = c(
    "model",
    "tail_level"
  ),
  value_variable =
    "cannibalisation_probability"
)

full_scarcity_summary <- summarise_bootstrap(
  x = full_bootstrap_probabilities,
  grouping_variables = c(
    "model",
    "tail_level"
  ),
  value_variable =
    "scarcity_probability"
)

full_difference_summary <- summarise_bootstrap(
  x = full_bootstrap_differences,
  grouping_variables = c(
    "tail_level"
  ),
  value_variable =
    "selected_minus_gaussian"
)

full_ratio_summary <- summarise_bootstrap(
  x = full_bootstrap_differences,
  grouping_variables = c(
    "tail_level"
  ),
  value_variable =
    "selected_to_gaussian"
)

full_selected_family_frequency <-
  unique(
    full_bootstrap_probabilities[
      model == "Selected",
      .(
        replication,
        family,
        family_name
      )
    ]
  )[
    ,
    .N,
    by = .(
      family,
      family_name
    )
  ][
    order(-N)
  ]

full_selected_family_frequency[
  ,
  share := N / sum(N)
]

# Add original estimates to summaries
original_full_cannibalisation <-
  original_full_probabilities[
    ,
    .(
      model,
      tail_level,
      original_estimate =
        cannibalisation_probability
    )
  ]

original_full_scarcity <-
  original_full_probabilities[
    ,
    .(
      model,
      tail_level,
      original_estimate =
        scarcity_probability
    )
  ]

original_full_difference <-
  original_full_differences[
    ,
    .(
      tail_level,
      original_difference =
        selected_minus_gaussian,
      original_ratio =
        selected_to_gaussian
    )
  ]

full_probability_summary <- merge(
  original_full_cannibalisation,
  full_probability_summary,
  by = c("model", "tail_level"),
  all.x = TRUE
)

full_scarcity_summary <- merge(
  original_full_scarcity,
  full_scarcity_summary,
  by = c("model", "tail_level"),
  all.x = TRUE
)

full_difference_summary <- merge(
  original_full_difference[
    ,
    .(
      tail_level,
      original_difference
    )
  ],
  full_difference_summary,
  by = "tail_level",
  all.x = TRUE
)

full_ratio_summary <- merge(
  original_full_difference[
    ,
    .(
      tail_level,
      original_ratio
    )
  ],
  full_ratio_summary,
  by = "tail_level",
  all.x = TRUE
)

# ------------------------------------------------------------
# 7. Regime-level q = 0.10 bootstrap
# ------------------------------------------------------------

regime_probability_results <- list()
regime_difference_results <- list()
regime_error_log <- list()

regime_numbers <- sort(unique(dat$regime))

for (r in regime_numbers) {
  
  cat("\nStarting regime", r, "\n")
  
  regime_key <- as.character(r)
  
  regime_data <- copy(
    dat[regime == r]
  )
  
  regime_probability_results[[regime_key]] <-
    vector("list", B)
  
  regime_difference_results[[regime_key]] <-
    vector("list", B)
  
  regime_error_log[[regime_key]] <-
    character(B)
  
  for (b in seq_len(B)) {
    
    if (b %% 25 == 0) {
      cat(
        "Regime",
        r,
        ":",
        b,
        "of",
        B,
        "\n"
      )
    }
    
    bootstrap_result <- tryCatch({
      
      index <- moving_block_indices(
        n = nrow(regime_data),
        block_length = block_length
      )
      
      bootstrap_data <- regime_data[index]
      
      fit_list <- fit_bootstrap_models(
        wind = bootstrap_data$wind,
        price = bootstrap_data$price
      )
      
      extracted <- extract_bootstrap_probabilities(
        fit_list = fit_list,
        replication = b,
        sample_label = paste("Regime", r),
        regime_number = r
      )
      
      # Retain q = 0.10 only for regime inference
      extracted$probabilities <-
        extracted$probabilities[
          tail_level == 0.10
        ]
      
      extracted$differences <-
        extracted$differences[
          tail_level == 0.10
        ]
      
      extracted
      
    }, error = function(e) {
      
      regime_error_log[[regime_key]][b] <<-
        conditionMessage(e)
      
      NULL
    })
    
    if (!is.null(bootstrap_result)) {
      
      regime_probability_results[[regime_key]][[b]] <- bootstrap_result$probabilities
      
      regime_difference_results[[regime_key]][[b]] <- bootstrap_result$differences
    }
  }
}

regime_bootstrap_probabilities <- rbindlist(
  lapply(
    regime_probability_results,
    function(x) {
      rbindlist(
        Filter(Negate(is.null), x),
        fill = TRUE
      )
    }
  ),
  fill = TRUE
)

regime_bootstrap_differences <- rbindlist(
  lapply(
    regime_difference_results,
    function(x) {
      rbindlist(
        Filter(Negate(is.null), x),
        fill = TRUE
      )
    }
  ),
  fill = TRUE
)

# ------------------------------------------------------------
# 8. Regime-level uncertainty summaries
# ------------------------------------------------------------

regime_difference_summary <- summarise_bootstrap(
  x = regime_bootstrap_differences,
  grouping_variables = c(
    "regime_number",
    "sample",
    "tail_level"
  ),
  value_variable =
    "selected_minus_gaussian"
)

regime_ratio_summary <- summarise_bootstrap(
  x = regime_bootstrap_differences,
  grouping_variables = c(
    "regime_number",
    "sample",
    "tail_level"
  ),
  value_variable =
    "selected_to_gaussian"
)

regime_selected_family_frequency <-
  unique(
    regime_bootstrap_probabilities[
      model == "Selected",
      .(
        regime_number,
        replication,
        family,
        family_name
      )
    ]
  )[
    ,
    .N,
    by = .(
      regime_number,
      family,
      family_name
    )
  ][
    ,
    share := N / sum(N),
    by = regime_number
  ][
    order(regime_number, -N)
  ]

# ------------------------------------------------------------
# 9. Original regime q = 0.10 estimates
# ------------------------------------------------------------

original_regime_results <- rbindlist(
  lapply(
    regime_numbers,
    function(r) {
      
      regime_data <- dat[regime == r]
      
      fits <- fit_bootstrap_models(
        wind = regime_data$wind,
        price = regime_data$price
      )
      
      extracted <- extract_bootstrap_probabilities(
        fit_list = fits,
        replication = 0L,
        sample_label = paste("Regime", r),
        regime_number = r
      )
      
      extracted$differences[
        tail_level == 0.10
      ]
    }
  ),
  fill = TRUE
)

regime_difference_summary <- merge(
  original_regime_results[
    ,
    .(
      regime_number,
      sample,
      tail_level,
      original_difference =
        selected_minus_gaussian
    )
  ],
  regime_difference_summary,
  by = c(
    "regime_number",
    "sample",
    "tail_level"
  ),
  all.x = TRUE
)

regime_ratio_summary <- merge(
  original_regime_results[
    ,
    .(
      regime_number,
      sample,
      tail_level,
      original_ratio =
        selected_to_gaussian
    )
  ],
  regime_ratio_summary,
  by = c(
    "regime_number",
    "sample",
    "tail_level"
  ),
  all.x = TRUE
)

# ------------------------------------------------------------
# 10. Error diagnostics
# ------------------------------------------------------------

full_errors <- data.table(
  replication = seq_len(B),
  error = full_error_log
)[
  nzchar(error)
]

regime_errors <- rbindlist(
  lapply(
    names(regime_error_log),
    function(r) {
      
      data.table(
        regime_number = as.integer(r),
        replication = seq_len(B),
        error = regime_error_log[[r]]
      )[
        nzchar(error)
      ]
    }
  ),
  fill = TRUE
)

cat("\nBootstrap diagnostics\n")
cat("---------------------\n")

cat(
  "Successful full-sample replications:",
  uniqueN(full_bootstrap_probabilities$replication),
  "of",
  B,
  "\n"
)

cat(
  "Failed full-sample replications:",
  nrow(full_errors),
  "\n"
)

print(
  regime_bootstrap_differences[
    ,
    .(
      successful_replications =
        uniqueN(replication)
    ),
    by = regime_number
  ]
)

# ------------------------------------------------------------
# 11. Save outputs
# ------------------------------------------------------------

fwrite(
  full_probability_summary,
  file.path(
    output_directory,
    "bootstrap_full_cannibalisation_summary.csv"
  )
)

fwrite(
  full_scarcity_summary,
  file.path(
    output_directory,
    "bootstrap_full_scarcity_summary.csv"
  )
)

fwrite(
  full_difference_summary,
  file.path(
    output_directory,
    "bootstrap_full_selected_minus_gaussian.csv"
  )
)

fwrite(
  full_ratio_summary,
  file.path(
    output_directory,
    "bootstrap_full_selected_to_gaussian.csv"
  )
)

fwrite(
  full_selected_family_frequency,
  file.path(
    output_directory,
    "bootstrap_full_family_frequencies.csv"
  )
)

fwrite(
  regime_difference_summary,
  file.path(
    output_directory,
    "bootstrap_regime_selected_minus_gaussian_q10.csv"
  )
)

fwrite(
  regime_ratio_summary,
  file.path(
    output_directory,
    "bootstrap_regime_selected_to_gaussian_q10.csv"
  )
)

fwrite(
  regime_selected_family_frequency,
  file.path(
    output_directory,
    "bootstrap_regime_family_frequencies_q10.csv"
  )
)

fwrite(
  full_errors,
  file.path(
    output_directory,
    "bootstrap_full_errors.csv"
  )
)

fwrite(
  regime_errors,
  file.path(
    output_directory,
    "bootstrap_regime_errors.csv"
  )
)

saveRDS(
  full_bootstrap_probabilities,
  file.path(
    output_directory,
    "bootstrap_full_raw_probabilities.rds"
  )
)

saveRDS(
  regime_bootstrap_differences,
  file.path(
    output_directory,
    "bootstrap_regime_raw_differences.rds"
  )
)

cat("\nBootstrap analysis complete.\n")
cat("Outputs saved to:\n")
cat(output_directory, "\n")


cat("\nEstimation complete.\n")
cat("Outputs saved to:\n", output_directory, "\n")