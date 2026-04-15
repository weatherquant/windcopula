# ============================================================
# Copula-Based Modelling of Wind–Price Risk
# Jason West
# Bureau of Meteorology, Australia
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

# ---- load data ----
setwd('Data/Wind/Copulas')
gen  <- fread("pricegen.csv")
wind <- fread("winds.csv")

gen[, date_time := mdy_hm(date_time)]
wind[, date_time := dmy_hm(date_time)]

dat <- merge(gen, wind, by = "date_time")

# sanity checks
stopifnot(nrow(dat) > 0)

# ---- select modelling variables ----
price <- dat$Regions_VIC_Price
genmw <- dat$PORTWF          # actual generation
windm <- dat$Loc1            # ERA5 / proxy wind speed

# ---- correlations (diagnostic only) ----
cor(price, genmw, use = "complete.obs")
cor(genmw, windm, use = "complete.obs")

# ---- joint plots ----
df_plot <- data.frame(price = price, wind = windm)

p <- ggplot(df_plot, aes(wind, price)) +
  geom_point(alpha = 0.35) +
  labs(x = "Wind speed", y = "Dispatch price") +
  labs(
    x = "Wind speed (m/s)",
    y = "Dispatch price (AUD/MWh)",
    title = "Observed Joint Behaviour of Wind Speed and Dispatch Prices"
  ) + ylim(c(-200,2000))

# ---- Figure 1 ----
ggMarginal(p, type = "density")

# ---- marginal fit gamma to wind ----
fit_wind  <- fitdist(windm, "gamma")

# ---- shift prices to positive support ----
price_shift <- price - min(price) + 1e-6

# ---- fit Weibull to shifted prices ----
fit_price <- fitdist(price_shift, "weibull")

# ------------------------------------------------------------
# Probability Integral Transforms
# ------------------------------------------------------------
u_wind <- pgamma(
  windm,
  shape = fit_wind$estimate["shape"],
  rate  = fit_wind$estimate["rate"]
)

u_price <- pweibull(
  price_shift,
  shape = fit_price$estimate["shape"],
  scale = fit_price$estimate["scale"]
)

U <- cbind(u_wind, u_price)

# ---- fit copula families ----
cop_fit <- BiCopSelect(
  U[,1], U[,2],
  familyset = c(1, 2, 6, 7, 10, 17),  # Gaussian, t, Clayton, Gumbel, Frank, BB8
  selectioncrit = "AIC",
  indeptest = TRUE
)

summary(cop_fit)

# ---- simulate joint distribution ----
set.seed(123)

U_sim <- BiCopSim(
  N = 100000,
  family = cop_fit$family,
  par = cop_fit$par,
  par2 = cop_fit$par2
)

# ---- inverse PIT ----
wind_sim <- qgamma(U_sim[,1],
                   shape = fit_wind$estimate["shape"],
                   rate  = fit_wind$estimate["rate"])

# ---- inverse PIT (price, shifted Weibull) ----
price_sim_shift <- qweibull(
  U_sim[,2],
  shape = fit_price$estimate["shape"],
  scale = fit_price$estimate["scale"]
)

# ---- shift prices back to original scale ----
price_sim <- price_sim_shift + min(price, na.rm = TRUE) - 1e-6

low_wind  <- quantile(wind_sim, 0.10)
high_wind <- quantile(wind_sim, 0.90)

low_price  <- quantile(price_sim, 0.10)
high_price <- quantile(price_sim, 0.90)

P_cannibal <- mean(wind_sim > high_wind & price_sim < low_price)
P_scarcity <- mean(wind_sim < low_wind  & price_sim > high_price)

# ---- PLOTS ----------
df_obs <- data.frame(wind = windm, price = price)

# Spare Figure Purpose: Show nonlinearity and tail clustering invisible to correlation.
ggplot(df_obs, aes(wind, price)) +
  geom_point(alpha = 0.42, size = 0.8) +
  geom_smooth(method = "loess", color = "red", se = FALSE) +
  labs(
    x = "Wind speed (m/s)",
    y = "Dispatch price (AUD/MWh)",
    title = "Observed Joint Behaviour of Wind Speed and Dispatch Prices"
  ) +
  theme_bw() + ylim(c(-200,2000))

# Fig 1 amended version
set.seed(1)
df_obs_samp <- df_obs[sample(nrow(df_obs), 6000), ]

p_obs <- ggplot(df_obs_samp, aes(wind, price)) +
  geom_point(alpha = 0.25, size = 0.6, colour = "grey30") +
  geom_smooth(
    method = "loess",
    colour = "firebrick",
    se = FALSE,
    linewidth = 0.9
  ) +
  labs(
    x = "Wind speed (m/s)",
    y = "Dispatch price (AUD/MWh)",
    title = "(a) Observed joint behaviour",
    subtitle = "Weak average relationship conceals nonlinear tail structure"
  ) +
  coord_cartesian(ylim = c(-200, 2000)) +
  theme_bw()

# Extract Kendall's tau from selected copula
tau_hat <- BiCopPar2Tau(
  cop_fit$family,
  cop_fit$par,
  cop_fit$par2
)

# Convert tau -> Gaussian correlation
rho_gauss <- sin(tau_hat * pi / 2)

# Simulate Gaussian copula
U_gauss <- BiCopSim(
  N = 50000,
  family = 1,
  par = rho_gauss
)

wind_gauss <- qgamma(
  U_gauss[,1],
  shape = fit_wind$estimate["shape"],
  rate  = fit_wind$estimate["rate"]
)

price_gauss_shift <- qweibull(
  U_gauss[,2],
  shape = fit_price$estimate["shape"],
  scale = fit_price$estimate["scale"]
)

price_gauss <- price_gauss_shift + min(price, na.rm = TRUE) - 1e-6

df_gauss <- data.frame(
  wind = wind_gauss,
  price = price_gauss,
  model = "Gaussian copula"
)

df_bb8 <- data.frame(
  wind = wind_sim,
  price = price_sim,
  model = "BB8 copula"
)

df_sim <- rbind(df_gauss, df_bb8)
df_sim <- df_sim[sample(nrow(df_sim), 12000), ]

# tail thresholds (compute once)
w90 <- quantile(windm, 0.9, na.rm = TRUE)
p10 <- quantile(price, 0.1, na.rm = TRUE)

p_sim <- ggplot(df_sim, aes(wind, price)) +
  geom_point(alpha = 0.18, size = 0.5) +
  geom_vline(xintercept = w90, linetype = "dashed",
             colour = "firebrick", linewidth = 0.4) +
  geom_hline(yintercept = p10, linetype = "dashed",
             colour = "firebrick", linewidth = 0.4) +
  facet_wrap(~model, ncol = 2) +
  labs(
    x = "Wind speed (m/s)",
    y = "Dispatch price (AUD/MWh)",
    title = "(b) Identical marginals, different dependence structures",
    subtitle = "Flexible copula captures asymmetric tail risk missed by Gaussian dependence"
  ) +
  coord_cartesian(ylim = c(-200, 2000)) +
  theme_bw()

library(patchwork)

# ---- Figure 2 ----
p_obs / p_sim +
  plot_annotation(
    title = "Wind-Price Dependence in a High Renewables Electricity Market",
    subtitle = "Correlation masks tail risk; flexible copulas recover asymmetric dependence",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13)
    )
  )

# Spare figure
P_gauss <- mean(
  wind_gauss > w90 & price_gauss < p10
)

P_bb8 <- mean(
  wind_sim > w90 & price_sim < p10
)

data.frame(
  Model = c("Gaussian", "BB8"),
  Probability = c(P_gauss, P_bb8)
)

# Spare figure - contours
ggplot(df_sim, aes(wind, price)) +
  stat_density_2d(
    aes(colour = after_stat(level)),
    bins = 6,
    linewidth = 0.7
  ) +
  facet_wrap(~model, ncol = 2) +
  labs(
    x = "Wind speed (m/s)",
    y = "Dispatch price (AUD/MWh)",
    title = "Joint density contours under alternative dependence structures",
    subtitle = "Differences emerge in the lower tail despite similar average dependence"
  ) +
  coord_cartesian(ylim = c(-200, 2000)) +
  theme_bw() +
  theme(legend.position = "none")

# Fig Purpose: Visually prove why Gaussian copulas fail.
# Gaussian copula
cop_gauss <- BiCopEst(U[,1], U[,2], family = 1)

# BB8 copula (already selected)
cop_bb8 <- cop_fit

grid <- seq(0.01, 0.99, length.out = 100)
dens_grid <- expand.grid(u = grid, v = grid)

dens_grid$gauss <- BiCopPDF(
  dens_grid$u, dens_grid$v,
  family = cop_gauss$family,
  par = cop_gauss$par
)

dens_grid$bb8 <- BiCopPDF(
  dens_grid$u, dens_grid$v,
  family = cop_bb8$family,
  par = cop_bb8$par,
  par2 = cop_bb8$par2
)

ggplot(dens_grid, aes(u, v, fill = gauss)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Gaussian Copula Density") +
  theme_bw()

ggplot(dens_grid, aes(u, v, fill = bb8)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Gaussian Copula Density") +
  theme_bw()

# ---- Figure 4 ----
# Show simulated joint extremes driving pricing
df_sim <- data.frame(wind = wind_sim, price = price_sim)

ggplot(df_sim, aes(wind, price)) +
  geom_point(alpha = 0.20, color = "steelblue") +
  geom_vline(xintercept = high_wind, linetype = "dashed") +
  geom_hline(yintercept = low_price, linetype = "dashed") +
  labs(
    x = "Simulated wind speed",
    y = "Simulated price",
    title = "Simulated Cannibalisation Region"
  ) +
  theme_bw() + annotate(
    "rect",
    xmin = high_wind, xmax = Inf,
    ymin = -Inf, ymax = low_price,
    fill = "red", alpha = 0.05
  )

# Tail Event Probabilities Under Each Model
simulate_tail_probs <- function(cop) {
  
  U_sim <- BiCopSim(
    N = 100000,
    family = cop$family,
    par = cop$par,
    par2 = cop$par2
  )
  
  wind_sim <- qgamma(
    U_sim[,1],
    shape = fit_wind$estimate["shape"],
    rate  = fit_wind$estimate["rate"]
  )
  
  price_sim <- qweibull(
    U_sim[,2],
    shape = fit_price$estimate["shape"],
    scale = fit_price$estimate["scale"]
  ) + min(price) - 1e-6
  
  list(
    cannibal = mean(wind_sim > high_wind & price_sim < low_price),
    scarcity = mean(wind_sim < low_wind & price_sim > high_price)
  )
}

P_gauss <- simulate_tail_probs(cop_gauss)
P_bb8   <- simulate_tail_probs(cop_bb8)

# Pricing Bias Table
pricing_bias <- data.frame(
  Model = c("Gaussian copula", "BB8 copula"),
  P_cannibalisation = c(P_gauss$cannibal, P_bb8$cannibal),
  P_scarcity = c(P_gauss$scarcity, P_bb8$scarcity)
)

pricing_bias$Cannibalisation_Bias =
  pricing_bias$P_cannibalisation /
  pricing_bias$P_cannibalisation[1]

pricing_bias

# Gaussian systematically underprices downside risk
# BB8 produces materially higher expected losses

# ---- Rolling Copula Estimation --------
window <- 365 * 2  # 2‑year rolling window
step   <- 30       # monthly update

dates <- dat$date_time
tau_roll <- c()
date_mid <- c()

for (i in seq(1, length(U[,1]) - window, by = step)) {
  
  U_sub <- U[i:(i+window), ]
  
  fit <- BiCopSelect(
    U_sub[,1], U_sub[,2],
    familyset = c(1,2,6,7,10,17),
    selectioncrit = "AIC"
  )
  
  tau_roll <- c(tau_roll, fit$tau)
  date_mid <- c(date_mid, dates[i + window/2])
}

df_tau <- data.frame(date = date_mid, tau = tau_roll)

ggplot(df_tau, aes(date, tau)) +
  geom_line(color = "darkred") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    y = "Kendall's Tau",
    title = "Time‑Varying Dependence Between Wind and Prices"
  ) +
  theme_bw()


# ---- Figure 3 ----
# ============================================================
# Cullen–Frey diagnostics: two‑panel figure
# ============================================================
library(fitdistrplus)

# ---- prepare shifted prices (required for Weibull support) ----
price_shift <- price - min(price, na.rm = TRUE) + 1e-6

# ---- set plotting layout ----
par(mfrow = c(1, 2),           # two panels, side‑by‑side
    mar = c(4, 4, 3, 1),       # margins
    oma = c(0, 0, 2, 0))       # outer margin for overall title
par(ann=FALSE) # Disables all titles and labels

# ---- Panel A: Wind speed ----
descdist(
  windm,
  discrete = FALSE,
  boot = 1000
)

title(main = "(a) Wind speed", cex.main = 1)

# ---- Panel B: Dispatch prices (shifted) ----
par(yaxt = "n")  # suppress default y-axis

descdist(
  price_shift,
  discrete = FALSE,
  boot = 1000
)

# add simple y-axis ticks (kurtosis)
axis(
  side = 2,
  at = c(2, 5, 10),
  labels = c("2", "5", "10")
)

title(main = "(b) Dispatch prices (shifted)", cex.main = 1)

# ---- overall figure title ----
mtext("Cullen–Frey Skewness–Kurtosis Diagnostics",
      outer = TRUE,
      cex = 1.1)

# ---- reset plotting defaults ----
par(mfrow = c(1, 1), yaxt = "s", ann=TRUE)

# ---- Figure 5 ----
# ============================================================
# Copula figures
# ============================================================
### Simulation function (reuse for both copulas)
simulate_joint <- function(cop, N = 1000) {
  
  U_sim <- BiCopSim(
    N = N,
    family = cop$family,
    par   = cop$par,
    par2  = cop$par2
  )
  
  wind_sim <- qgamma(
    U_sim[,1],
    shape = fit_wind$estimate["shape"],
    rate  = fit_wind$estimate["rate"]
  )
  
  price_sim <- qweibull(
    U_sim[,2],
    shape = fit_price$estimate["shape"],
    scale = fit_price$estimate["scale"]
  ) + min(price) - 1e-6
  
  data.frame(wind = wind_sim, price = price_sim)
}

# Generate simulations
sim_gauss <- simulate_joint(cop_gauss)
sim_bb8   <- simulate_joint(cop_fit)

# ---- define common tail thresholds ----
high_wind  <- quantile(c(sim_gauss$wind,  sim_bb8$wind),  0.90)
low_price  <- quantile(c(sim_gauss$price, sim_bb8$price), 0.10)

# Two‑panel copula outcome plot
library(ggplot2)
library(patchwork)

common_limits <- list(
  x = range(c(sim_gauss$wind, sim_bb8$wind)),
  y = range(c(sim_gauss$price, sim_bb8$price))
)

common_limits$y[2]<-min(common_limits$y[2],1500) # curtail for figure

p1 <- ggplot(sim_gauss, aes(wind, price)) +
  geom_point(alpha = 0.75, size = 0.4, colour = "grey30") +
  # highlight tail region (same in both panels)
  annotate(
    "rect",
    xmin = high_wind, xmax = Inf,
    ymin = -Inf, ymax = low_price,
    fill = "firebrick",
    alpha = 0.08
  ) +
  geom_density_2d(colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = common_limits$x) +
  scale_y_continuous(limits = common_limits$y) +
  labs(
    title = "(a) Gaussian copula",
    x = "Wind speed (m/s)",
    y = "Dispatch price (AUD/MWh)"
  ) +
  theme_bw()

p2 <- ggplot(sim_bb8, aes(wind, price)) +
  geom_point(alpha = 0.75, size = 0.4, colour = "steelblue") +
  annotate(
    "rect",
    xmin = high_wind, xmax = Inf,
    ymin = -Inf, ymax = low_price,
    fill = "firebrick",
    alpha = 0.08
  ) +
  geom_density_2d(colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = common_limits$x) +
  scale_y_continuous(limits = common_limits$y) +
  labs(
    title = "(b) BB8 copula",
    x = "Wind speed (m/s)",
    y = "Dispatch price (AUD/MWh)"
  ) +
  theme_bw()

p2 <- p2 +
  annotate(
    "text",
    x = Inf, y = -Inf,
    hjust = 1.1, vjust = -0.5,
    label = paste0(
      "P(high wind & low price) = ",
      round(P_bb8$cannibal, 3)
    ),
    size = 3
  )

p1 + p2


# ---- Figure 6 ----
# ============================================================
# 3D copula density plot (BB8)
# ============================================================
library(VineCopula)
library(plot3D)

# ---- grid in copula space ----
u <- seq(0.01, 0.99, length.out = 60)
v <- seq(0.01, 0.99, length.out = 60)

grid <- expand.grid(u = u, v = v)

# ---- copula density ----
dens_bb8 <- BiCopPDF(
  grid$u,
  grid$v,
  family = cop_fit$family,
  par    = cop_fit$par,
  par2   = cop_fit$par2
)

# reshape to matrix for plotting
dens_matrix <- matrix(dens_bb8, nrow = length(u), ncol = length(v))

# ---- density quantile thresholds ----
q50 <- quantile(dens_bb8, 0.50)
q90 <- quantile(dens_bb8, 0.90)
q95 <- quantile(dens_bb8, 0.95)


# ---- 3D surface (muted colours) ----
persp3D(
  x = u,
  y = v,
  z = dens_matrix,
  theta = 35,
  phi = 25,
  expand = 0.6,
  col = "grey85",          # neutral surface
  border = "grey60",
  ticktype = "simple",
  nticks = 4,
  xlab = "u (Wind)",
  ylab = "v (Price)",
  zlab = "Copula density",
  main = "BB8 Copula Density with Tail Quantile Contours"
)

# ---- add 95% contour (thin, dashed, dark grey) ----
contour3D(
  x = u,
  y = v,
  z = dens_matrix,
  colvar = dens_matrix,
  level = q95,
  add = TRUE,
  col = "black",
  lwd = 1.5,
  lty = 2,                # dashed
  drawlabels = FALSE
)

legend(
  "topleft",
  legend = c("95% density contour"),
  col = c("black"),
  lty = c(2, 1),
  lwd = c(2, 3),
  bty = "n"
)

# Spare figure
# Dynamic Version
library(rgl)

persp3d(
  x = u,
  y = v,
  z = dens_matrix,
  col = "lightblue",
  alpha = 0.9,
  xlab = "u (Wind)",
  ylab = "v (Price)",
  zlab = "Density"
)

### Section 6 Discussion - Hedging example:
# Indicator payoff: protection pays $1 when cannibalisation occurs
payoff_gauss <- as.numeric(
  wind_gauss > high_wind & price_gauss < low_price
)

payoff_bb8 <- as.numeric(
  wind_sim > high_wind & price_sim < low_price
)

price_gauss_contract <- mean(payoff_gauss)
price_bb8_contract   <- mean(payoff_bb8)

# Present results as relative mispricing
pricing_example <- data.frame(
  Model = c("Gaussian copula", "BB8 copula"),
  Contract_Value = c(price_gauss_contract, price_bb8_contract)
)

pricing_example$Relative_Value <-
  pricing_example$Contract_Value /
  pricing_example$Contract_Value[1]

# Hedge plot
ggplot(pricing_example, aes(Model, Contract_Value)) +
  geom_col(fill = "steelblue") +
  labs(
    y = "Model‑implied contract value",
    title = "Effect of dependence assumptions on cannibalisation protection pricing"
  ) +
  theme_bw()

### Different plot for payoff structure
# create payoff indicators
df_gauss_plot <- data.frame(
  wind = wind_gauss,
  price = price_gauss,
  payoff = wind_gauss > high_wind & price_gauss < low_price,
  model = "Gaussian copula"
)

df_bb8_plot <- data.frame(
  wind = wind_sim,
  price = price_sim,
  payoff = wind_sim > high_wind & price_sim < low_price,
  model = "BB8 copula"
)

df_payoff <- rbind(df_gauss_plot, df_bb8_plot)

# thin for plotting
set.seed(1)
df_payoff <- df_payoff[sample(nrow(df_payoff), 12000), ]

ggplot(df_payoff, aes(wind, price)) +
  geom_point(color = "grey70", alpha = 0.35, size = 0.6) +
  geom_point(
    data = subset(df_payoff, payoff),
    color = "firebrick",
    alpha = 0.6,
    size = 0.8
  ) +
  geom_vline(xintercept = high_wind, linetype = "dashed") +
  geom_hline(yintercept = low_price, linetype = "dashed") +
  facet_wrap(~model, ncol = 2) +
  labs(
    x = "Wind speed",
    y = "Dispatch price",
    title = "Payoff-relevant joint tail events under alternative dependence assumptions",
    subtitle = "Red points indicate simulated outcomes that trigger cannibalisation protection"
  ) +
  coord_cartesian(ylim = c(-900, 2000)) +
  theme_bw()

# function to compute cumulative payoff contributions
tail_curve <- function(wind, price, model) {
  q_seq <- seq(0.01, 0.10, by = 0.005)
  
  data.frame(
    quantile = q_seq,
    value = sapply(q_seq, function(q) {
      mean(wind > high_wind & price < quantile(price, q))
    }),
    model = model
  )
}

df_curve <- rbind(
  tail_curve(wind_gauss, price_gauss, "Gaussian copula"),
  tail_curve(wind_sim, price_sim, "BB8 copula")
)

ggplot(df_curve, aes(quantile, value, color = model)) +
  geom_line(linewidth = 1) +
  geom_point() +
  labs(
    x = "Lower price quantile threshold",
    y = "Cumulative contract value",
    title = "Contribution of joint tail events to contract valuation"
  ) +
  theme_bw()










