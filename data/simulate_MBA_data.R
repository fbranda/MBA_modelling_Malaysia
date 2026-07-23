# Simulated dataset: data_MBA_modelling
library(tidyverse)
set.seed(42)
n <- 2000

# 1) Core functions

# Parameters of the 2-component log-normal mixture for one variable
mix_params <- function(min_val, q1, median_val, mean_val, p_high, sdlog_high) {
  meanlog_low  <- log(median_val)
  sdlog_low    <- (log(median_val) - log(q1)) / 0.6745
  mean_low     <- exp(meanlog_low + sdlog_low^2 / 2)
  mean_high    <- (mean_val - (1 - p_high) * mean_low) / p_high
  mean_high    <- max(mean_high, mean_low * 2)
  meanlog_high <- log(mean_high) - sdlog_high^2 / 2
  list(meanlog_low = meanlog_low, sdlog_low = sdlog_low,
       meanlog_high = meanlog_high, sdlog_high = sdlog_high)
}

# Quantile function (inverse CDF) of the mixture, evaluated at a vector of
# uniforms u (which may themselves be correlated across variables)
qmix <- function(u, min_val, max_val, meanlog_low, sdlog_low,
                 meanlog_high, sdlog_high, p_high, n_grid = 20000) {
  lo_grid  <- max(min_val, 1e-6)
  x_grid   <- exp(seq(log(lo_grid), log(max_val), length.out = n_grid))
  cdf_grid <- (1 - p_high) * plnorm(x_grid, meanlog_low, sdlog_low) +
    p_high  * plnorm(x_grid, meanlog_high, sdlog_high)
  x <- approx(cdf_grid, x_grid, xout = u, rule = 2)$y
  pmin(pmax(x, min_val), max_val)
}

# Independent draw (used only for pop.density)
simulate_skewed_var <- function(n, min_val, q1, median_val, mean_val, q3, max_val,
                                p_high = 0.05, sdlog_high = 1.0) {
  p <- mix_params(min_val, q1, median_val, mean_val, p_high, sdlog_high)
  x <- qmix(runif(n), min_val, max_val, p$meanlog_low, p$sdlog_low,
            p$meanlog_high, p$sdlog_high, p_high)
  x[which.min(x)] <- min_val
  x[which.max(x)] <- max_val
  x
}

# Mildly-skewed bounded variables (coordinates, age) - unchanged
simulate_bounded_var <- function(n, min_val, q1, mean_val, q3, max_val) {
  range_val <- max_val - min_val
  mu <- (mean_val - min_val) / range_val
  mu <- min(max(mu, 0.02), 0.98)
  iqr_norm <- (q3 - q1) / range_val
  sd_norm  <- iqr_norm / 1.349
  var_norm <- sd_norm^2
  k <- mu * (1 - mu) / var_norm - 1
  k <- max(k, 2)
  a <- mu * k
  b <- (1 - mu) * k
  x01 <- rbeta(n, a, b)
  x <- min_val + x01 * range_val
  x[which.min(x)] <- min_val
  x[which.max(x)] <- max_val
  x
}

# 2) Parameters for skewed variables: 12 antigens + population density

skewed_specs <- data.frame(
  variable = c("Toxo.SAG2A", "Giardia.VSP3", "Giardia.VSP5", "Lf.Brugia.BmR1",
               "Lf.Wuchereria.Wb123", "Lf.Brugia.Bm14", "Lf.Brugia.Bm33",
               "Trachoma.pgp3", "Trachoma.ct694", "Strongyloides.NIE",
               "Yaws.rp17", "Yaws.TmpA", "pop.density"),
  min    = c(0.00,     0.0,    0.00,    0.0,     0.00,   0.0,   0.0,
             0.0,      0.0,    0.0,     0.50,    0.0,    0.0),
  q1     = c(24.50,   45.0,   43.00,  258.5,    63.50,  66.0, 144.2,
             51.0,   551.4,   91.5,    48.19,   65.0,   13.0),
  median = c(58.75,  107.0,   88.62,  438.5,    90.71, 101.0, 206.2,
             122.0,   955.0,  140.0,    75.00,  102.5,   34.0),
  mean   = c(1582.74, 544.7,  498.69, 1022.9,   121.80, 309.4, 539.0,
             3427.9,  2501.8,  779.0,   380.85, 258.7,  184.8),
  q3     = c(559.12,  340.0,  277.19,  784.1,   134.50, 169.5, 344.0,
             898.3,  1887.5,  263.0,   127.00, 195.8,  203.0),
  max    = c(33418.00, 32343.8, 32470.75, 43905.0, 11646.00, 30413.5, 34534.8,
             45446.0, 40468.5, 38823.0,  31837.00, 18641.0,  1259.0),
  p_high     = c(0.2167, 0.0723, 0.1186, 0.0319, 0.0627, 0.0501, 0.0911,
                 0.1954, 0.0592, 0.1042, 0.0472, 0.0914, 0.1527),
  sdlog_high = c(0.3398, 0.6740, 0.6029, 1.5805, 1.2020, 0.9482, 0.7519,
                 0.8360, 0.9152, 1.0435, 0.7719, 0.4572, 0.4729),
  stringsAsFactors = FALSE
)

# 3) Target correlation matrix between the 12 antigens (read from Fig. S1)

ab_vars <- c("Toxo.SAG2A", "Giardia.VSP3", "Giardia.VSP5", "Lf.Brugia.BmR1",
             "Lf.Wuchereria.Wb123", "Lf.Brugia.Bm14", "Lf.Brugia.Bm33",
             "Trachoma.pgp3", "Trachoma.ct694", "Strongyloides.NIE",
             "Yaws.rp17", "Yaws.TmpA")

R <- matrix(c(
  1.00,0.14,0.12,0.14,0.18,0.16,0.15,0.11,0.12,0.17,0.15,0.14,
  0.14,1.00,0.84,0.20,0.26,0.22,0.21,0.12,0.13,0.19,0.19,0.17,
  0.12,0.84,1.00,0.20,0.26,0.19,0.19,0.11,0.13,0.18,0.19,0.17,
  0.14,0.20,0.20,1.00,0.39,0.43,0.40,0.16,0.26,0.25,0.30,0.26,
  0.18,0.26,0.26,0.39,1.00,0.42,0.39,0.16,0.24,0.36,0.35,0.37,
  0.16,0.22,0.19,0.43,0.42,1.00,0.42,0.20,0.28,0.26,0.30,0.26,
  0.15,0.21,0.19,0.40,0.39,0.42,1.00,0.17,0.20,0.26,0.29,0.27,
  0.11,0.12,0.11,0.16,0.16,0.20,0.17,1.00,0.67,0.17,0.25,0.15,
  0.12,0.13,0.13,0.26,0.24,0.28,0.20,0.67,1.00,0.15,0.28,0.22,
  0.17,0.19,0.18,0.25,0.36,0.26,0.26,0.17,0.15,1.00,0.28,0.24,
  0.15,0.19,0.19,0.30,0.35,0.30,0.29,0.25,0.28,0.28,1.00,0.38,
  0.14,0.17,0.17,0.26,0.37,0.26,0.27,0.15,0.22,0.24,0.38,1.00
), nrow = 12, byrow = TRUE, dimnames = list(ab_vars, ab_vars))

# A correlation matrix read off a plot and rounded to 2 decimals is
# sometimes not exactly positive definite: project it onto the nearest
# valid correlation matrix before using it in mvrnorm()
R_pd <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)

# 4) Correlated generation of the 12 antigens via Gaussian copula

Z <- MASS::mvrnorm(n, mu = rep(0, 12), Sigma = R_pd)
U <- pnorm(Z)                     # correlated uniforms, one column per antigen
colnames(U) <- ab_vars

antigen_list <- lapply(ab_vars, function(v) {
  s <- skewed_specs[skewed_specs$variable == v, ]
  p <- mix_params(s$min, s$q1, s$median, s$mean, s$p_high, s$sdlog_high)
  x <- qmix(U[, v], s$min, s$max, p$meanlog_low, p$sdlog_low,
            p$meanlog_high, p$sdlog_high, s$p_high)
  x[which.min(x)] <- s$min
  x[which.max(x)] <- s$max
  x
})
names(antigen_list) <- ab_vars
antigen_df <- as.data.frame(antigen_list)

# pop.density: no reported correlation with the antigens, so it is still
# generated independently, as in the original version
s <- skewed_specs[skewed_specs$variable == "pop.density", ]
pop.density <- simulate_skewed_var(n, s$min, s$q1, s$median, s$mean, s$q3, s$max,
                                   p_high = s$p_high, sdlog_high = s$sdlog_high)

skewed_df <- cbind(antigen_df, pop.density)

skewed_df[ab_vars]    <- round(skewed_df[ab_vars], 2)
skewed_df$pop.density <- round(skewed_df$pop.density)

# 5) Coordinates and age (unchanged)

x   <- round(simulate_bounded_var(n, 447573, 468594, 484645, 503522, 534260))
y   <- round(simulate_bounded_var(n, 652927, 676045, 721060, 749936, 813749))
age <- round(simulate_bounded_var(n, 0.2164, 10.0000, 29.2249, 45.4212, 105.0000), 4)

# 6) Categorical variables (unchanged)

wealth <- factor(
  sample(c("Low-income", "Lower-middle-income", "Middle-income", "Wealthy"),
         n, replace = TRUE, prob = c(1330, 1594, 1600, 1644)),
  levels = c("Low-income", "Lower-middle-income", "Middle-income", "Wealthy")
)

occupation <- factor(
  sample(c("housewife", "none", "other", "outside activities", "student"),
         n, replace = TRUE, prob = c(1166, 1776, 401, 1173, 1652)),
  levels = c("housewife", "none", "other", "outside activities", "student")
)

gender <- factor(
  sample(c("Female", "Male"), n, replace = TRUE, prob = c(3240, 2928)),
  levels = c("Female", "Male")
)

veg.stratum <- factor(
  sample(c("Dense-vegetation", "Moderate-vegetation", "Sparse-vegetation"),
         n, replace = TRUE, prob = c(1789, 2276, 2103)),
  levels = c("Dense-vegetation", "Moderate-vegetation", "Sparse-vegetation")
)


# 7) Final assembly (20 columns, same order/names as the original)

data_MBA_modelling <- data.frame(
  skewed_df,
  x = x, y = y, age = age,
  wealth = wealth, occupation = occupation,
  gender = gender, veg.stratum = veg.stratum
)

data_MBA_modelling$Units <- seq(1, dim(data_MBA_modelling)[1], 1)
data_MBA_modelling <- data_MBA_modelling %>% relocate(Units, .before = 1)

write.csv(data_MBA_modelling, "data_MBA_modelling.csv", row.names = FALSE)
