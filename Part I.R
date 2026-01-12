############################################################
# Case study: Scheduling MRI facilities — PART I (R)
# Statistical analysis + bootstrap uncertainty
# Data file: ScanRecords.csv
############################################################

rm(list=ls())

# Packages
install.packages(c(
  "dplyr",
  "tidyr",
  "lubridate",
  "readr",
  "purrr",
  "MASS"
), repos='http://cran.us.r-project.org')

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(readr)
  library(purrr)
  library(MASS)
})

# ----------------------------
# 0) Settings you can choose
# ----------------------------
setwd("Git/MRI-Scan-Project-Computational")
DATA_PATH <- "ScanRecords.csv"

B <- 5000                 # bootstrap replications
set.seed(1)

# Service level for slot length:
# e.g. p_slot=0.90 means "90% of scans finish within the slot"
p_slot <- 0.90

# Round slots to nearest 5 minutes (hospital-friendly)
round_minutes <- 5

# ----------------------------
# 1) Load 
# ----------------------------
df <- read_csv(DATA_PATH, show_col_types = FALSE) %>%
  mutate(
    Date = as.Date(Date),
    Time = as.numeric(Time),
    Duration = as.numeric(Duration),
    PatientType = as.factor(PatientType)
  )

stopifnot(all(c("Date","Time","Duration","PatientType") %in% names(df)))

# Keep only working-hour calls (the case says calls 8:00–17:00 on working days)
df <- df %>% filter(Time >= 8, Time <= 17)

day_index <- df %>%
  distinct(Date) %>%
  arrange(Date) %>%
  mutate(day_id = row_number())

df <- df %>%
  left_join(day_index, by="Date") %>%
  mutate(
    # working-hours since start of first day (0..), day length = 9 hours
    t_work = (day_id - 1) * 9 + (Time - 8)
  )

# Split types
df1 <- df %>% filter(PatientType == "Type 1")
df2 <- df %>% filter(PatientType == "Type 2")

# ----------------------------
# 2) Helper functions
# ----------------------------

# Round hours to nearest X minutes (default 5 minutes)
round_slot_hours <- function(hours, minutes=5) {
  step <- minutes / 60
  step * round(hours / step)
}

# Bootstrap CI helper (percentile CI)
ci_pct <- function(x, level=0.95) {
  alpha <- (1 - level)/2
  unname(quantile(x, probs=c(alpha, 1-alpha), na.rm=TRUE))
}

# Daily arrival counts for a given data frame
daily_counts <- function(dfx) {
  dfx %>%
    count(Date, name="n") %>%
    arrange(Date)
}

# Interarrival times on the working-time axis (hours)
interarrival <- function(dfx) {
  dfx %>%
    arrange(t_work) %>%
    mutate(delta = t_work - lag(t_work)) %>%
    filter(!is.na(delta)) %>%
    pull(delta)
}

# ----------------------------
# 3) TYPE 1: arrivals/day ~ Poisson(lambda), durations ~ Normal(mu, sigma)
# ----------------------------

# ---- Arrivals/day: estimate lambda (mean arrivals/day)
c1 <- daily_counts(df1)
lambda_hat <- mean(c1$n)
D1 <- nrow(c1)

# Parametric bootstrap for lambda via "days of Poisson counts"
boot_lambda <- replicate(B, {
  mean(rpois(D1, lambda_hat))
})
lambda_ci <- ci_pct(boot_lambda)

# ---- Durations: estimate mu, sigma, and slot quantile q_p
dur1 <- df1$Duration
n1 <- length(dur1)
mu_hat <- mean(dur1)
sd_hat <- sd(dur1)

# Parametric bootstrap for (mu, sd, q_p)
boot_mu <- numeric(B)
boot_sd <- numeric(B)
boot_qp <- numeric(B)

for (b in 1:B) {
  sim <- rnorm(n1, mean=mu_hat, sd=sd_hat)
  boot_mu[b] <- mean(sim)
  boot_sd[b] <- sd(sim)
  boot_qp[b] <- unname(quantile(sim, probs=p_slot))
}

mu_ci <- ci_pct(boot_mu)
sd_ci <- ci_pct(boot_sd)
q1_hat <- unname(quantile(dur1, probs=p_slot))
q1_ci <- ci_pct(boot_qp)

slot1_hours <- round_slot_hours(q1_hat, round_minutes)

risk1_hat <- 1 - pnorm(slot1_hours, mean=mu_hat, sd=sd_hat)
boot_risk1 <- 1 - pnorm(slot1_hours, mean=boot_mu, sd=boot_sd)
risk1_ci <- ci_pct(boot_risk1)

# ----------------------------
# 4) TYPE 2: unknown arrivals + unknown durations
#    Use nonparametric bootstrap ("resample what we saw")
# ----------------------------

# ---- Arrivals/day: treat daily counts as i.i.d. but distribution unknown
c2 <- daily_counts(df2)
D2 <- nrow(c2)
arr2_mean_hat <- mean(c2$n)

boot_arr2_mean <- replicate(B, {
  mean(sample(c2$n, size=D2, replace=TRUE))
})
arr2_mean_ci <- ci_pct(boot_arr2_mean)

# ---- Durations: unknown distribution => plug-in estimators + bootstrap uncertainty
dur2 <- df2$Duration
n2 <- length(dur2)

dur2_mean_hat <- mean(dur2)
dur2_median_hat <- median(dur2)
q2_hat <- unname(quantile(dur2, probs=p_slot))
dur2_sd_hat <- sd(dur2)

boot_dur2_mean <- numeric(B)
boot_dur2_median <- numeric(B)
boot_dur2_qp <- numeric(B)
boot_dur2_sd <- numeric(B)

for (b in 1:B) {
  samp <- sample(dur2, size=n2, replace=TRUE)
  boot_dur2_mean[b] <- mean(samp)
  boot_dur2_median[b] <- median(samp)
  boot_dur2_qp[b] <- unname(quantile(samp, probs=p_slot))
  boot_dur2_sd[b] <- sd(samp)
}

dur2_mean_ci <- ci_pct(boot_dur2_mean)
dur2_median_ci <- ci_pct(boot_dur2_median)
q2_ci <- ci_pct(boot_dur2_qp)

slot2_hours <- round_slot_hours(q2_hat, round_minutes)
risk2_hat <- mean(dur2 > slot2_hours)
boot_risk2 <- replicate(B, {
  mean(sample(dur2, size=n2, replace=TRUE) > slot2_hours)
})
risk2_ci <- ci_pct(boot_risk2)

# ECDF for type 2 scan
ecdf_dur2 <- ecdf(dur2)
ecdf_arr2 <- ecdf(c2$n)

# plotting ecdf
plot(ecdf_dur2,
     main = "ECDF - Type 2 Scan Durations",
     xlab = "Duration (hours)",
     ylab = "F(t)")

plot(ecdf_arr2,
     main = "ECDF - Type 2 Daily Arrivals",
     xlab = "Number of patients per day",
     ylab = "F(n)")

# ----------------------------
# 5) Interarrival times (optional, but useful for Part II simulation)
# ----------------------------
ia1 <- interarrival(df1)  # should look exponential-ish if Poisson process
ia2 <- interarrival(df2)  # unknown

ia1_mean_hat <- mean(ia1)
ia2_mean_hat <- mean(ia2)

# Bootstrap CI for mean interarrival 
boot_ia1_mean <- replicate(B, mean(sample(ia1, replace=TRUE)))
boot_ia2_mean <- replicate(B, mean(sample(ia2, replace=TRUE)))
ia1_mean_ci <- ci_pct(boot_ia1_mean)
ia2_mean_ci <- ci_pct(boot_ia2_mean)

# ----------------------------
# 6) Output summary (print + objects for report)
# ----------------------------

cat("\n==================== PART I RESULTS ====================\n")

cat("\n--- TYPE 1 (known families: Poisson arrivals, Normal durations) ---\n")
cat(sprintf("Arrivals/day lambda_hat = %.3f | 95%% boot CI [%.3f, %.3f]\n",
            lambda_hat, lambda_ci[1], lambda_ci[2]))
cat(sprintf("Duration mu_hat = %.3f | 95%% boot CI [%.3f, %.3f]\n",
            mu_hat, mu_ci[1], mu_ci[2]))
cat(sprintf("Duration sd_hat = %.3f | 95%% boot CI [%.3f, %.3f]\n",
            sd_hat, sd_ci[1], sd_ci[2]))
cat(sprintf("p=%.2f quantile of duration q_hat = %.3f | 95%% boot CI [%.3f, %.3f]\n",
            p_slot, q1_hat, q1_ci[1], q1_ci[2]))
cat(sprintf("Recommended slot length (rounded to %d min): %.3f hours (%.0f minutes)\n",
            round_minutes, slot1_hours, slot1_hours*60))
cat(sprintf("Risk(Duration > slot) under fitted Normal: %.3f | 95%% boot CI [%.3f, %.3f]\n",
            risk1_hat, risk1_ci[1], risk1_ci[2]))

cat("\n--- TYPE 2 (unknown: nonparametric plug-in + bootstrap) ---\n")
cat(sprintf("Arrivals/day mean = %.3f | 95%% boot CI [%.3f, %.3f]\n",
            arr2_mean_hat, arr2_mean_ci[1], arr2_mean_ci[2]))
cat(sprintf("Duration mean = %.3f | 95%% boot CI [%.3f, %.3f]\n",
            dur2_mean_hat, dur2_mean_ci[1], dur2_mean_ci[2]))
cat(sprintf("Duration median = %.3f | 95%% boot CI [%.3f, %.3f]\n",
            dur2_median_hat, dur2_median_ci[1], dur2_median_ci[2]))
cat(sprintf("p=%.2f quantile of duration q_hat = %.3f | 95%% boot CI [%.3f, %.3f]\n",
            p_slot, q2_hat, q2_ci[1], q2_ci[2]))
cat(sprintf("Recommended slot length (rounded to %d min): %.3f hours (%.0f minutes)\n",
            round_minutes, slot2_hours, slot2_hours*60))
cat(sprintf("Risk(Duration > slot) empirical: %.3f | 95%% boot CI [%.3f, %.3f]\n",
            risk2_hat, risk2_ci[1], risk2_ci[2]))

cat("\n--- Interarrival means (working-time axis, optional) ---\n")
cat(sprintf("Type 1 mean interarrival (hours): %.3f | 95%% boot CI [%.3f, %.3f]\n",
            ia1_mean_hat, ia1_mean_ci[1], ia1_mean_ci[2]))
cat(sprintf("Type 2 mean interarrival (hours): %.3f | 95%% boot CI [%.3f, %.3f]\n",
            ia2_mean_hat, ia2_mean_ci[1], ia2_mean_ci[2]))

# ----------------------------
# 7) Monte Carlo simulation to check type II
# ----------------------------
M <- 1000 # monte carlo replications
B_MC <- 199

# existing stats from the boostrapping of Type 2
dur2_mean_hat <- mean(dur2)
q2_hat <- unname(quantile(dur2, probs=p_slot))
risk2_hat <- mean(dur2 > slot2_hours)

# storing true values that fall in the stats
mean_in <- logical(M)
quantile_in <- logical(M)
risk_in <- logical(M)

# store MC for bias estimates
mc_mean <- numeric(M)
mc_quantile <- numeric(M)
mc_risk <- numeric(M)

# sample generating loop
for (m in 1:M) {
  new_sample <- sample(dur2, size = n2, replace = TRUE) # assumes ecdf as the real
  
  # montecarlo stats
  mc_mean[m] <- mean(new_sample)
  mc_quantile[m] <- unname(quantile(new_sample, probs = p_slot))
  mc_risk[m] <- mean(new_sample > slot2_hours)
  
  # store MC boot stats
  mcboot_mean <- numeric(B_MC)
  mcboot_quantile <- numeric(B_MC)
  mcboot_risk <- numeric(B_MC)
  
  for (b in 1:B_MC) {
    newsample_boot <- sample(new_sample, size = n2, replace = TRUE)
    
    # boot stats from resampling
    mcboot_mean[b] <- mean(newsample_boot)
    mcboot_quantile[b]  <- quantile(newsample_boot, p_slot)
    mcboot_risk[b] <- mean(newsample_boot > slot2_hours)
  }
  
  # confidence intervals
  mcci_mean <- ci_pct(mcboot_mean)
  mcci_quantile <- ci_pct(mcboot_quantile)
  mcci_risk <- ci_pct(mcboot_risk)
  
  #do the means fall into the real range
  mean_in[m] <- (mcci_mean[1] <= dur2_mean_hat && dur2_mean_hat <= mcci_mean[2])
  quantile_in[m]   <- (mcci_quantile[1] <= q2_hat && q2_hat <= mcci_quantile[2])
  risk_in[m] <- (mcci_risk[1] <= risk2_hat && risk2_hat<= mcci_risk[2])
}

# printing results
cat(sprintf("Mean: bias = %.4f\n", mean(mc_mean) - dur2_mean_hat))
cat(sprintf("Q%.0f: bias = %.4f\n", 100*p_slot, mean(mc_quantile) - q2_hat))
cat(sprintf("Risk: bias = %.4f\n", mean(mc_risk) - risk2_hat))


cat("\n=== Bootstrap CI coverage under ECDF model for Type 2 patients ===\n")
cat(sprintf("Target coverage: 0.95 | M=%d, B_in=%d\n", M, B_MC))
cat(sprintf("Mean coverage: %.3f\n", mean(mean_in)))
cat(sprintf("Q%.0f coverage: %.3f\n", 100*p_slot, mean(quantile_in)))
cat(sprintf("Risk coverage: %.3f\n", mean(risk_in)))


# ----------------------------
# 8) Export bootstrapped parameters to CSV for DES
# ----------------------------

# Create a data frame with all bootstrapped parameter values
# Each row is one bootstrap replicate for probabilistic sensitivity analysis

des_params <- data.frame(
  replicate = 1:B,

  # Type 1 parameters (Poisson arrivals, Normal durations)
  t1_lambda = boot_lambda,
  t1_duration_mean = boot_mu,
  t1_duration_sd = boot_sd,
  t1_duration_q90 = boot_qp,
  t1_overtime_risk = boot_risk1,

  # Type 2 parameters (nonparametric)
  t2_arrivals_mean = boot_arr2_mean,
  t2_duration_mean = boot_dur2_mean,
  t2_duration_median = boot_dur2_median,
  t2_duration_q90 = boot_dur2_qp,
  t2_overtime_risk = boot_risk2,

  # Interarrival times
  t1_interarrival_mean = boot_ia1_mean,
  t2_interarrival_mean = boot_ia2_mean
)

# Export to CSV
write_csv(des_params, "des_bootstrap_params.csv")

# Also create a summary table with point estimates and CIs
des_summary <- data.frame(
  parameter = c(
    "t1_lambda", "t1_duration_mean", "t1_duration_sd", "t1_duration_q90",
    "t1_slot_hours", "t1_overtime_risk", "t1_interarrival_mean",
    "t2_arrivals_mean", "t2_duration_mean", "t2_duration_median",
    "t2_duration_q90", "t2_slot_hours", "t2_overtime_risk", "t2_interarrival_mean"
  ),
  estimate = c(
    lambda_hat, mu_hat, sd_hat, q1_hat, slot1_hours, risk1_hat, ia1_mean_hat,
    arr2_mean_hat, dur2_mean_hat, dur2_median_hat, q2_hat, slot2_hours, risk2_hat, ia2_mean_hat
  ),
  ci_lower = c(
    lambda_ci[1], mu_ci[1], sd_ci[1], q1_ci[1], NA, risk1_ci[1], ia1_mean_ci[1],
    arr2_mean_ci[1], dur2_mean_ci[1], dur2_median_ci[1], q2_ci[1], NA, risk2_ci[1], ia2_mean_ci[1]
  ),
  ci_upper = c(
    lambda_ci[2], mu_ci[2], sd_ci[2], q1_ci[2], NA, risk1_ci[2], ia1_mean_ci[2],
    arr2_mean_ci[2], dur2_mean_ci[2], dur2_median_ci[2], q2_ci[2], NA, risk2_ci[2], ia2_mean_ci[2]
  )
)

write_csv(des_summary, "des_params_summary.csv")

cat("\n--- Exported DES parameters ---\n")
cat(sprintf("Bootstrap samples (B=%d): des_bootstrap_params.csv\n", B))
cat("Summary with point estimates and 95%% CIs: des_params_summary.csv\n")

