## ============================================================
## Pier 25 (USGS-01376520) Salinity Gap-Filling
## from Pier 84 (USGS-01376515) regression
## ============================================================

#needed packages
library(arrow)       
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(dataRetrieval)  
# =============================================================
#Reading in the continuous data from USGS sites using 
# package dataRetrieval. Can only grab so much at a time so 
# having to split it up for the download and then merge.
#==============================================================

#continuous data
#pier_25_sal_cont.1 <- read_waterdata_continuous(
#  monitoring_location_id = "USGS-01376520",
#  parameter_code = "90860",
#  time = c("2022-01-01", "2023-12-31")
#)
#pier_25_sal_cont.2 <- read_waterdata_continuous(
#  monitoring_location_id = "USGS-01376520",
#  parameter_code = "90860",
#  time = c("2024-01-01", "2025-12-31")
#)

#pier_25_sal_cont_full <- bind_rows(pier_25_sal_cont.1, pier_25_sal_cont.2)

#easier to read in later
#write_parquet(
#  pier_25_sal_cont_full,
#  "pier25_salinity_2022_25.parquet"
#)

#pier_84_sal_cont.1 <- read_waterdata_continuous(
#  monitoring_location_id = "USGS-01376515",
#  parameter_code = "90860",
#  time = c("2022-01-01", "2023-12-31")
#)
#pier_84_sal_cont.2 <- read_waterdata_continuous(
#  monitoring_location_id = "USGS-01376515",
#  parameter_code = "90860",
#  time = c("2024-01-01", "2025-12-31")
#)

#pier_84_sal_cont_full <- bind_rows(pier_84_sal_cont.1, pier_84_sal_cont.2)

#write_parquet(
#  pier_84_sal_cont_full,
#  "pier84_salinity_2022_25.parquet"
#)

# =======================================================================
# PARAMETERS read in from above and setting up parquet for output.
# setting up the minimum number of observations needed to trust the 
#regression as a safety catch. Setting up an interpolate residual 
# gap to bring to make sure there are no gaps over 1 hour that would
# be interpolated
#========================================================================

PARQUET_P25  <- "pier25_salinity_2022_25.parquet"
PARQUET_P84  <- "pier84_salinity_2022_25.parquet"
OUTPUT_FILE  <- "pier25_salinity_gapfilled.parquet"
MIN_OBS      <- 100   # minimum concurrent obs required to trust regression
INTERP_MAX   <- 4     # max gap length (# of 15-min steps = 1 hr) to linearly
                      # interpolate residual gaps; set to 0 to disable

#loading and checking the parquet files with stream gage data
p25_raw <- read_parquet(PARQUET_P25)
p84_raw <- read_parquet(PARQUET_P84)

# Inspect column names and adjust if parquets are different 
cat("Pier 25 columns:", paste(names(p25_raw), collapse = ", "), "\n") 
cat("Pier 84 columns:", paste(names(p84_raw), collapse = ", "), "\n")
#both match so all good

#=======================================================================
# Claude response for: Standardize to a clean datetime + value pair
# Adjust column names below if parquet has different field names.
#=======================================================================

clean_station <- function(df, sal_col = "value", time_col = "time",
                          approval_col = "result_approval") {
  out <- df |>
    rename(datetime = all_of(time_col),
           salinity = all_of(sal_col)) |>
    mutate(
      # Round to nearest 15-min to ensure clean join spine
      datetime = floor_date(as.POSIXct(datetime, tz = "America/New_York"),
                            unit = "15 minutes"),
      salinity = as.numeric(salinity)
    )

  # If an approval/grade column exists, flag provisional data
  if (approval_col %in% names(df)) {
    out <- out |>
      rename(approval = all_of(approval_col))
  }

  out |>
    select(datetime, salinity, any_of("approval")) |>
    filter(!is.na(salinity)) |>
    # Keep one value per timestamp if duplicates exist (take mean)
    group_by(datetime) |>
    summarise(salinity = mean(salinity, na.rm = TRUE), .groups = "drop") |>
    arrange(datetime)
}

p25 <- clean_station(p25_raw)
p84 <- clean_station(p84_raw)

cat("\nPier 25: ", nrow(p25), "valid 15-min readings\n")
cat("Pier 84: ", nrow(p84), "valid 15-min readings\n")
#Pier 25:  106802 valid 15-min readings
#Pier 84:  91419 valid 15-min readings
#===================================================================
# Build a full 15-min time spine 
#===================================================================
t_start <- min(c(p25$datetime, p84$datetime))
t_end   <- max(c(p25$datetime, p84$datetime))

spine <- tibble(
  datetime = seq(t_start, t_end, by = "15 min")
)
#check length of time spline
cat("\nFull spine:", nrow(spine), "15-min steps (",
    round(nrow(spine) / (4 * 24 * 365.25), 1), "years)\n")
#Full spine: 139807 15-min steps ( 4 years)

# oin both stations onto spine 
combined <- spine |>
  left_join(p25 |> rename(sal_p25 = salinity), by = "datetime") |>
  left_join(p84 |> rename(sal_p84 = salinity), by = "datetime")

# Summarise data availability
avail <- combined |>
  summarise(
    n_total          = n(),
    n_p25_present    = sum(!is.na(sal_p25)),
    n_p84_present    = sum(!is.na(sal_p84)),
    n_both_present   = sum(!is.na(sal_p25) & !is.na(sal_p84)),
    n_p25_gap        = sum(is.na(sal_p25)),
    n_p84_can_fill   = sum(is.na(sal_p25) & !is.na(sal_p84)),
    pct_p25_missing  = round(100 * sum(is.na(sal_p25)) / n(), 1)
  )

cat("\n=== Data Availability ===\n")
print(avail)


#==================================================================
#Identify gap periods in Pier 25 to be filled
#==================================================================
gap_runs <- combined |>
  mutate(
    is_gap = is.na(sal_p25),
    run_id = cumsum(is_gap != lag(is_gap, default = FALSE))
  ) |>
  filter(is_gap) |>
  group_by(run_id) |>
  summarise(
    gap_start    = min(datetime),
    gap_end      = max(datetime),
    gap_steps    = n(),
    gap_hours    = round(n() / 4, 1),
    p84_coverage = round(100 * mean(!is.na(sal_p84)), 1),
    .groups      = "drop"
  ) |>
  arrange(desc(gap_steps))

cat("\n=== Pier 25 Gap Periods (longest first) ===\n")
print(gap_runs, n = 30)
#=================================================================
#Regression: sal_p25 ~ sal_p84 on concurrent data 
#=================================================================

#setting up training data where concurrent availible
train <- combined |>
  filter(!is.na(sal_p25), !is.na(sal_p84))

cat("\n=== Regression Training Set ===\n")
cat("Concurrent observations:", nrow(train), "\n")
# 60924
#making sure over set limit earlier
if (nrow(train) < MIN_OBS) {
  stop(paste("Only", nrow(train), "concurrent observations —",
             "not enough to build a reliable regression. Check your data."))
}


mod <- lm(sal_p25 ~ sal_p84, data = train)
r2   <- summary(mod)$r.squared
rmse <- sqrt(mean(residuals(mod)^2))

cat("\n=== Regression: sal_p25 ~ sal_p84 ===\n")
cat("Intercept:", round(coef(mod)[1], 4), "\n")
cat("Slope:    ", round(coef(mod)[2], 4), "\n")
cat("R²:       ", round(r2, 4), "\n")
cat("RMSE:     ", round(rmse, 4), "PSU\n")

#=================================================================
#Look at plots to ensure regression appropriate
#=================================================================

## Scatter: concurrent observations
p_scatter <- ggplot(train, aes(x = sal_p84, y = sal_p25)) +
  geom_hex(bins = 60) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.7) +
  geom_smooth(method = "lm", color = "firebrick", se = TRUE) +
  scale_fill_viridis_c(option = "magma") +
  labs(
    title = "Pier 25 vs Pier 84 — Concurrent Salinity (15-min)",
    subtitle = paste0("n = ", nrow(train),
                      " | R² = ", round(r2, 3),
                      " | RMSE = ", round(rmse, 3), " PSU"),
    x = "Pier 84 Salinity (PSU)",
    y = "Pier 25 Salinity (PSU)"
  ) +
  theme_bw()

#ggsave("pier25_84_scatter.png", p_scatter, width = 7, height = 5.5, dpi = 150)

## Residuals over time to check for temporal drift
train_resid <- train |>
  mutate(residual = residuals(mod))

p_resid <- ggplot(train_resid, aes(x = datetime, y = residual)) +
  geom_point(alpha = 0.1, size = 0.4, color = "steelblue") +
  geom_smooth(method = "loess", span = 0.1, color = "firebrick", se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Regression Residuals Over Time",
    subtitle = "Systematic drift would indicate the pier relationship changes seasonally",
    x = NULL, y = "Residual (PSU)"
  ) +
  theme_bw()

#ggsave("pier25_84_residuals.png", p_resid, width = 9, height = 4, dpi = 150)

#=========================================================================
# Apply fill from regression
#=========================================================================

filled <- combined |>
  mutate(
    # Predict sal_p25 from sal_p84 using the regression
    sal_p25_pred = predict(mod, newdata = tibble(sal_p84 = sal_p84)),

    # Fill: use observed where available, regression where gap + p84 present
    sal_filled  = case_when(
      !is.na(sal_p25)                    ~ sal_p25,
      is.na(sal_p25) & !is.na(sal_p84)  ~ sal_p25_pred,
      TRUE                               ~ NA_real_
    ),

    fill_source = case_when(
      !is.na(sal_p25)                    ~ "observed_p25",
      is.na(sal_p25) & !is.na(sal_p84)  ~ "regression_p84",
      TRUE                               ~ "still_missing"
    )
  )

#Fill summary 
fill_summary <- filled |>
  count(fill_source) |>
  mutate(
    pct       = round(100 * n / nrow(filled), 1),
    est_hours = round(n / 4, 0)
  )

print(fill_summary)

# ---- 11. Time series plot of filled record ----------------------------------
p_ts <- filled |>
  slice_sample(n = min(nrow(filled), 50000)) |>  # downsample for plot speed
  ggplot(aes(x = datetime, y = sal_filled, color = fill_source)) +
  geom_point(size = 0.3, alpha = 0.4) +
  scale_color_manual(
    values = c(
      "observed_p25"   = "steelblue",
      "regression_p84" = "firebrick",
      "still_missing"  = "grey80",
      paste0("interpolated_<=", INTERP_MAX * 15, "min") = "darkorange"
    )
  ) +
  labs(
    title = "Pier 25 Salinity — Observed and Gap-Filled",
    x = NULL, y = "Salinity (PSU)", color = "Source"
  ) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

#ggsave("pier25_filled_timeseries.png", p_ts, width = 12, height = 4.5, dpi = 150)

#Export
filled |>
  select(datetime, sal_p25_observed = sal_p25, sal_p84 = sal_p84,
         sal_p25_filled = sal_filled, fill_source) |>
  write_parquet(OUTPUT_FILE)


# getting single daily mean salinity value, max, min, and other info for

filled <- read_parquet("pier25_salinity_gapfilled.parquet")

daily_salinity <- filled |>
  mutate(Date = as.Date(datetime, tz = "America/New_York")) |>
  group_by(Date) |>
  summarise(
    sal_mean        = mean(sal_p25_filled, na.rm = TRUE),
    sal_sd          = sd(sal_p25_filled, na.rm = TRUE),
    sal_min         = min(sal_p25_filled, na.rm = TRUE),
    sal_max         = max(sal_p25_filled, na.rm = TRUE),
    n_obs           = sum(!is.na(sal_p25_filled)),        # number of 15-min readings used
    n_observed      = sum(fill_source == "observed", na.rm = TRUE),
    n_regression    = sum(fill_source == "regression_p84", na.rm = TRUE),
    n_interpolated  = sum(fill_source == "interpolated", na.rm = TRUE),
    pct_observed    = round(100 * n_observed / n_obs, 1),
    day_complete    = n_obs >= 88,   # TRUE if ≥ 88 of 96 possible readings present (>= 91.7%)
    .groups = "drop"
  ) |>
  # Set daily mean to NA if fewer than half the day's readings are available
  mutate(
    sal_mean = if_else(n_obs < 48, NA_real_, sal_mean)
  )
write_parquet(daily_salinity, "daily_mean_sal.parquet")

cat("=== Daily Salinity Summary ===\n")
cat("Days total:        ", nrow(daily_salinity), "\n")
cat("Days with data:    ", sum(!is.na(daily_salinity$sal_mean)), "\n")
cat("Days missing:      ", sum(is.na(daily_salinity$sal_mean)), "\n")
cat("Days 100% observed:", sum(daily_salinity$pct_observed == 100, na.rm = TRUE), "\n")
cat("Mean salinity:     ", round(mean(daily_salinity$sal_mean, na.rm = TRUE), 3), "PSU\n")

print(head(daily_salinity, 10))
