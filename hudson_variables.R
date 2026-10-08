# =============================================================================
# Hudson River Harbor — Environmental Covariate Data Frame
# Study period: 2022-03-30 to 2025-10-16 (for now, possilby more detections later)
# =============================================================================
# Variables included:
#   1. Water temperature — NOAA Battery (°C, daily mean)
#   2. River discharge — USGS Troy (m³/s, 4-day lag, daily mean)
#   3. Tidal metrics — NOAA Battery MHHW (m, daily max/min/range)
#   4. Lunar phase — lunar package
#   5. Photoperiod — based on coordinate latitude in cneter of Upper Bay
#   6. Salinity — Pier 25 HRECOS w/ regress filled days from p84 — found in other code(pier25_gap_fill_15_min.R)
#   7. Barometric pressure — NOAA Battery station (daily mean)
# =============================================================================

#necessary packages
library(dplyr)
library(lubridate)
library(tidyr)
library(dataRetrieval)   # USGS streamgage retrieval
library(lunar)           # lunar cycle
#library(chillR)          # photoperiod / day length //piece of junk that would not work for me
library(readr)  
library(arrow)          # for reading parquet files and directing internally 

#setwd("~/Data/hudson_project/water_quality")
# -----------------------------------------------------------------------------
# STUDY PERIOD SPINE
#   Every covariate will be left-joined onto this so gaps are explicit NAs
#  rather than silently dropped rows.
# -----------------------------------------------------------------------------
# setting start date as 2 weeks prior to first detection
study_start <- as.Date("2022-03-16")
study_end   <- as.Date("2025-10-16")

study_days <- tibble(
  Date = seq(study_start, study_end, by = "day")
)


# =============================================================================
#   WATER TEMPERATURE — NOAA Battery (Station 8518750)
#   Source CSVs: hourly water temp from CO-OPS export
# =============================================================================

temp_bttry_22 <- read.csv("2022_battery_wtr_tmps.csv")
temp_bttry_23 <- read.csv("2023_battery_wtr_tmps.csv")
temp_bttry_24 <- read.csv("2024_battery_wtr_tmps.csv")
temp_bttry_25 <- read.csv("2025_battery_wtr_tmps.csv")

all_temps <- rbind(temp_bttry_22, temp_bttry_23,
                       temp_bttry_24, temp_bttry_25)

all_temps <- all_temps %>%
  mutate(
    Date      = as.Date(Date, format = "%Y/%m/%d"),
    temp_f    = suppressWarnings(as.numeric(Water.Temp...F.)),
    temp_c    = (temp_f - 32) * (5 / 9)
  ) %>%
  filter(!is.na(temp_c))

#getting the mean daily water temperature from the hourly readings 
#since all data will be brought in at a daily scale
daily_temp <- all_temps %>%
  group_by(Date) %>%
  summarise(temp_c = mean(temp_c, na.rm = TRUE), .groups = "drop")

#adding in lag response since change in temperature may be more important 
#than just the raw temperature
daily_temp <- daily_temp %>%
  mutate(
    delta_temp_1d = temp_c - lag(temp_c, n =1),
    delta_temp_3d = temp_c - lag(temp_c, n =3),
    delta_temp_7d = temp_c - lag(temp_c, n =7)
  )


# Getting degree days 
# this has to start prior to our study period (march 2022) for proper cumulative
#calculation that year
degree_days <- daily_temp %>%
  arrange(Date) %>%
  mutate(year = year(Date)) %>%
  group_by(year) %>%
  mutate(
    degree_days = cumsum(temp_c)
  ) %>%
  ungroup() %>%
  select(-year) %>%
  select(
    Date,
    degree_days
  )

head(degree_days)
schema(degree_days)
# =============================================================================
#   RIVER DISCHARGE — USGS Troy (Station 01358000)
#   Parameter 00060 = streamflow (ft³/s), statistic 00003 = daily mean
#
#     LAG RATIONALE:
#       Troy is nearest flow station but far upstream. INput from flow will   
#       be lagged between there and our study area. 
#       Looking at coorelations between discharge and salinity there was 
#       maximum negative coorelation with a lag of one day and second of 2 
#       days so I used 1 day lag.
# =============================================================================

#completed via trailing code but written into parquet for quicker entry and 
# in case USGS stops being operable with current code
daily_discharge <- read_parquet("daily_discharge_troy.parquet")

#LAG_DAYS <- 1L

#hudson_troy_raw <- read_waterdata_daily(
#  monitoring_location_id = "USGS-01358000",
#  parameter_code         = "00060",
#  unit_of_measure        = "ft^3/s",
#  statistic_id           = "00003",
#  time                   = "2022-01-01T00:00:00Z/2026-01-01T00:00:00Z"
#)

# We pull from 2022-01-01 (not study start) so the lag doesn't create NAs
# for the first LAG_DAYS rows of the study period.
#daily_discharge <- hudson_troy_raw %>%
#  mutate(
#    Date          = as.Date(time),
#    discharge_m3s = value * 0.028317          # ft³/s → m³/s
#  ) %>%
#  select(Date, discharge_m3s) %>%
#  arrange(Date) %>%
#  mutate(
#    discharge_lag1_m3s = lag(discharge_m3s, n = LAG_DAYS),
#    delta_discharge_1d = discharge_lag1_m3s - lag(discharge_lag1_m3s, n=1),
#    delta_discharge_3d = discharge_lag1_m3s - lag(discharge_lag1_m3s, n=3),
#    delta_discharge_7d = discharge_lag1_m3s - lag(discharge_lag1_m3s, n=7)
#  )

#daily_dischrg_1 <- daily_discharge %>% 
#  sf::st_drop_geometry() 

#write_parquet(daily_dischrg_1, "daily_discharge_troy.parquet")

# head(daily_discharge, 10)

# =============================================================================
#   TIDAL METRICS — NOAA Battery MHHW High-Low data
#   Source CSVs: CO-OPS High/Low export, verified water level (m)
# =============================================================================

#setting up function to bring out csv data correctly
read_tides <- function(f) read.csv(f, stringsAsFactors = FALSE)

all_tides <- bind_rows(
  read_tides("2022_MHHW_Battery_HL_CO-OPS_8518750_met (3).csv"),
  read_tides("2023_MHHW_Battery_HL_CO-OPS_8518750_met (3).csv"),
  read_tides("2024_MHHW_Battery_HL_CO-OPS_8518750_met (3).csv"),
  read_tides("2025_MHHW_Battery_HL_CO-OPS_8518750_met.csv")
) %>%
  mutate(
    date_time     = as.POSIXct(paste(Date, Time..GMT.),
                               format = "%Y/%m/%d %H:%M", tz = "GMT"),
    verified_tide = suppressWarnings(as.numeric(Verified..m.))
  ) %>%
  filter(!is.na(verified_tide), !is.na(date_time))

#getting the maximum (verified) and minimum tides as well as the flux (difference between them)
daily_tides <- all_tides %>%
  group_by(Date = as.Date(date_time)) %>%
  summarise(
    tide_max   = max(verified_tide,  na.rm = TRUE),
    tide_min   = min(verified_tide,  na.rm = TRUE),
    tide_range = tide_max - tide_min,
    .groups    = "drop"
  )


# =============================================================================
#   LUNAR PHASE — {lunar} package
#   lunar.phase() returns radians (0 = new moon, π = full moon).
#   lunar_illumination: proportion illuminated (0–1), useful for models
#   that want a continuous, biologically meaningful predictor.
#   lunar_cycle: integer 0-based synodic cycle index (for grouping).
# =============================================================================

lunar_data <- study_days %>%
  mutate(
    lunar_phase       = lunar.phase(Date, name = FALSE),          # radians
    lunar_illumination = (1 - cos(lunar_phase)) / 2,             # 0–1
    lunar_cycle       = floor(lunar.phase(Date, name = FALSE) /
                                (2 * pi) * 29.5)                  # index
  )


# =============================================================================
#   PHOTOPERIOD — chillR did not work. Instead used 
#   mathematical model based on the latitude for daily photoperiod
#   Lat from middle of Upper Bay
#   Hours of daylight computed from model related to position and Julian day.
# =============================================================================

#central to our study area
SITE_LAT <- 40.67059

photo_data <- study_days %>%
  mutate(
    jday = yday(Date),
    photoperiod = 24 - (24/pi) * acos(
      pmin(pmax(
        -tan(SITE_LAT * pi/180) *
        tan(23.44 * pi/180 * cos(2*pi*(jday+10)/365)),
        -1), 1)
    )
  ) %>%
  select(Date, jday, photoperiod)


# =============================================================================
#   SALINITY — USGS and HRECOS station Pier 25 (40.72074	-74.0163)
#   Some days are fileld from regression using nearby pier 84 (40.76462	
#   -74.0032). Mean values were derived from 15 minute interval readings 
#   with a minimum of 48 daily readins (50%) required for mean inclusion.
#   Salinity in psu. 
# =============================================================================

daily_salinity_raw <- read_parquet("daily_mean_sal.parquet")
summary(daily_salinity_raw)
#pulling our salinity values of interest
daily_salinity <- daily_salinity_raw %>% select(
  Date,
  sal_max,
  sal_min,
  sal_mean
)

# =============================================================================
#   BAROMETRIC PRESSURE — NOAA Battery station (8518750)
#    product = air_pressure; daily mean (mb)
#   delta_pressure_mb = 1-day change; falling pressure could be a migration cue
# =============================================================================

#similar to flow, folling code was run and its output was saved here
daily_pressure <- read_parquet("daily_pressure.parquet")

#fetch_coops_met <- function(begin_date, end_date, product,
#                            station = "8518750") {
#  base <- "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter"
#  url  <- paste0(
#    base,
#    "?begin_date=", gsub("-", "", begin_date),
#    "&end_date=",   gsub("-", "", end_date),
#    "&station=",    station,
#    "&product=",    product,
#    "&time_zone=GMT&units=metric",
#    "&application=web_services&format=csv"
#  )
#  df <- tryCatch(
#    read_csv(url, show_col_types = FALSE),
#    error = function(e) {
#      message(product, " API call failed: ", e$message); NULL
#    }
#  )
#  df
#}

#make_31day_chunks <- function(start_date, end_date, max_days = 31) {
#  starts <- seq(as.Date(start_date), as.Date(end_date), by = max_days)
#  ends   <- pmin(starts + (max_days - 1), as.Date(end_date))
#  Map(c, format(starts, "%Y-%m-%d"), format(ends, "%Y-%m-%d"))
#}

#date_chunks <- make_31day_chunks("2022-03-01", "2025-10-31")

# -- Barometric pressure -------------------------------------------------------
#pressure_raw <- bind_rows(Filter(Negate(is.null),
#  lapply(date_chunks, \(x) fetch_coops_met(x[1], x[2], "air_pressure"))
#))



# CO-OPS pressure columns: Date Time, Pressure, Quality
#daily_pressure <- pressure_raw %>%
#  mutate(
#    date_time = as.POSIXct(`Date Time`, format = "%Y-%m-%d %H:%M", tz = "GMT"),
#    pressure  = suppressWarnings(as.numeric(Pressure))
#  ) %>%
#  filter(!is.na(pressure)) %>%
#  group_by(Date = as.Date(date_time)) %>%
#  summarise(
#    pressure_mb   = mean(pressure, na.rm = TRUE),
    # 1-day pressure change: negative = falling (potential movement cue)
#    .groups        = "drop"
#  ) %>%
#  arrange(Date) %>%
#  mutate(
#    delta_pressure_1d = pressure_mb - lag(pressure_mb, 1),
#    delta_pressure_3d = pressure_mb - lag(pressure_mb, n=3),
#    delta_pressure_7d = pressure_mb - lag(pressure_mb, n=7)
#  )

#write_parquet(daily_pressure, "daily_pressure.parquet")


# =============================================================================
#   ASSEMBLE MASTER COVARIATE DATA FRAME
#   Left join everything onto the study spine — gaps become NA so you
#   can see exactly what dates are missing
# =============================================================================

covariates <- study_days %>%
  left_join(daily_temp,      by = "Date") %>%
  left_join(degree_days, by = "Date") %>%
  left_join(daily_discharge %>%
              select(Date, discharge_lag1_m3s, delta_discharge_1d, delta_discharge_3d, delta_discharge_7d),
            by = "Date") %>%
  left_join(daily_tides,     by = "Date") %>%
  left_join(lunar_data,      by = "Date") %>%
  left_join(photo_data,      by = "Date") %>%
  left_join(daily_salinity,  by = "Date") %>%
  left_join(daily_pressure,  by = "Date")

# Quick completeness check
cat("\n--- Covariate completeness (% non-NA) ---\n")
covariates %>%
  summarise(across(everything(), \(x) round(mean(!is.na(x)) * 100, 1))) %>%
  tidyr::pivot_longer(everything(), names_to = "variable", values_to = "pct_complete") %>%
  print(n = Inf)

summary(covariates)
#there are some infinite values currenlty listed for salinity so getting those to 
# show as NA instead
covariates.1 <- covariates %>%
  mutate(
    sal_max = ifelse(is.infinite(sal_max), NA, sal_max),
    sal_min = ifelse(is.infinite(sal_min), NA, sal_min)
  )
summary(covariates.1)
covariates.1 %>%
  summarise(across(everything(), \(x) round(mean(!is.na(x)) * 100, 1))) %>%
  tidyr::pivot_longer(everything(), names_to = "variable", values_to = "pct_complete") %>%
  print(n = Inf)
# sal max = 99.8, sal min = 99.8 and sal mean = 99.5, all else 100

#write_parquet(covariates.1, "covariate_list.parquet")


#==============================================================================
# Looking at correlation of the covariates since lots of them are gong to be
# Not sure exactly what this will bring me but that is why we are here

#==============================================================================

covariates <- read_parquet("covariate_list.parquet")
library(corrplot)
num_vars <- covariates %>% select(where(is.numeric))
cor_matrix <- cor(num_vars, use = "pairwise.complete.obs")
corrplot(cor_matrix, method = "pie", type = "upper", tl.cex = 1.4)

#creating exportable 
#png("correlation_plot_0831.png", coor_plot_var, width = 10, height = 8, units = "in", res = 300)
#corrplot(cor_matrix, method = "pie", type = "upper", tl.cex = 1.4)
#dev.off()

#daily_dets <- read_parquet("daily_fish_counts.parquet")
#dets_and_variables <- left_join(daily_dets, covariates, by = "Date")
#write_parquet(dets_and_variables, "det_variables_insert_0831.parquet")


# going to get rid of a few variables to check this again
list(covariates)
limited_covarites <- covariates %>%
  select(
   Date,
   jday,
   temp_c,
   delta_temp_1d,
   delta_temp_3d,
   delta_temp_7d,
   discharge_lag1_m3s,
   delta_discharge_1d,
   delta_discharge_3d,
   delta_discharge_7d,
   tide_range,
   lunar_illumination,
   lunar_cycle,
   photoperiod,
   sal_mean,
   pressure_mb,
   delta_pressure_1d,
   delta_pressure_3d,
   delta_pressure_7d
  )

# checking coorelation for these more limited variables
num_vars.1 <- limited_covarites %>% select(where(is.numeric))
cor_matrix.1 <- cor(num_vars.1, use = "pairwise.complete.obs")
corrplot(cor_matrix.1, method = "number", type = "upper", tl.cex = 1.0)

#write_parquet(limited_covarites, "limited_covariates.parquet")
# salinity and flow are highly correlated but mechanistically serve two
# different purposes in this as salinity may be environmental predictor for 
# finer scale occupancy while flow could be more of a migration cue. Temperature
# and photoperiod are also fairly coorelated (0.61) as well as temperature and 
# discharge (0.51)

data <- read_parquet("limited_covariates.parquet")
summary(data)
schema(data)
################removed variables###################
#tide min and max, salinity min and max, lunar phase


# =============================================================================
#   SCALING FOR BAYESIAN MODELLING (to have comparable effect size across 
#   different predictors)
#   All continuous predictors scaled to mean = 0, SD = 1.
#   Raw columns retained with suffix _raw.
# =============================================================================

# EXCLUDED:
#   - lunar_cycle: discrete synodic index (0-29), not a continuous measurement
#   - jday: excluded from the model matrix entirely and keeping in photoperiod.
#   Kept raw in the data frame for reference/plotting only.
vars_to_scale <- c(
  "temp_c",
  "delta_temp_1d", "delta_temp_3d", "delta_temp_7d",
  "discharge_lag1_m3s",
  "delta_discharge_1d", "delta_discharge_3d", "delta_discharge_7d",
  "tide_range",
  "lunar_illumination",
  "photoperiod",
  "sal_mean",
  "pressure_mb",
  "delta_pressure_1d", "delta_pressure_3d", "delta_pressure_7d"
)


# Grab mean/sd for every variable before scaling so they can be stored
# and reused (e.g., for scaling new data or back-transforming coefficients)
scaling_params <- data %>%
  summarise(across(all_of(vars_to_scale),
                   list(mean = \(x) mean(x, na.rm = TRUE),
                        sd   = \(x) sd(x,   na.rm = TRUE)),
                   .names = "{.col}__{.fn}")) %>%
  pivot_longer(everything(),
               names_to = c("variable", "stat"),
               names_sep = "__",
               values_to = "value") %>%
  pivot_wider(names_from = stat, values_from = value)

print(scaling_params, n = Inf)


# Create scaled (_z) versions, keeping originals untouched
data_scaled <- data %>%
  mutate(across(all_of(vars_to_scale),
                \(x) as.numeric(scale(x)),
                .names = "{.col}_z"))

# Sanity check: scaled vars should have mean ~0, sd ~1 (ignoring NAs)
data_scaled %>%
  summarise(across(ends_with("_z"),
                   list(mean = \(x) round(mean(x, na.rm = TRUE), 3),
                        sd   = \(x) round(sd(x,   na.rm = TRUE), 3)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  print(n = Inf)

summary(data_scaled)
schema(data_scaled)

# Save both the scaled data and the scaling parameters (needed if you
# ever scale new detection-period data using these same means/sds)
write_parquet(data_scaled, "covariates_scaled.parquet")
write_parquet(scaling_params, "scaling_params.parquet")

#all information kept on my OneDrive with back up on external hard drive
