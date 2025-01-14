# Functions for reading and plotting GEMS data

#' Calculate nitrogen saturation using GSW functions
#'
#' Calculate nitrogen solubility (ml/L) using Hamme and Emerson (2004) equations.
#'
#' @param sal Salinity (PSU)
#' @param temp Temperature (°C)
#'
#' @return Nitrogen saturation (μmol/kg)
#'
#' @examples
#' 
#' sat <- cal_n2_sat(35, 10)
#' sat
#' all.equal(sat, 500.885)
#' 
#' @export
cal_n2_sat <- function(sal, temp) {
  # scaled temp
  ts <- log((298.15 - temp)/(273.15 + temp))
  
  # Constants from Hamme and Emerson (2004)
  A0 <- 6.42931
  A1 <- 2.92704
  A2 <- 4.32531
  A3 <- 4.69149
  B0 <- -7.44129e-3
  B1 <- -8.02566e-3
  B2 <- -1.46775e-2
  
  # Calculate N2 solubility in umol/kg
  lnC <- A0 + A1 * ts + A2 * ts^2 + A3 * ts^3 + sal * (B0 + B1 * ts + B2 * ts^2)
  exp(lnC)
}

#' Calculate argon saturation using GSW functions
#'
#' This function calculates argon saturation in μmol/kg
#' using the equation from Hamme and Emerson (2004) 
#' and an empirical formula for Ar solubility.
#'
#' @param sal Salinity (PSU)
#' @param temp Temperature (°C)
#'
#' @return Argon saturation (μmol/kg)
#'
#' @examples
#' 
#' sat <- cal_ar_sat(35, 10)
#' sat
#' all.equal(sat, 13.4622)
#' 
#' @export
# all.equal(cal_ar_sat(35, 10), 13.4622)
cal_ar_sat <- function(sal, temp) {
  # Calculate argon solubility (umol/kg) using 
  # scaled temp
  ts <- log((298.15 - temp)/(273.15 + temp))
  
  # Constants from Hamme and Emerson (2004)
  A0 <- 2.79150
  A1 <- 3.17609
  A2 <- 4.13116
  A3 <- 4.90379
  B0 <- -6.96233e-3
  B1 <- -7.66670e-3
  B2 <- -1.16888e-2
  
  # Calculate Ar solubility in umol/kg
  lnC <- A0 + A1 * ts + A2 * ts^2 + A3 * ts^3 + sal * (B0 + B1 * ts + B2 * ts^2)
  exp(lnC)
}

process_gems <- function(df, bg_et_range, nit_sat_umol) {
  # get mean timestamps and widen data
  df_wide <- rga_wider(df)
  # determine backgrounds
  df_bg <- df_wide |> 
    filter(et > bg_et_range[1],
           et < bg_et_range[2])
  
  bg_29 <- df_bg |> 
    pull(mass_29) |> 
    mean()
  
  bg_30 <- df_bg |> 
    pull(mass_30) |> 
    mean()
  
  # add ratios and concentrations
  norm_rga(df_wide, bg_29, bg_30, nit_sat_umol)
  
}

#' Calculate rate based on slope of linear fit
#'
#' @param df a dataframe with 2 columns: elapsed time (sec) and the measurement of interest
#' @param et_range 2 element vector with the range of elapsed times to calculate rate over
#'
#' @return A rate in units of measurement per day
#' @export
calc_rate <- function(df, et_center){
  sdf <- df %>% 
    select(et, umol_30) %>% 
    filter(et > et_center - 1000,
           et < et_center + 1000)
  
  coef(lm(sdf[[2]] ~ sdf[[1]]))[2] * 3600 * 24
}

rga_wider <- function(df) {
  df %>% 
    mutate(cycle = cumsum(mass == 18)) %>% 
    group_by(cycle) %>% 
    mutate(cycle_ts = mean(timestamp)) %>% 
    ungroup() %>% 
    select(timestamp = cycle_ts, mass, pressure) %>% 
    pivot_wider(names_from = mass, names_prefix = "mass_",
                values_from = pressure)
}

rga_ratios <- function(df) {
  df %>% 
    mutate(mass_28_18 = mass_28 / mass_18,
           mass_29_18 = mass_29 / mass_18,
           mass_30_18 = mass_30 / mass_18,
           mass_28_40 = mass_28 / mass_40,
           mass_29_40 = mass_29 / mass_40,
           mass_30_40 = mass_30 / mass_40,
           mass_40_28 = mass_40 / mass_28)
}

# generate roxygen documentation for this function
#' Normalize RGA data
#' 
#' Calculates the molar concentration of N2 and Ar.
#' Background pressures are subtracted from masses 29 and 30,
#' assuming concentration of 29 and 30 are 0 at the start of the measurement.
#' 
#' Ar background is assumed to be saturated at the start of the measurement.
#' 
#' @param df A dataframe of RGA data in long format.
#' @param bg_29 Background pressure for mass 29
#' @param bg_30 Background pressure for mass 30
#' @param bg_40_28 Background pressure for mass 40/28
#' @param nit_sat_umol Nitrogen saturation in μmol/kg
#' @param ar_sat_umol Argon saturation in μmol/kg
#' @param t0 Timestamp of the first measurement
#' 
#' @return A dataframe of RGA data in long format with additional columns for molar
#' concentration of N2 and Ar, and ratios of the masses to mass 28 and mass 40.
norm_rga <- function(df, 
                     bg_29, # raw background
                     bg_30, 
                     bg_40_28,
                     nit_sat_umol, 
                     ar_sat_umol,
                     t0 = df$timestamp[1]) {
  df %>% 
    mutate(et = as.integer(timestamp - t0),
           # assumes no 29 and 30 before label
           mass_29_28 = ( mass_29 - bg_29 ) / mass_28, 
           mass_30_28 = ( mass_30 - bg_30 ) / mass_28,
           # molar concentration based on nitrogen saturation
           umol_29 = mass_29_28 * nit_sat_umol,
           umol_30 = mass_30_28 * nit_sat_umol,
           umol_40 = mass_40_28 * (ar_sat_umol/bg_40_28))
}

read_gems <- function(file) {
  raw_file <- read_lines(file)
  gems_raw <- raw_file[str_starts(raw_file, "(R:)*\\d+,")]
  gems_data <- str_remove(gems_raw, "^R:")
  read_csv(I(gems_data),
           col_names = c("hour", "min", "sec", "month", 
                         "day", "year", "mass", "current")) |> 
    mutate(mass = as.factor(mass),
           current = current*1E-16,
           pressure = current/0.0801,
           timestamp = lubridate::make_datetime(year, month, day, hour, min, sec, 
                                                tz = "America/New_York")) |> 
    # hack for bad first file
    mutate(timestamp = ifelse(timestamp > "2024-08-18",
                              timestamp,
                              timestamp + 177775450),
         timestamp = as.POSIXct(timestamp)) |> 
    select(timestamp, mass, current, pressure)
}

plot_wide <- function(df_wide) {
  df_wide %>% 
    select(timestamp, mass_28, mass_29, mass_30, mass_32, mass_40) %>% 
    dygraph() %>% 
    dyOptions(logscale = TRUE)
}

plot_gems <- function(data, log = TRUE) {
  stopifnot(is.logical(log))
  p <- data |> 
    ggplot(aes(timestamp, pressure, color = mass)) +
    geom_line() +
    ylab("Pressure (Torr)") +
    theme(axis.title.x = element_blank())
  
  if (log) {
    p +
      scale_y_log10()
  } else {
    p
  }
}

#' Read SRS RGASoft files
#' 
#' Old RGA software tabular format. Watch out for missing spaces causing NA's
#'
#' @param filename File to read.
#'
#' @return A dataframe of RGA data in long format.
#' @export
#'
read_rgasoft <- function(filename) {
  
  # get timestamp
  ts <- read_lines(filename, n_max = 1)[[1]] |> 
    lubridate::mdy_hms()
  #masses <- read_table(file)
  # rowlist <- read_lines(filename, skip = 30, n_max = 1) |> 
  #   str_split("\\s+") 
  # column_names <- rowlist[[1]]
  #read_lines(filename, skip = 32) |> 
  #  str_replace("(\\d)(-)", "$1 $2") |> 
  read_csv(filename,
           col_names = c("time_s", "mass_18", "mass_28", "mass_29", 
                         "mass_30", "mass_32", "mass_40", "mass_44", 
                         "mass_45", "mass_46", "X"), 
           col_types = cols(
             time_s = col_number(),
             mass_18 = col_double(),
             mass_28 = col_double(),
             mass_29 = col_double(),
             mass_30 = col_double(),
             mass_32 = col_double(),
             mass_40 = col_double(),
             mass_44 = col_double(),
             mass_45 = col_double(),
             mass_46 = col_double(),
             X = col_logical()
           ), skip = 32) |> 
    select(-X) |> 
    mutate(timestamp = ts + time_s) |> 
    select(!time_s) |> 
    pivot_longer(cols = starts_with("mass"),
                 names_to = "mass",
                 names_pattern = "mass_(.+)",
                 values_to = "pressure") |> 
    mutate(mass = as.factor(mass))
}

#' Read fixed width RGASoft files
#'
#' @param filename 
#'
#' @return
#' @export
#'
#' @examples
read_rgasoft_ascii <- function(filename) {
  
  # get timestamp
  ts <- read_lines(filename, n_max = 1)[[1]] |> 
    mdy_hms()
  #masses <- read_table(file)
  # rowlist <- read_lines(filename, skip = 30, n_max = 1) |> 
  #   str_split("\\s+") 
  # column_names <- rowlist[[1]]
  #read_lines(filename, skip = 32) |> 
  #  str_replace("(\\d)(-)", "$1 $2") |> 
  read_table(filename,
             col_names = c("time_s", "mass_18", "mass_28", "mass_29", "mass_30", "mass_32", "mass_40", "mass_44", "mass_45", "mass_46", "X"), 
             col_types = cols(
  time_s = col_number(),
  mass_18 = col_double(),
  mass_28 = col_double(),
  mass_29 = col_double(),
  mass_30 = col_double(),
  mass_32 = col_double(),
  mass_40 = col_double(),
  mass_44 = col_double(),
  mass_45 = col_double(),
  mass_46 = col_double(),
  X = col_logical()
), skip = 32) |> 
    select(-X) |> 
    mutate(timestamp = ts + time_s)
}
