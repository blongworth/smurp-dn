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

#' Read GEMS data
#' 
#' @param filename Path of a file in GEMS format
#' 
#' @return A dataframe of GEMS data in long format
read_gems <- function(filename) {
  raw_file <- read_lines(filename)
  gems_raw <- raw_file[str_starts(raw_file, "(R:)*\\d+,")]
  gems_data <- str_remove(gems_raw, "^R:")
  read_csv(I(gems_data),
           col_names = c("hour", "min", "sec", "month", 
                         "day", "year", "mass", "current"),
           col_types = list(
             hour = col_integer(),
             min = col_integer(),
             sec = col_integer(),
             month = col_integer(),
             day = col_integer(),
             year = col_integer(),
             mass = col_integer(),
             current = col_double())) |> 
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

#' Pivot RGA data to wide format with a timestamp for each cycle
#' 
#' @param df A dataframe of RGA data in long format
#' 
#' @return A dataframe of RGA data in wide format
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

#' Calculate mass ratios to normalize to mass 28
#' 
#' @param df A dataframe of RGA data in wide format
#' 
#' @return A dataframe of with additional columns for mass ratios
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
  
  coef(lm(sdf[[2]] ~ sdf[[1]]))[2] * 3600
}

#' Interactive plot of RGA data
#' 
#' @param df_wide A dataframe of RGA data in wide format
#' @param title Title of the plot
#' 
#' @return A dygraph object
plot_wide <- function(df_wide, title = "GEMS Data") {
  df_wide %>% 
    select(timestamp, mass_28, mass_29, mass_30, mass_32, mass_40) %>% 
    dygraph(main = title) %>% 
    dyAxis("y", label = "Pressure (Torr)") %>%
    dyOptions(logscale = TRUE) |> 
    dyRangeSelector() |>
    dyLegend()
}

#' plot of RGA data
#' 
#' @param df_wide A dataframe of RGA data in wide format
#' @param title Title of the plot
#' 
#' @return A dygraph object
plot_rga <- function(df_wide, title = "GEMS Data") {
  df_wide %>% 
    select(timestamp, mass_28, mass_29, mass_30, mass_32, mass_40) %>% 
    pivot_longer(cols = starts_with("mass"), names_to = "mass", names_prefix = "mass_", values_to = "pressure") |> 
    ggplot(aes(timestamp, pressure, color = mass)) +
      geom_line() +
      labs(title = title,
           x = "",
           y = "Partial pressure [log Pa]") +
      scale_y_log10()
}

#' plot of RGA data
#' 
#' @param df_wide A dataframe of RGA data in wide format
#' @param title Title of the plot
#' 
#' @return A dygraph object
plot_rga_conc <- function(df_wide, title = "GEMS Data") {
  df_wide %>% 
    select(et, umol_29, umol_30, umol_40) %>% 
    pivot_longer(cols = starts_with("umol"), names_to = "species", names_prefix = "umol_", values_to = "conc") |> 
    mutate(Species = case_when(species == 29 ~ "14/15N2",
                               species == 30 ~ "15/15N2",
                               species == 40 ~ "Argon")) |>
    ggplot(aes(et/3600, conc, color = Species)) +
      geom_line() +
      labs(title = title,
           x = "Elapsed time [h]",
           y = "Concentration [log umol]") +
      scale_y_log10()
}