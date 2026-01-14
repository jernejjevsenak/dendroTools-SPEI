# ----------------------------
# 1) Load required libraries
# ----------------------------
library(dplyr)
library(lubridate)
library(reshape2)    # note: largely superseded by tidyr
library(dendroTools)
library(SPEI)
library(zoo)
library(dplR)
library(lmomco)
library(ggplot2)

# -----------------------------------------
# 2) User inputs / key analysis parameters
# -----------------------------------------
lat <- 51.0386          # Site latitude in decimal degrees (positive = Northern Hemisphere)
lower_limit <- 21       # Minimum aggregation window (days) for water deficit in SPEI
upper_limit <- 270      # Maximum aggregation window (days) for water deficit in SPEI
previous_year <- TRUE   # If TRUE, include previous-year climate in response analyses

# -----------------------------------------
# 3) Read dendro data (TRWi chronology)
# -----------------------------------------
temp_rwl <- read.crn("my_chron.crn")  # dplR chronology format

# Remove sample depth column if present (not needed for later steps)
if ("samp.depth" %in% names(temp_rwl)) {
  temp_rwl$samp.depth <- NULL
}

# -----------------------------------------
# 4) Read daily/monthly climate inputs
#    Expectation: these CSVs contain a time column (or year/month/day)
#    plus a value column (e.g., Tmean, Tmax, Tmin, Psum).
# -----------------------------------------
Tavg <- read.csv("my_chron_Tmean.csv")
Tmax <- read.csv("my_chron_Tmax.csv")
Tmin <- read.csv("my_chron_Tmin.csv")
Prec <- read.csv("my_chron_Psum.csv")

# -----------------------------------------------------------
# 5) Helper function: Daily Potential Evapotranspiration (PET)
#    Method: Hargreaves (commonly used when only temperature
#    is available). Radiation terms follow SPEI package logic.
#
#    Inputs:
#      date : Date vector
#      tavg : mean daily temperature (°C)
#      tdif : daily temperature range (Tmax - Tmin) (°C), must be >= 0
#      lat  : latitude in decimal degrees
#
#    Output:
#      PET in mm/day (non-negative)
# -----------------------------------------------------------
HargreavesPET_daily <- function(date, tavg, tdif, lat) {
  
  # Basic input checks (helps catch silent problems early)
  stopifnot(inherits(date, "Date"))
  stopifnot(is.numeric(tavg), is.numeric(tdif), is.numeric(lat))
  stopifnot(lat >= -90 && lat <= 90)
  
  # Day of year (1..365/366 depending on lubridate’s yday)
  doy <- yday(date)
  
  # Convert latitude to radians
  phi <- pi / 180 * lat
  
  # Solar declination (radians) and inverse relative distance Earth-Sun
  # (Using 365 in both terms for consistency; adjust if you handle leap years explicitly)
  delta <- 0.409 * sin(2 * pi / 365 * doy - 1.39)
  dr    <- 1 + 0.033 * cos(2 * pi / 365 * doy)
  
  # Sunset hour angle (ws). Numerical issues can yield NaN for extreme latitudes/doys.
  ws <- suppressWarnings(acos(-tan(phi) * tan(delta)))
  ws[is.na(ws) | is.nan(ws)] <- 0.1  # small fallback to avoid crashing downstream
  
  # Extraterrestrial radiation (Ra) in MJ m^-2 day^-1
  Ra <- (24 * 60 / pi) * 0.0820 * dr * (
    (ws * sin(phi) * sin(delta)) +
      (cos(phi) * cos(delta) * sin(ws))
  )
  
  # Hargreaves PET (mm/day). tdif must be non-negative; clamp just in case.
  tdif <- pmax(tdif, 0)
  
  PET <- 0.0023 * (tavg + 17.8) * sqrt(tdif) * Ra
  
  # PET cannot be negative
  PET <- pmax(PET, 0)
  
  return(PET)
}

# ------------------------------------------------------------
# SPEI_daily(): Daily SPEI from a daily water deficit series (wd)
#
# - Aggregates wd using a kernel over 'scale' days (e.g., 21, 30, 90...)
# - Fits a distribution per seasonal group (here: month-of-year)
# - Transforms to standardized normal (SPEI-like index)
#
# Credits: logic is based on SPEI package approach (distribution fit + CDF->qnorm).
# ------------------------------------------------------------
SPEI_daily <- function(wd_data, scale = 21, kernel = list(type = "rectangular", shift = 0),
                       distribution = "log-Logistic", fit = "ub-pwm", na.rm = FALSE,
                       ref.start = NULL, ref.end = NULL, x = FALSE, params = NULL) {
  
  # -------------------------
  # 1) Basic input checks
  # -------------------------
  if (!is.data.frame(wd_data)) stop("wd_data must be a data.frame.")
  if (!("wd" %in% names(wd_data))) stop("wd_data must contain a column named 'wd'.")
  
  scale <- as.numeric(scale)
  na.rm <- as.logical(na.rm)
  x <- as.logical(x)
  
  if (!(distribution %in% c("log-Logistic", "Gamma", "PearsonIII"))) {
    stop('distribution must be one of "log-Logistic", "Gamma" or "PearsonIII"')
  }
  if (!(fit %in% c("max-lik", "ub-pwm", "pp-pwm"))) {
    stop('fit must be one of "ub-pwm" (default), "pp-pwm" or "max-lik"')
  }
  if ((!is.null(ref.start) && length(ref.start) != 2) || (!is.null(ref.end) && length(ref.end) != 2)) {
    stop("ref.start and ref.end must be numeric vectors of length 2: c(year, month).")
  }
  
  # -------------------------
  # 2) Decide how to define seasons (groups)
  # -------------------------
  has_date <- "date" %in% names(wd_data)
  
  if (has_date) {
    # Preferred path for daily data: month-of-year grouping (1..12)
    wd_data$date <- as.Date(wd_data$date)
    if (anyNA(wd_data$date)) stop("wd_data$date could not be converted with as.Date().")
    
    wd_data <- wd_data[order(wd_data$date), , drop = FALSE]
    
    water_deficit <- wd_data$wd
    dates <- wd_data$date
    
    cyc_all <- as.integer(format(dates, "%m"))  # seasonal group for every observation
    fr <- 12                                     # number of seasonal groups
  } else {
    # Fallback: keep behavior for no-date input using a simple 12-cycle
    water_deficit <- wd_data[, "wd"]
    tmp_ts <- ts(as.matrix(water_deficit), frequency = 12)
    cyc_all <- cycle(tmp_ts)
    fr <- frequency(tmp_ts)
  }
  
  # Represent the series as a (n x m) matrix-like object
  # NOTE: if there's only one series, m = 1
  if (has_date) {
    df_water_deficit <- ts(as.matrix(water_deficit))          # frequency not needed when using dates
  } else {
    df_water_deficit <- ts(as.matrix(water_deficit), frequency = 12)  # ensure cycle() works
  }
  
  m <- ncol(df_water_deficit)
  
  # -------------------------
  # 3) Prepare coefficient container
  # -------------------------
  coef <- switch(distribution,
                 "Gamma" = array(NA, c(2, m, fr), list(par = c("alpha", "beta"), colnames(df_water_deficit), NULL)),
                 "log-Logistic" = array(NA, c(3, m, fr), list(par = c("xi", "alpha", "kappa"), colnames(df_water_deficit), NULL)),
                 "PearsonIII" = array(NA, c(3, m, fr), list(par = c("mu", "sigma", "gamma"), colnames(df_water_deficit), NULL))
  )
  
  dim_one <- ifelse(distribution == "Gamma", 2, 3)
  
  if (!is.null(params)) {
    if (dim(params)[1] != dim_one || dim(params)[2] != m || dim(params)[3] != fr) {
      stop(sprintf("params must have dimensions (%d, %d, %d)", dim_one, m, fr))
    }
  }
  
  # -------------------------
  # 4) Reference period handling (optional)
  #    Interprets ref.start/ref.end as c(year, month), like SPEI monthly convention.
  # -------------------------
  if (has_date && !is.null(ref.start) && !is.null(ref.end)) {
    ref_start_date <- as.Date(sprintf("%04d-%02d-01", ref.start[1], ref.start[2]))
    ref_end_date   <- as.Date(sprintf("%04d-%02d-01", ref.end[1], ref.end[2]))
    ref_end_date   <- as.Date(format(ref_end_date + 31, "%Y-%m-01")) - 1
    
    fit_idx <- which(dates >= ref_start_date & dates <= ref_end_date)
  } else {
    fit_idx <- seq_along(water_deficit)
  }
  
  # Output container for standardized values (same shape as input ts)
  std <- df_water_deficit * NA
  
  # -------------------------
  # 5) Main loop over series (columns)
  # -------------------------
  for (s in 1:m) {
    
    acu <- as.numeric(df_water_deficit[, s])
    
    # -------------------------
    # 5a) Aggregate (kernel sum) over 'scale' days
    # -------------------------
    if (scale > 1) {
      wgt <- kern(scale, kernel$type, kernel$shift)
      
      # Fit subset aggregation
      acu_fit <- acu[fit_idx]
      if (length(acu_fit) >= scale) {
        tmp <- rowSums(embed(acu_fit, scale) * wgt, na.rm = na.rm)
        acu_fit_agg <- rep(NA, length(acu_fit))
        acu_fit_agg[scale:length(acu_fit)] <- tmp
      } else {
        acu_fit_agg <- rep(NA, length(acu_fit))
      }
      
      # Full-series aggregation
      if (length(acu) >= scale) {
        tmp2 <- rowSums(embed(acu, scale) * wgt, na.rm = na.rm)
        acu_agg <- rep(NA, length(acu))
        acu_agg[scale:length(acu)] <- tmp2
      } else {
        acu_agg <- rep(NA, length(acu))
      }
      
    } else {
      acu_fit_agg <- acu[fit_idx]
      acu_agg <- acu
    }
    
    # Build seasonal grouping vectors for fit and prediction portions
    if (has_date) {
      cyc_fit  <- cyc_all[fit_idx]
      cyc_pred <- cyc_all
    } else {
      cyc_pred <- cyc_all
      cyc_fit  <- cyc_all[fit_idx]
    }
    
    # -------------------------
    # 5b) Fit per seasonal group (month-of-year)
    # -------------------------
    for (c in 1:fr) {
      
      f_local <- which(cyc_fit == c)   # indices within fit subset
      ff      <- which(cyc_pred == c)  # indices within full series
      
      fit_vals <- acu_fit_agg[f_local]
      fit_vals <- fit_vals[!is.na(fit_vals)]
      
      if (length(fit_vals) == 0) {
        std[ff, s] <- NA
        coef[, s, c] <- NA
        next()
      }
      
      # (Kept from your original) Sorting doesn't change estimation but is harmless
      month <- sort.default(fit_vals, method = "quick")
      
      # Guard against degenerate inputs
      month_sd <- sd(month, na.rm = TRUE)
      if (is.na(month_sd) || month_sd == 0) {
        std[ff, s] <- NA
        coef[, s, c] <- NA
        next()
      }
      
      # Zero handling for Gamma/PearsonIII
      pze <- 0
      if (distribution != "log-Logistic") {
        pze <- sum(month == 0) / length(month)
        month <- month[month > 0]
      }
      
      if (length(month) < 4) {
        std[ff, s] <- NA
        coef[, s, c] <- NA
        next()
      }
      
      # -------------------------
      # 5c) Parameter estimation (or use provided params)
      # -------------------------
      if (is.null(params)) {
        
        pwm <- switch(fit,
                      "pp-pwm"  = pwm.pp(month, -0.35, 0, nmom = 3),
                      "ub-pwm"  = TLMoments::PWM(month, order = 0:2),
                      "max-lik" = TLMoments::PWM(month, order = 0:2)  # used as starting point
        )
        
        lmom_obj <- pwm2lmom(pwm)
        if (!are.lmom.valid(lmom_obj) || anyNA(lmom_obj$lambdas) || any(is.nan(lmom_obj$lambdas))) {
          std[ff, s] <- NA
          coef[, s, c] <- NA
          next()
        }
        
        fortran_vec <- c(lmom_obj$lambdas[1:2], lmom_obj$ratios[3])
        
        f_params <- switch(distribution,
                           "log-Logistic" = tryCatch(lmom::pelglo(fortran_vec), error = function(e) lmom::parglo(lmom_obj)$para),
                           "Gamma"        = tryCatch(lmom::pelgam(fortran_vec), error = function(e) lmom::pargam(lmom_obj)$para),
                           "PearsonIII"   = tryCatch(lmom::pelpe3(fortran_vec), error = function(e) lmom::parpe3(lmom_obj)$para)
        )
        
        # Optional ML refinement for log-Logistic
        if (distribution == "log-Logistic" && fit == "max-lik") {
          f_params <- parglo.maxlik(month, f_params)$para
        }
        
      } else {
        f_params <- as.vector(params[, s, c])
      }
      
      # -------------------------
      # 5d) CDF -> qnorm transform for this seasonal group
      # -------------------------
      pred_vals <- acu_agg[ff]
      
      cdf_res <- switch(distribution,
                        "log-Logistic" = lmom::cdfglo(pred_vals, f_params),
                        "Gamma"        = lmom::cdfgam(pred_vals, f_params),
                        "PearsonIII"   = lmom::cdfpe3(pred_vals, f_params)
      )
      
      # Avoid exact 0/1 (infinite qnorm)
      cdf_res <- pmin(pmax(cdf_res, 1e-6), 1 - 1e-6)
      
      std[ff, s] <- qnorm(cdf_res)
      coef[, s, c] <- f_params
      
      # Zero adjustment for Gamma/PearsonIII
      if (distribution != "log-Logistic") {
        std[ff, s] <- qnorm(pze + (1 - pze) * pnorm(std[ff, s]))
      }
    }
  }
  
  colnames(std) <- colnames(df_water_deficit)
  
  # -------------------------
  # 6) Build return object (same structure you had)
  # -------------------------
  z <- list(
    call = match.call(expand.dots = FALSE),
    fitted = std,
    coefficients = coef,
    scale = scale,
    kernel = list(
      type = kernel$type,
      shift = kernel$shift,
      values = kern(scale, kernel$type, kernel$shift)
    ),
    distribution = distribution,
    fit = fit,
    na.action = na.rm
  )
  
  if (x) z$df_water_deficit <- df_water_deficit
  if (!is.null(ref.start)) z$ref.period <- rbind(ref.start, ref.end)
  
  wd_data$SPEI <- as.numeric(z$fitted)
  return(wd_data)
}


################################################################################
# 1) Prepare daily climate inputs (Tavg/Tmax/Tmin/Prec must be aligned row-wise)
################################################################################

# --- Build a daily Date vector from Y/M/D columns (assumes Tavg has these columns)
dates <- as.Date(sprintf("%04d-%02d-%02d", Tavg[, "Y"], Tavg[, "M"], Tavg[, "D"]))

# --- Daily temperature range (Tmax - Tmin); clamp negatives to 0
T_diff <- as.numeric(Tmax[, "Tmax"]) - as.numeric(Tmin[, "Tmin"])
T_diff <- pmax(T_diff, 0)  # if Tmin > Tmax (bad/odd input), use 0

# --- Calculate PET (mm/day) using Hargreaves
PET <- HargreavesPET_daily(
  date = dates,
  tavg = as.numeric(Tavg[, "Tmean"]),
  tdif = T_diff,
  lat  = lat
)

# --- Water deficit (WD): precipitation - PET (mm/day)
WD_values <- as.numeric(Prec[, "Prec"]) - PET

# --- Build the WD data frame for SPEI_daily()
# Keep same structure as original: reusing Tavg’s Y/M/D fields
WD <- Tavg
WD$Tmean <- NULL
WD$wd <- WD_values

################################################################################
# 2) For each time scale (days), compute SPEI and correlate with tree-ring series
################################################################################

# List to store correlation “heatmap rows” for each scale length
temporal_matrix_list <- list()
place_holder <- 1

for (ij in lower_limit:upper_limit) {
  
  # Keep the original chronology subset behavior
  temp_rwl_subset <- temp_rwl
  
  # --- Compute daily SPEI aggregated over 'ij' days
  SPEI_temp <- SPEI_daily(wd_data = WD, scale = ij)
  
  # --- IMPORTANT CORRECTION (kept): shift SPEI so it aligns with end-of-window
  # (This preserves your existing alignment convention.)
  SPEI_temp$SPEI <- c(
    SPEI_temp$SPEI[ij:length(SPEI_temp$SPEI)],
    rep(NA, ij - 1)
  )
  
  # --- Prepare daily SPEI in wide "year x day-of-year" matrix format
  SPEI_temp$date <- paste0(SPEI_temp$Y, "-", SPEI_temp$M, "-", SPEI_temp$D)
  SPEI_temp <- data_transform(SPEI_temp[, c("date", "SPEI")])  # wide df: rows=years, cols=DOY
  
  # --- Fix missing DOY 366 using DOY 365 (kept)
  # NOTE: This assumes your wide format has at least 366 columns.
  SPEI_temp[, 366] <- ifelse(is.na(SPEI_temp[, 366]), SPEI_temp[, 365], SPEI_temp[, 366])
  
  # --- If previous year should be considered, concatenate [prev-year DOY cols] + [current-year DOY cols]
  if (previous_year == TRUE) {
    
    SPEI_temp$temp_year <- row.names(SPEI_temp)
    SPEI_temp <- dplyr::arrange(SPEI_temp, desc(temp_year))
    SPEI_temp <- years_to_rownames(SPEI_temp, "temp_year")
    
    SPEI_temp_previous <- SPEI_temp[-1, , drop = FALSE]
    SPEI_temp_current  <- SPEI_temp[-nrow(SPEI_temp), , drop = FALSE]
    
    row_names_current <- row.names(SPEI_temp_current)
    
    SPEI_temp <- cbind(SPEI_temp_previous, SPEI_temp_current)
    SPEI_temp <- data.frame(SPEI_temp)
    row.names(SPEI_temp) <- row_names_current
  }
  
  # --- Merge tree-ring chronology (temp_rwl_subset) with SPEI matrix by year
  ncol_temp_rwl_subset <- ncol(temp_rwl_subset)
  colnames_temp_rwl_subset <- colnames(temp_rwl_subset)
  
  SPEI_temp$temp_year <- row.names(SPEI_temp)
  temp_rwl_subset$temp_year <- row.names(temp_rwl_subset)
  
  temporal_data <- merge(temp_rwl_subset, SPEI_temp, by = "temp_year")
  
  # Split merged object back into TRW and SPEI parts (kept)
  temp_rwl_subset <- data.frame(
    temporal_data[, c(2:(1 + ncol_temp_rwl_subset))],
    row.names = temporal_data$temp_year
  )
  colnames(temp_rwl_subset) <- colnames_temp_rwl_subset
  
  SPEI_temp <- data.frame(
    temporal_data[, c((1 + ncol_temp_rwl_subset + 1):ncol(temporal_data))],
    row.names = temporal_data$temp_year
  )
  
  # --- Correlation storage: 1 row per scale; columns represent DOY (and prev-year DOY if enabled)
  temporal_matrix <- matrix(NA, nrow = 1, ncol = ncol(SPEI_temp))
  
  # We only compute correlations where a full ij-day window can end at that column index.
  # Your original indexing sets results at column (j + ij).
  max_j <- ncol(SPEI_temp) - ij
  
  for (j in 0:max_j) {
    
    # Column representing the window end day (j+ij), per your convention
    x <- SPEI_temp[, (j + 1)]
    x <- matrix(x, nrow = nrow(SPEI_temp), ncol = 1)
    
    temporal_correlation <- cor(
      temp_rwl_subset[, 1], x[, 1],
      method = "pearson",
      use = "pairwise.complete.obs"
    )
    
    temporal_matrix[1, j + ij] <- temporal_correlation
  }
  
  # Save this scale’s correlation row into the list
  temporal_matrix_list[[place_holder]] <- temporal_matrix
  place_holder <- place_holder + 1
  
  print(ij)
}

# Combine list into a single data frame: rows = scales, cols = DOY (and prev-year DOY if used)
daily_SPEI_correlations <- data.frame(do.call(rbind, temporal_matrix_list))

################################################################################
# 3) Visualization example (ggplot2 heatmap)
################################################################################

daily_SPEI_correlations$season_length <- seq(lower_limit, upper_limit)
melted <- melt(daily_SPEI_correlations, id.vars = c("season_length"))

# Hide small correlations (threshold kept)
melted$value <- ifelse(abs(melted$value) < 0.25, NA, melted$value)

ggplot(melted, aes_(x = ~as.numeric(variable), y = ~season_length, fill = ~value)) +
  geom_tile() +
  xlab("Day of Year") +
  ylab("Season Length") +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = c(20, 150, 300, 450, 584, 725),
    labels = c("jan*", "may*", "oct*", "MAR", "AUG", "NOV")
  ) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(30, 270, by = 30)) +
  geom_vline(xintercept = 366, linetype = "dashed") +
  scale_fill_gradient2(
    low = "red", mid = "white", high = "blue",
    na.value = "gray97",
    midpoint = 0, limits = c(-0.8, 0.8)
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title.y = element_text(size = 18),
    text = element_text(size = 18),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 16),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(3, "line"),
    panel.background = element_rect(
      fill = "gray97",
      colour = "gray80",
      size = 0.5,
      linetype = "solid"
    )
  )

ggsave("SPEI_example3.png", height = 7, width = 10)
