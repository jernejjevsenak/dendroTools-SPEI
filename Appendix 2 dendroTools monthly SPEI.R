# ----------------------------
# 1) Load required libraries
# ----------------------------
library(dplyr)
library(lubridate)
library(reshape2)
library(dendroTools)
library(SPEI)
library(zoo)
library(dplR)
library(lmomco)
library(ggplot2)

################################################################################
# 2) Inputs and data loading
################################################################################

lat <- 51.0386          # site latitude (decimal degrees)
lower_limit <- 1        # minimum aggregation window (months)
upper_limit <- 12       # maximum aggregation window (months)
previous_year <- FALSE  # if TRUE, include previous-year months in predictors

# --- Tree-ring chronology
temp_rwl <- read.crn("my_chron.crn")
if ("samp.depth" %in% names(temp_rwl)) temp_rwl$samp.depth <- NULL

# --- Daily climate data (assumes columns include Y, M and variable)
Tavg <- read.csv("my_chron_Tmean.csv")
Tmax <- read.csv("my_chron_Tmax.csv")
Tmin <- read.csv("my_chron_Tmin.csv")
Prec <- read.csv("my_chron_Psum.csv")

# --- Convert daily to monthly (mean temp, sum precip)
Tavg <- Tavg %>% group_by(Y, M) %>% summarise(Tmean = mean(Tmean, na.rm = TRUE), .groups = "drop")
Tmax <- Tmax %>% group_by(Y, M) %>% summarise(Tmax  = mean(Tmax,  na.rm = TRUE), .groups = "drop")
Tmin <- Tmin %>% group_by(Y, M) %>% summarise(Tmin  = mean(Tmin,  na.rm = TRUE), .groups = "drop")
Prec <- Prec %>% group_by(Y, M) %>% summarise(Prec  = sum(Prec,  na.rm = TRUE), .groups = "drop")

################################################################################
# 3) SPEI_monthly(): monthly-scale SPEI from monthly water deficit series
################################################################################
# Aggregates WD over 'scale' months using a kernel, fits a distribution by month-of-year,
# and converts CDF probabilities to standardized normal values (SPEI-like index).
#
# Notes:
# - Expects wd_data to contain column "wd" and ideally be ordered by time.
# - Uses monthly seasonality (fr = 12) via ts frequency.
# - ref.start/ref.end follow ts window() conventions: c(year, month).
SPEI_monthly <- function(wd_data, scale = 6, kernel = list(type = "rectangular", shift = 0),
                         distribution = "log-Logistic", fit = "ub-pwm", na.rm = FALSE,
                         ref.start = NULL, ref.end = NULL, x = FALSE, params = NULL) {
  
  # -------------------------
  # 1) Basic checks + coerce
  # -------------------------
  if (!is.data.frame(wd_data)) stop("wd_data must be a data.frame or tibble.")
  wd_data <- as.data.frame(wd_data)  # important: avoids tibble subsetting surprises
  
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
  # 2) Build a monthly ts object
  # -------------------------
  # Use [[ ]] so we always get a vector (works for data.frame AND tibble)
  water_deficit <- as.numeric(wd_data[["wd"]])
  
  if (all(c("Y", "M") %in% names(wd_data))) {
    
    # Order by year-month to ensure correct time sequence
    ord <- order(wd_data[["Y"]], wd_data[["M"]])
    wd_data <- wd_data[ord, , drop = FALSE]
    water_deficit <- as.numeric(wd_data[["wd"]])
    
    start_ym <- c(wd_data[["Y"]][1], wd_data[["M"]][1])
    df_water_deficit <- ts(as.matrix(water_deficit), frequency = 12, start = start_ym)
    
  } else {
    # Fallback when Y/M not available (still assumes regular monthly sequence)
    df_water_deficit <- ts(as.matrix(water_deficit), frequency = 12)
  }
  
  m <- ncol(df_water_deficit)
  fr <- frequency(df_water_deficit)  # should be 12
  
  # -------------------------
  # 3) Parameter container
  # -------------------------
  coef <- switch(distribution,
                 "Gamma"        = array(NA, c(2, m, fr), list(par = c("alpha", "beta"), colnames(df_water_deficit), NULL)),
                 "log-Logistic" = array(NA, c(3, m, fr), list(par = c("xi", "alpha", "kappa"), colnames(df_water_deficit), NULL)),
                 "PearsonIII"   = array(NA, c(3, m, fr), list(par = c("mu", "sigma", "gamma"), colnames(df_water_deficit), NULL))
  )
  
  dim_one <- ifelse(distribution == "Gamma", 2, 3)
  
  if (!is.null(params)) {
    if (dim(params)[1] != dim_one || dim(params)[2] != m || dim(params)[3] != fr) {
      stop(sprintf("params must have dimensions (%d, %d, %d)", dim_one, m, fr))
    }
  }
  
  # -------------------------
  # 4) Reference period window
  # -------------------------
  df_fit <- if (!is.null(ref.start) && !is.null(ref.end)) {
    window(df_water_deficit, ref.start, ref.end)
  } else {
    df_water_deficit
  }
  
  std <- df_water_deficit * NA
  
  # -------------------------
  # 5) Main loop
  # -------------------------
  for (s in 1:m) {
    
    acu_fit  <- as.numeric(df_fit[, s])
    acu_pred <- as.numeric(df_water_deficit[, s])
    
    # --- aggregation over 'scale' months
    if (scale > 1) {
      wgt <- kern(scale, kernel$type, kernel$shift)
      
      if (length(acu_fit) >= scale) {
        tmp <- rowSums(embed(acu_fit, scale) * wgt, na.rm = na.rm)
        acu_fit_agg <- rep(NA, length(acu_fit))
        acu_fit_agg[scale:length(acu_fit)] <- tmp
      } else {
        acu_fit_agg <- rep(NA, length(acu_fit))
      }
      
      if (length(acu_pred) >= scale) {
        tmp2 <- rowSums(embed(acu_pred, scale) * wgt, na.rm = na.rm)
        acu_pred_agg <- rep(NA, length(acu_pred))
        acu_pred_agg[scale:length(acu_pred)] <- tmp2
      } else {
        acu_pred_agg <- rep(NA, length(acu_pred))
      }
      
    } else {
      acu_fit_agg  <- acu_fit
      acu_pred_agg <- acu_pred
    }
    
    # keep cycle() stable by creating ts objects
    acu_fit_ts  <- ts(acu_fit_agg,  frequency = fr, start = start(df_fit))
    acu_pred_ts <- ts(acu_pred_agg, frequency = fr, start = start(df_water_deficit))
    
    for (c in 1:fr) {
      
      f  <- which(cycle(acu_fit_ts) == c)
      ff <- which(cycle(acu_pred_ts) == c)
      
      fit_vals <- acu_fit_ts[f]
      fit_vals <- fit_vals[!is.na(fit_vals)]
      
      if (length(fit_vals) == 0) {
        std[ff, s] <- NA
        coef[, s, c] <- NA
        next()
      }
      
      month <- sort.default(as.numeric(fit_vals), method = "quick")
      
      # --- estimate parameters
      if (is.null(params)) {
        
        month_sd <- sd(month, na.rm = TRUE)
        if (is.na(month_sd) || month_sd == 0) {
          std[ff, s] <- NA
          coef[, s, c] <- NA
          next()
        }
        
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
        
        pwm <- switch(fit,
                      "pp-pwm"  = pwm.pp(month, -0.35, 0, nmom = 3),
                      "ub-pwm"  = TLMoments::PWM(month, order = 0:2),
                      "max-lik" = TLMoments::PWM(month, order = 0:2)
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
        
        if (distribution == "log-Logistic" && fit == "max-lik") {
          f_params <- parglo.maxlik(month, f_params)$para
        }
        
      } else {
        f_params <- as.vector(params[, s, c])
        pze <- 0
      }
      
      # --- CDF -> qnorm
      pred_vals <- as.numeric(acu_pred_ts[ff])
      
      cdf_res <- switch(distribution,
                        "log-Logistic" = lmom::cdfglo(pred_vals, f_params),
                        "Gamma"        = lmom::cdfgam(pred_vals, f_params),
                        "PearsonIII"   = lmom::cdfpe3(pred_vals, f_params)
      )
      
      cdf_res <- pmin(pmax(cdf_res, 1e-6), 1 - 1e-6)
      
      std[ff, s] <- qnorm(cdf_res)
      coef[, s, c] <- f_params
      
      if (distribution != "log-Logistic") {
        std[ff, s] <- qnorm(pze + (1 - pze) * pnorm(std[ff, s]))
      }
    }
  }
  
  colnames(std) <- colnames(df_water_deficit)
  
  z <- list(
    call = match.call(expand.dots = FALSE),
    fitted = std,
    coefficients = coef,
    scale = scale,
    kernel = list(type = kernel$type, shift = kernel$shift, values = kern(scale, kernel$type, kernel$shift)),
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
# 4) Monthly PET, WD, SPEI scales, and correlation heatmap
################################################################################

# Monthly temperature range (clamp negatives)
T_diff <- as.numeric(Tmax$Tmax) - as.numeric(Tmin$Tmin)
T_diff <- pmax(T_diff, 0)

# Monthly PET:
# You used SPEI::hargreaves() (fine for monthly); keep that behavior.
PET <- hargreaves(Tmin$Tmin, Tmax$Tmax, lat = lat, na.rm = TRUE)

# Monthly water deficit
WD_values <- Prec$Prec - PET

# WD data frame for SPEI_monthly (keep Y/M structure)
WD <- Tavg
WD$Tmean <- NULL
WD$wd <- WD_values

# Compute correlations for multiple month scales
temporal_matrix_list <- list()
place_holder <- 1

for (ij in lower_limit:upper_limit) {
  
  temp_rwl_subset <- temp_rwl
  
  SPEI_temp <- SPEI_monthly(wd_data = WD, scale = ij)
  
  # Convert to wide: rows = year, cols = month (1..12)
  SPEI_temp <- dcast(formula = Y ~ M, value.var = "SPEI", data = SPEI_temp)
  row.names(SPEI_temp) <- SPEI_temp$Y
  SPEI_temp$Y <- NULL
  
  # Optionally include previous year months
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
  
  # Merge TRW and SPEI by year
  ncol_temp_rwl_subset <- ncol(temp_rwl_subset)
  colnames_temp_rwl_subset <- colnames(temp_rwl_subset)
  
  SPEI_temp$temp_year <- row.names(SPEI_temp)
  temp_rwl_subset$temp_year <- row.names(temp_rwl_subset)
  
  temporal_data <- merge(temp_rwl_subset, SPEI_temp, by = "temp_year")
  
  temp_rwl_subset <- data.frame(
    temporal_data[, c(2:(1 + ncol_temp_rwl_subset))],
    row.names = temporal_data$temp_year
  )
  colnames(temp_rwl_subset) <- colnames_temp_rwl_subset
  
  SPEI_temp <- data.frame(
    temporal_data[, c((1 + ncol_temp_rwl_subset + 1):ncol(temporal_data))],
    row.names = temporal_data$temp_year
  )
  
  # Correlation vector for this scale: one value per (month column position)
  temporal_matrix <- matrix(NA, nrow = 1, ncol = ncol(SPEI_temp))
  
  for (j in 0:(ncol(SPEI_temp) - ij)) {
    
    x <- SPEI_temp[, (j + 1)]
    x <- matrix(x, nrow = nrow(SPEI_temp), ncol = 1)
    
    temporal_correlation <- cor(
      temp_rwl_subset[, 1], x[, 1],
      method = "pearson",
      use = "pairwise.complete.obs"
    )
    
    temporal_matrix[1, j + ij] <- temporal_correlation
  }
  
  temporal_matrix_list[[place_holder]] <- temporal_matrix
  place_holder <- place_holder + 1
}

monthly_SPEI_correlations <- data.frame(do.call(rbind, temporal_matrix_list))

################################################################################
# 5) Plot heatmap (ggplot2)
################################################################################

monthly_SPEI_correlations$season_length <- seq(lower_limit, upper_limit)
melted <- melt(monthly_SPEI_correlations, id.vars = c("season_length"))

# Remove weak correlations
melted$value <- ifelse(abs(melted$value) < 0.25, NA, melted$value)

ggplot(melted, aes_(x = ~as.numeric(variable), y = ~season_length, fill = ~value)) +
  geom_tile() +
  xlab("Month") +
  ylab("Season Length") +
  scale_x_continuous(expand = c(0, 0), breaks = 1:12) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(1, 24, by = 2)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", na.value = "gray97", midpoint = 0) +
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
    panel.background = element_rect(fill = "gray97", colour = "gray80", size = 0.5, linetype = "solid")
  )

ggsave("SPEI_example_monthly2.png", height = 7, width = 10)
