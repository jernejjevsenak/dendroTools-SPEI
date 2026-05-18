# Open the required libraries
library("dplyr")
library("lubridate")
library("reshape2")
library("dendroTools")
library("SPEI")
library("zoo")
library("dplR")
library("lmomco")
library("ggplot2")
library("TLMoments")
library("lmom")

################################################################################
# Open your data and define the key arguments
################################################################################

lat <- 51.0386       # latitude of your site
lower_limit <- 1     # lowest number of months used to aggregate water deficit into SPEI
upper_limit <- 12    # maximum number of months used to aggregate water deficit into SPEI
previous_year <- FALSE # should previous year be considered?

################################################################################
# Open TRWi data
################################################################################

temp_rwl <- read.crn("my_chron.crn")
temp_rwl$samp.depth <- NULL # remove the sample depth if present

################################################################################
# Open climate data
################################################################################

Tavg <- read.table("my_chron_Tmean.csv", sep = ",", header = TRUE)
Tmax <- read.table("my_chron_Tmax.csv", sep = ",", header = TRUE)
Tmin <- read.table("my_chron_Tmin.csv", sep = ",", header = TRUE)
Prec <- read.table("my_chron_Psum.csv", sep = ",", header = TRUE)

################################################################################
# Convert daily climate data to monthly climate data
################################################################################

Tavg <- Tavg %>%
  group_by(Y, M) %>%
  summarise(Tmean = mean(Tmean, na.rm = TRUE), .groups = "drop")

Tmax <- Tmax %>%
  group_by(Y, M) %>%
  summarise(Tmax = mean(Tmax, na.rm = TRUE), .groups = "drop")

Tmin <- Tmin %>%
  group_by(Y, M) %>%
  summarise(Tmin = mean(Tmin, na.rm = TRUE), .groups = "drop")

Prec <- Prec %>%
  group_by(Y, M) %>%
  summarise(Prec = sum(Prec, na.rm = TRUE), .groups = "drop")

# Sort all monthly climate data by year and month
Tavg <- arrange(Tavg, Y, M)
Tmax <- arrange(Tmax, Y, M)
Tmin <- arrange(Tmin, Y, M)
Prec <- arrange(Prec, Y, M)

# Check that all climate data frames have identical year-month structure
stopifnot(
  all(Tavg$Y == Tmax$Y, Tavg$M == Tmax$M),
  all(Tavg$Y == Tmin$Y, Tavg$M == Tmin$M),
  all(Tavg$Y == Prec$Y, Tavg$M == Prec$M)
)

################################################################################
# SPEI_monthly function
################################################################################

SPEI_monthly <- function(wd_data,
                         scale = 6,
                         kernel = list(type = "rectangular", shift = 0),
                         distribution = "log-Logistic",
                         fit = "ub-pwm",
                         na.rm = FALSE,
                         ref.start = NULL,
                         ref.end = NULL,
                         x = FALSE,
                         params = NULL) {
  
  # ---------------------------------------------------------------------------
  # Basic checks and preparation
  # ---------------------------------------------------------------------------
  
  if (!all(c("Y", "M", "wd") %in% names(wd_data))) {
    stop("wd_data must contain columns: Y, M, and wd")
  }
  
  wd_data <- wd_data[order(wd_data$Y, wd_data$M), ]
  
  if (nrow(wd_data) %% 12 != 0) {
    warning("Number of rows is not a multiple of 12. Check whether all years have complete monthly data.")
  }
  
  water_deficit <- wd_data[, "wd"]
  
  df_water_deficit <- ts(
    as.matrix(water_deficit),
    frequency = 12,
    start = c(wd_data$Y[1], wd_data$M[1])
  )
  
  scale <- as.numeric(scale)
  na.rm <- as.logical(na.rm)
  x <- as.logical(x)
  
  if (!(distribution %in% c("log-Logistic", "Gamma", "PearsonIII"))) {
    stop('Distrib must be one of "log-Logistic", "Gamma" or "PearsonIII"')
  }
  
  if (!(fit %in% c("max-lik", "ub-pwm", "pp-pwm"))) {
    stop('Method must be one of "ub-pwm", "pp-pwm" or "max-lik"')
  }
  
  if ((!is.null(ref.start) && length(ref.start) != 2) |
      (!is.null(ref.end) && length(ref.end) != 2)) {
    stop("Start and end of the reference period must be a numeric vector of length two.")
  }
  
  m <- ncol(df_water_deficit)
  fr <- frequency(df_water_deficit)
  
  # ---------------------------------------------------------------------------
  # Prepare coefficient array
  # ---------------------------------------------------------------------------
  
  coef <- switch(
    distribution,
    "Gamma" = array(
      NA,
      c(2, m, fr),
      list(par = c("alpha", "beta"), colnames(df_water_deficit), NULL)
    ),
    "log-Logistic" = array(
      NA,
      c(3, m, fr),
      list(par = c("xi", "alpha", "kappa"), colnames(df_water_deficit), NULL)
    ),
    "PearsonIII" = array(
      NA,
      c(3, m, fr),
      list(par = c("mu", "sigma", "gamma"), colnames(df_water_deficit), NULL)
    )
  )
  
  dim_one <- ifelse(distribution == "Gamma", 2, 3)
  
  if (!is.null(params)) {
    if (dim(params)[1] != dim_one | dim(params)[2] != m | dim(params)[3] != fr) {
      stop(paste0(
        "parameters array should have dimensions (",
        dim_one, ", ", m, ", ", fr, ")"
      ))
    }
  }
  
  # ---------------------------------------------------------------------------
  # Reference period
  # ---------------------------------------------------------------------------
  
  if (!is.null(ref.start) && !is.null(ref.end)) {
    df_water_deficit.fit <- window(df_water_deficit, ref.start, ref.end)
  } else {
    df_water_deficit.fit <- df_water_deficit
  }
  
  std <- df_water_deficit * NA
  
  # ---------------------------------------------------------------------------
  # Loop through series
  # ---------------------------------------------------------------------------
  
  for (s in 1:m) {
    
    acu <- df_water_deficit.fit[, s]
    acu.pred <- df_water_deficit[, s]
    
    # Aggregate water deficit over selected scale
    if (scale > 1) {
      wgt <- kern(scale, kernel$type, kernel$shift)
      
      acu[scale:length(acu)] <- rowSums(embed(acu, scale) * wgt, na.rm = na.rm)
      acu[1:(scale - 1)] <- NA
      
      acu.pred[scale:length(acu.pred)] <- rowSums(embed(acu.pred, scale) * wgt, na.rm = na.rm)
      acu.pred[1:(scale - 1)] <- NA
    }
    
    # Loop through months
    for (c in 1:fr) {
      
      f <- which(cycle(acu) == c)
      f <- f[!is.na(acu[f])]
      
      ff <- which(cycle(acu.pred) == c)
      ff <- ff[!is.na(acu.pred[ff])]
      
      month <- sort.default(acu[f], method = "quick")
      
      if (length(month) == 0) {
        std[ff, s] <- NA
        next
      }
      
      if (is.null(params)) {
        
        month_sd <- sd(month, na.rm = TRUE)
        
        if (is.na(month_sd) || month_sd == 0) {
          std[ff, s] <- NA
          next
        }
        
        if (distribution != "log-Logistic") {
          pze <- sum(month == 0) / length(month)
          month <- month[month > 0]
        }
        
        if (length(month) < 4) {
          std[ff, s] <- NA
          coef[, s, c] <- NA
          next
        }
        
        pwm <- switch(
          fit,
          "pp-pwm" = pwm.pp(month, -0.35, 0, nmom = 3),
          "ub-pwm" = TLMoments::PWM(month, order = 0:2),
          "max-lik" = TLMoments::PWM(month, order = 0:2)
        )
        
        lmom <- pwm2lmom(pwm)
        
        if (!are.lmom.valid(lmom) ||
            anyNA(lmom[[1]]) ||
            any(is.nan(lmom[[1]]))) {
          std[ff, s] <- NA
          coef[, s, c] <- NA
          next
        }
        
        fortran_vec <- c(lmom$lambdas[1:2], lmom$ratios[3])
        
        f_params <- switch(
          distribution,
          "log-Logistic" = tryCatch(
            lmom::pelglo(fortran_vec),
            error = function(e) {
              parglo(lmom)$para
            }
          ),
          "Gamma" = tryCatch(
            lmom::pelgam(fortran_vec),
            error = function(e) {
              pargam(lmom)$para
            }
          ),
          "PearsonIII" = tryCatch(
            lmom::pelpe3(fortran_vec),
            error = function(e) {
              parpe3(lmom)$para
            }
          )
        )
        
        if (distribution == "log-Logistic" && fit == "max-lik") {
          f_params <- parglo.maxlik(month, f_params)$para
        }
        
      } else {
        
        f_params <- as.vector(params[, s, c])
        
      }
      
      cdf_res <- switch(
        distribution,
        "log-Logistic" = lmom::cdfglo(acu.pred[ff], f_params),
        "Gamma" = lmom::cdfgam(acu.pred[ff], f_params),
        "PearsonIII" = lmom::cdfpe3(acu.pred[ff], f_params)
      )
      
      # Avoid Inf values from exact 0 or 1 probabilities
      cdf_res <- pmin(pmax(cdf_res, .Machine$double.eps), 1 - .Machine$double.eps)
      
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
    kernel = list(
      type = kernel$type,
      shift = kernel$shift,
      values = kern(scale, kernel$type, kernel$shift)
    ),
    distribution = distribution,
    fit = fit,
    na.action = na.rm
  )
  
  if (isTRUE(x)) {
    z$df_water_deficit <- df_water_deficit
  }
  
  if (!is.null(ref.start)) {
    z$ref.period <- rbind(ref.start, ref.end)
  }
  
  wd_data$SPEI <- as.numeric(z$fitted)
  
  if (isTRUE(x)) {
    z$wd_data <- wd_data
    return(z)
  } else {
    return(wd_data)
  }
}

################################################################################
# Calculate PET and water deficit
################################################################################

# Create monthly time-series objects for Hargreaves PET
Tmin_ts <- ts(
  Tmin$Tmin,
  frequency = 12,
  start = c(Tmin$Y[1], Tmin$M[1])
)

Tmax_ts <- ts(
  Tmax$Tmax,
  frequency = 12,
  start = c(Tmax$Y[1], Tmax$M[1])
)

# Calculate PET using the Hargreaves method
PET <- hargreaves(Tmin_ts, Tmax_ts, lat = lat, na.rm = TRUE)
PET <- as.numeric(PET)

# Calculate water deficit
WD_values <- Prec$Prec - PET

# Create water-deficit data frame
WD <- Tavg
WD$Tmean <- NULL
WD$wd <- WD_values

################################################################################
# Optional diagnostic check
################################################################################

test <- SPEI_monthly(WD, scale = 3, x = TRUE)
print(frequency(test$df_water_deficit)) # should be 12

################################################################################
# Loop from lower to upper window and calculate monthly SPEI correlations
################################################################################

temporal_matrix_list <- list()
place_holder <- 1

for (ij in lower_limit:upper_limit) {
  
  temp_rwl_subset <- temp_rwl
  
  SPEI_temp <- SPEI_monthly(wd_data = WD, scale = ij)
  
  SPEI_temp <- dcast(formula = Y ~ M, value.var = "SPEI", data = SPEI_temp)
  row.names(SPEI_temp) <- SPEI_temp$Y
  SPEI_temp$Y <- NULL
  
  # Here we rearrange water deficit data if the previous year is considered
  if (previous_year == TRUE) {
    
    SPEI_temp$temp_year <- row.names(SPEI_temp)
    SPEI_temp <- dplyr::arrange(SPEI_temp, desc(temp_year))
    SPEI_temp <- years_to_rownames(SPEI_temp, "temp_year")
    
    SPEI_temp_previous <- SPEI_temp[-1, , F]
    SPEI_temp_current <- SPEI_temp[-nrow(SPEI_temp), , F]
    
    row_names_current <- row.names(SPEI_temp_current)
    
    SPEI_temp <- cbind(SPEI_temp_previous, SPEI_temp_current)
    SPEI_temp <- data.frame(SPEI_temp)
    row.names(SPEI_temp) <- row_names_current
  }
  
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
  
  temporal_matrix <- matrix(NA, nrow = 1, ncol = ncol(SPEI_temp))
  
  for (j in 0:(ncol(SPEI_temp) - ij)) {
    
    # Use the END month of the SPEI window.
    # Example:
    # ij = 3 and j = 0 means Jan-Feb-Mar SPEI is taken from March column
    # and stored under March.
    x <- SPEI_temp[, (j + ij)]
    x <- matrix(x, nrow = nrow(SPEI_temp), ncol = 1)
    
    temporal_correlation <- cor(
      temp_rwl_subset[, 1],
      x[, 1],
      method = "pearson",
      use = "pairwise.complete.obs"
    )
    
    # Store the correlation at the END month of the climate window
    temporal_matrix[1, j + ij] <- temporal_correlation
  }
  
  temporal_matrix_list[[place_holder]] <- temporal_matrix
  place_holder <- place_holder + 1
  
  print(ij)
}

monthly_SPEI_correlations <- data.frame(do.call(rbind, temporal_matrix_list))

################################################################################
# Visualization with ggplot2
################################################################################

monthly_SPEI_correlations$season_length <- seq(lower_limit, upper_limit)
melted <- melt(monthly_SPEI_correlations, id.vars = c("season_length"))

# Optional: remove weak correlations
# melted$value <- ifelse(abs(melted$value) < 0.25, NA, melted$value)

ggplot(melted, aes_(x = ~as.numeric(variable), y = ~season_length, fill = ~value)) +
  geom_tile() +
  xlab("Month") +
  ylab("Season Length") +
  scale_x_continuous(expand = c(0, 0), breaks = seq(1, 12)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(1, 24, by = 2)) +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    na.value = "gray97",
    midpoint = 0
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

ggsave("SPEI_example_monthly.png", height = 7, width = 10)