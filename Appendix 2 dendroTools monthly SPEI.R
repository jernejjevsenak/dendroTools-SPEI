# Open the required libraries
library("dplyr")
library("lubridate")
library("reshape2")
library("dendroTools")
library("SPEI")
library("zoo")
library("dplR")

##############################################################################################
# Open your data and define the key arguments
lat = 51.0386 # this is latitude of your site
lower_limit = 1 # this is the lowest number of days that will be used to aggregate water deficit into SPEI
upper_limit = 24  # this is maximum number of days that will be used to aggregate water deficit into SPEI
previous_year = TRUE # should previous year be considered?

# Open TRWi data
temp_rwl <- read.crn("my_chron.crn")
temp_rwl$samp.depth <- NULL # remove the sample depth if present

# Open climate data
Tavg <- read.table("my_chron_Tmean.csv", sep = ",", header = TRUE)
Tmax <- read.table("my_chron_Tmax.csv", sep = ",", header = TRUE)
Tmin <- read.table("my_chron_Tmin.csv", sep = ",", header = TRUE)
Prec <- read.table("my_chron_Psum.csv", sep = ",", header = TRUE)

# convert to monthly
Tavg <- dplyr::group_by(Tavg, Y, M) %>% summarise(Tmean = mean(Tmean, na.rm = TRUE))
Tmax <- dplyr::group_by(Tmax, Y, M) %>% summarise(Tmax = mean(Tmax, na.rm = TRUE))
Tmin <- dplyr::group_by(Tmin, Y, M) %>% summarise(Tmin = mean(Tmin, na.rm = TRUE))
Prec <- dplyr::group_by(Prec, Y, M) %>% summarise(Prec = sum(Prec, na.rm = TRUE))

################################################################################
################################################################################

# SPEI_daily function will aggregate the water deficit (wd) series into SPEI
# The argument scale is used to define the number of days to be used for the aggregation
# All credits go to the authors of SPEI R package (https://cran.r-project.org/web/packages/SPEI/index.html)

SPEI_monthly <- function(wd_data, scale = 6, kernel=list(type='rectangular',shift=0),
					   distribution='log-Logistic', fit='ub-pwm', na.rm=FALSE, ref.start=NULL,
					   ref.end=NULL, x=FALSE, params=NULL){

  water_deficit <- wd_data[,"wd" ]
  water_deficit = as.ts(zooreg(water_deficit))
  df_water_deficit = water_deficit

  scale <- as.numeric(scale)
  na.rm <- as.logical(na.rm)
  x <- as.logical(x)

  if (!(distribution %in% c('log-Logistic', 'Gamma', 'PearsonIII'))) {
    stop('Distrib must be one of "log-Logistic", "Gamma" or "PearsonIII"')
  }
  if (!(fit %in% c('max-lik', 'ub-pwm', 'pp-pwm'))) {
    stop('Method must be one of "ub-pwm" (default), "pp-pwm" or "max-lik"')
  }
  if ( (!is.null(ref.start) && length(ref.start)!=2) | (!is.null(ref.end) && length(ref.end)!=2) ) {
    stop('Start and end of the reference period must be a numeric vector of length two.')
  }

  if (!is.ts(df_water_deficit)) {
    # frequency and start can be random
    df_water_deficit <- ts(as.matrix(df_water_deficit), frequency = 20, start = c(1950, 5))
  } else {
    df_water_deficit <- ts(as.matrix(df_water_deficit), frequency=frequency(df_water_deficit), start=start(df_water_deficit))
  }
  
  m <- ncol(df_water_deficit)
  fr <- frequency(df_water_deficit)

  coef = switch(distribution,
                "Gamma" = array(NA,c(2,m,fr),list(par=c('alpha','beta'),colnames(df_water_deficit),NULL)),
                "log-Logistic" = array(NA,c(3,m,fr),list(par=c('xi','alpha','kappa'),colnames(df_water_deficit),NULL)),
                "PearsonIII" = coef <- array(NA,c(3,m,fr),list(par=c('mu','sigma','gamma'),colnames(df_water_deficit),NULL))
  )

  dim_one = ifelse(distribution == "Gamma", 2, 3)

  if (!is.null(params)) {
    if (dim(params)[1]!=dim_one | dim(params)[2]!=m | dim(params)[3]!=12) {
      stop(paste('parameters array should have dimensions (', dim_one, ', ', m, ', 12)',sep=' '))
    }
  }

  # Loop through series (columns in df_water_deficit)
  if (!is.null(ref.start) && !is.null(ref.end)) {
    df_water_deficit.fit <- window(df_water_deficit,ref.start,ref.end)
  } else {
    df_water_deficit.fit <- df_water_deficit
  }
  std <- df_water_deficit*NA
  for (s in 1:m) {
    # Cumulative series (acu)
    acu <- df_water_deficit.fit[,s]
    acu.pred <- df_water_deficit[,s]
    if (scale>1) {
      wgt <- kern(scale,kernel$type,kernel$shift)
      acu[scale:length(acu)] <- rowSums(embed(acu,scale)*wgt,na.rm=na.rm)
      acu[1:(scale-1)] <- NA
      acu.pred[scale:length(acu.pred)] <- rowSums(embed(acu.pred,scale)*wgt,na.rm=na.rm)
      acu.pred[1:(scale-1)] <- NA
    }

    # Loop through the months/days
    for (c in (1:fr)) {
      # Filter month m, excluding NAs
      f <- which(cycle(acu)==c)
      f <- f[!is.na(acu[f])]
      ff <- which(cycle(acu.pred)==c)
      ff <- ff[!is.na(acu.pred[ff])]

      # Daily series, sorted
      month <- sort.default(acu[f], method="quick")

      if (length(month)==0) {
        std[f] <- NA
        next()
      }

      if (is.null(params)) {
        month_sd = sd(month,na.rm=TRUE)
        if (is.na(month_sd) || (month_sd == 0)) {
          std[f] <- NA
          next
        }

        if(distribution != "log-Logistic"){
          pze <- sum(month==0)/length(month)
          month = month[month > 0]
        }

        # Stop early and assign NAs if month's df_water_deficit is length < 4
        if(length(month) < 4){
          std[ff,s] = NA
          coef[,s,c] <- NA
          next
        }

        # Calculate probability weighted moments based on fit with lmomco or TLMoments
        pwm = switch(fit,
                     "pp-pwm" = pwm.pp(month,-0.35,0, nmom=3),
                     #pwm.ub(month, nmom=3)
                     TLMoments::PWM(month, order=0:2)
        )

        # Check L-moments validity
        lmom <- pwm2lmom(pwm)
        if ( !are.lmom.valid(lmom) || anyNA(lmom[[1]]) || any(is.nan(lmom[[1]])) ){
          next
        }

        # lmom fortran functions need specific inputs L1, L2, T3
        # this is handled by lmomco internally with lmorph
        fortran_vec = c(lmom$lambdas[1:2], lmom$ratios[3])

        # Calculate parameters based on distribution with lmom then lmomco
        f_params = switch(distribution,
                          "log-Logistic" = tryCatch(lmom::pelglo(fortran_vec), error = function(e){ parglo(lmom)$para }),
                          "Gamma" = tryCatch(lmom::pelgam(fortran_vec), error = function(e){ pargam(lmom)$para }),
                          "PearsonIII" = tryCatch(lmom::pelpe3(fortran_vec), error = function(e){ parpe3(lmom)$para })
        )

        # Adjust if user chose log-Logistic and max-lik
        if(distribution == 'log-Logistic' && fit=='max-lik'){
          f_params = parglo.maxlik(month, f_params)$para
        }
      } else {

        f_params = as.vector(params[,s,c])

      }

      # Calculate cdf based on distribution with lmom
      cdf_res = switch(distribution,
                       "log-Logistic" = lmom::cdfglo(acu.pred[ff], f_params),
                       "Gamma" = lmom::cdfgam(acu.pred[ff], f_params),
                       "PearsonIII" = lmom::cdfpe3(acu.pred[ff], f_params)
      )

      std[ff,s] = qnorm(cdf_res)
      coef[,s,c] <- f_params

      # Adjust if user chose Gamma or PearsonIII
      if(distribution != 'log-Logistic'){
        std[ff,s] = qnorm(pze + (1-pze)*pnorm(std[ff,s]))
      }

    } # next c (month)
  } # next s (series)
  colnames(std) <- colnames(df_water_deficit)

  z <- list(call=match.call(expand.dots=FALSE),
            fitted=std,coefficients=coef,scale=scale,
            kernel=list(type=kernel$type, shift=kernel$shift,values=kern(scale,kernel$type,kernel$shift)),
            distribution=distribution,fit=fit,na.action=na.rm)
  
  if (x) z$df_water_deficit <- df_water_deficit
  if (!is.null(ref.start)) z$ref.period <- rbind(ref.start,ref.end)

  wd_data$SPEI <- as.numeric(z$fitted)
  
  return(wd_data)
}



# Calculate T_diff
T_diff <- as.numeric(Tmax$Tmax) - as.numeric(Tmin$Tmin)
T_diff = ifelse(T_diff < 0, 0, T_diff) # If Tmin is greater than Tmax, set the difference to 0

# Calculate PET using the Hargreaves method

PET <- hargreaves(Tmin$Tmin, Tmax$Tmax, Ra = NA, lat = lat, Pre = Prec$Prec, na.rm = TRUE)

key_df <- data.frame(year_org = unique(Tmin$Y),
                     year = seq(1, length(unique(Tmin$Y)), by = 1)
)


PET <- data.frame(Y=as.matrix(PET), date=as.Date(as.yearmon(time(PET)))) %>%
  mutate(year = year(date), month = month(date)) %>% left_join(key_df, by = "year") %>%
  select(ET0_har, year_org, month) %>% rename(year = year_org)


# Calculate Water Deficit (WD)
WD_values = Prec$Prec - PET$ET0_har

# We take a random data frame (here Tavg) and convert it into data frame for WD
WD <- Tavg
WD$Tmean <- NULL
WD$wd<- WD_values

# Loop from lower to upper window and aggregate water deficit in SPEI of different scales
temporal_matrix_list <- list() # this is a list which will be used to save calculated SPEI of different lengths (scales)
place_holder = 1

for (ij in (lower_limit:upper_limit)){

  temp_rwl_subset <- temp_rwl
  
  SPEI_temp <- SPEI_monthly(wd_data = WD, scale = ij)
  SPEI_temp <- dcast(formula = Y ~ M, value.var = "SPEI", data = SPEI_temp) # convert into wide data frame
  row.names(SPEI_temp) <- SPEI_temp$Y
  SPEI_temp$Y <- NULL
  # Here we rearrange water deficit data if the previous year is considered
  if (previous_year == TRUE) {

    SPEI_temp$temp_year <- row.names(SPEI_temp)
    SPEI_temp <- dplyr::arrange(SPEI_temp, desc(temp_year))
    SPEI_temp <- years_to_rownames(SPEI_temp, "temp_year")
    SPEI_temp_previous <- SPEI_temp[-1, , F]
    SPEI_temp_current <- SPEI_temp[-nrow(SPEI_temp), ,F]
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

  temp_rwl_subset <- data.frame(temporal_data[, c(2:(1 + ncol_temp_rwl_subset))], row.names = temporal_data$temp_year)
  colnames(temp_rwl_subset) <- colnames_temp_rwl_subset

  SPEI_temp <- data.frame(temporal_data[, c((1 + ncol_temp_rwl_subset + 1):ncol(temporal_data))],
                         row.names = temporal_data$temp_year)

  temporal_matrix <- matrix(NA, nrow = 1, ncol = (ncol(SPEI_temp)))


  for (j in 0: (ncol(SPEI_temp) - ij)) {

    x <- SPEI_temp[, (j + 1)]
    x <- matrix(x, nrow = nrow(SPEI_temp), ncol = 1)

    temporal_correlation <- cor(temp_rwl_subset[, 1], x[, 1], method = 'pearson')
    temporal_matrix[1, j + ij] <- temporal_correlation

  }

  # save the calculations in temporal_matrix_list
  temporal_matrix_list[[place_holder]] <- temporal_matrix
  place_holder = place_holder + 1

}

monthly_SPEI_correlations <- data.frame(do.call(rbind, temporal_matrix_list))

################################################################################

# An example code for visualization with ggplot2

monthly_SPEI_correlations$season_length <- seq(lower_limit, upper_limit)
melted <- melt(monthly_SPEI_correlations, id.vars = c("season_length"))

# delete lower values
# melted$value <- ifelse(abs(melted$value) < 0.25, NA, melted$value)

ggplot(melted,aes_(x = ~as.numeric(variable), y = ~season_length, fill = ~value)) +
  geom_tile() +
  xlab("Month") +
  ylab("Season Length") +
  scale_y_continuous(expand=c(0,0), breaks = seq(1, 24, by = 2)) +
  geom_vline(xintercept = 12) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", na.value = 'gray97', midpoint = 0) +
  theme_minimal() +  theme(axis.text = element_text(size = 10),
                           axis.title.y = element_text(size = 18), text = element_text(size = 18),
                           axis.title.x = element_blank(),
                           plot.title = element_text(size = 16),
                           legend.title = element_blank(),
                           legend.position = "bottom", legend.key.width = unit(3, "line"),
                           panel.background = element_rect(fill = "gray97",
                                                           colour = "gray80",
                                                           size = 0.5, linetype = "solid"))

ggsave(paste0("SPEI_example_monthly_updated.png"), height = 7, width = 10)

