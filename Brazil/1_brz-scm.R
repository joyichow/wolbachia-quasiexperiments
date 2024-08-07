# Brazil data SCM

identifiers <- c("CACHOEIRA", "CHARITAS", "SAO FRANCISCO", "CAFUBA", "CAMBOINHAS", "ENGENHO DO MATO", 
                 "ITACOATIARA","ITAIPU","JACARE",
                 "JARDIM IMBUI","MARAVISTA","PIRATININGA","SANTO ANTONIO","SERRA GRANDE","BOA VIAGEM",
                 "CENTRO","FATIMA","GRAGOATA","ICARAI",
                 "INGA","MORRO DO ESTADO","PE PEQUENO",
                 "PONTA DA AREIA","SANTA ROSA","SAO DOMINGOS","VIRADOURO","VITAL BRASIL", 
                 "CUBANGO", "ENGENHOCA","FONSECA","SANTANA", "SAO LOURENCO")

source("/your/local/filepath/CSCM_helper_functions.R")

library(dplyr)
library(lubridate)
brazil.den <- read.csv("/your/local/filepath/Niteroi_NeighbourhoodData_ExclJurujuba_Apr2021.csv")
brazil.den <- brazil.den[c(1:3, 7)]

# Site_IDs
brazil.den <- brazil.den %>%
  left_join(brazil.den %>%
              distinct(neighborhood) %>%
              mutate(site_ID = row_number()), 
            by = "neighborhood")

brazil.den <- brazil.den %>%
  mutate(yearmonth = as.Date(paste0(yearmonth, "_01"), format="%Y_%m_%d"),
         quarter = paste(year(yearmonth), quarter(yearmonth), sep = "_"))

brazil.den <- brazil.den %>%
  group_by(neighborhood, quarter, site_ID, zonecode) %>%
  summarize(notified = sum(notified), .groups = 'drop')

brazil.den <- brazil.den %>%
  group_by(neighborhood) %>%
  mutate(notified_norm = (notified - min(notified)) / (max(notified) - min(notified)))

site_maxmins <- brazil.den %>%
  filter(site_ID %in% 1:51) %>%
  group_by(site_ID) %>%
  summarise(min_notified = min(notified), max_notified = max(notified), .groups = 'drop')

brazil.den <- brazil.den %>%
  mutate(time_ID = as.integer(factor(quarter, levels = unique(quarter))))

brazil.den$site_ID <- as.numeric(brazil.den$site_ID)

brazil.den <- as.data.frame(brazil.den)

brazil.den.placebo <- subset(brazil.den, site_ID >= 33)

# treatment identifiers 
treatment_ids <- c(1:32) # treatment + in-time placbeo
# treatment_ids <- c(33:51) # for in-space placebo

#treatment times
t_ints_real <- c(rep(42, 3), rep(43, 11), rep(45, 13), rep(51, 5), rep(45, 19))       #treatment + in-space placebo
# t_ints_real <- c(rep(42-5, 3), rep(43-5, 11), rep(45-5, 13), rep(51-5, 5), rep(45, 19)) #in-time placebo

registerDoParallel(cores = 8)

# parallelising SCM
brazil_results <- foreach(i = treatment_ids, .combine = "rbind", .export = c("countSynth", "csynth_inner", "factor.init", "factor.sim", 
                                                                                   "holdout.csynth", "brazil.interventions"), 
                                 .packages = c('dplyr', 'zoo', 'readxl', 'ggplot2', 'Synth', 'glmnet', 'dplyr', 
                                               'osqp', 'optimx', 'tidyr', 'lubridate', 'zoo', 'data.table')) %dopar% {
                                                 
                                                 source("/your/local/filepath/CSCM_helper_functions.R")
                                                 # tryCatch({
                                                   sink("brazil.txt", append=TRUE)
                                                   test <- brazil.interventions(i, brazil.den, t_ints_real)
                                                 # }, 
                                                 # error=function(...) NULL)
                                                 return(test)
                                               }

stopImplicitCluster()

brz_error_df <- do.call("rbind", brazil_results[,1])

write_xlsx(brz_error_df, "/your/local/filepath/transfer.xlsx")

# Run confidence intervals first (brazil-confints.R)
library(dplyr)
agg.IE.df <- do.call("rbind", brazil_results[,3])

agg.IE.df <- agg.IE.df %>%
  group_by(IE_df_x, IE_df_method) %>%
  summarise(mean_value = mean(IE_df_y, na.rm = TRUE))

mean_to_df <- function(df) {
  mean_values <- colMeans(df, na.rm = TRUE)
  df_mean <- as.data.frame(t(mean_values))
  return(df_mean)
}

max_cols_SCM <- max(sapply(confidence_intervals_SCM_normalized, ncol))
max_cols_CSCM <- max(sapply(confidence_intervals_CSCM_normalized, ncol))

SCM.ci.top <- do.call(rbind, lapply(confidence_intervals_SCM_normalized, extract_row, row_number = 1, max_cols = max_cols_SCM))
SCM.ci.bottom <- do.call(rbind, lapply(confidence_intervals_SCM_normalized, extract_row, row_number = 2, max_cols = max_cols_SCM))
CSCM.ci.top <- do.call(rbind, lapply(confidence_intervals_CSCM_normalized, extract_row, row_number = 1, max_cols = max_cols_CSCM))
CSCM.ci.bottom <- do.call(rbind, lapply(confidence_intervals_CSCM_normalized, extract_row, row_number = 2, max_cols = max_cols_CSCM))

SCM.ci.top.mean <- mean_to_df(SCM.ci.top)
SCM.ci.bottom.mean <- mean_to_df(SCM.ci.bottom)
CSCM.ci.top.mean <- mean_to_df(CSCM.ci.top)
CSCM.ci.bottom.mean <- mean_to_df(CSCM.ci.bottom)

mean(colMeans(SCM.ci.top, na.rm = TRUE))
mean(colMeans(SCM.ci.bottom, na.rm = TRUE))
mean(colMeans(CSCM.ci.top, na.rm = TRUE))
mean(colMeans(CSCM.ci.bottom, na.rm = TRUE))

# Aggregated intervention efficacy for Brazil
sum_cfdf_y_SCM <- 0
sum_cfdf_y_CSCM <- 0
sum_cdf_y_SCM <- 0
sum_cdf_y_CSCM <- 0

# Loop through the list of data frames
for (df in brazil_results[,2]) {
  # Filter the data for SCM and CSCM
  df_SCM <- df[df$cfdf_method == "SCM" & df$cfdf_x > 0, ]
  df_CSCM <- df[df$cfdf_method == "CSCM" & df$cfdf_x > 0, ]
  
  # Sum cfdf_y and cdf_y for each method
  sum_cfdf_y_SCM <- sum_cfdf_y_SCM + sum(df_SCM$cfdf_y)
  sum_cfdf_y_CSCM <- sum_cfdf_y_CSCM + sum(df_CSCM$cfdf_y)
  sum_cdf_y_SCM <- sum_cdf_y_SCM + sum(df_SCM$cdf_y)
}

cat("Total cfdf_y for SCM:", sum_cfdf_y_SCM, "\n")
cat("Total cfdf_y for CSCM:", sum_cfdf_y_CSCM, "\n")
cat("Total cdf_y for SCM:", sum_cdf_y_SCM, "\n")

IE_SCM <- ((sum_cfdf_y_SCM - sum_cdf_y_SCM)*100)/sum_cfdf_y_SCM
IE_CSCM <- ((sum_cfdf_y_CSCM - sum_cdf_y_SCM)*100)/sum_cfdf_y_CSCM


# Combine the rows into a new dataframe
brazil.sc.ie <- do.call(rbind, ie.extract)
brazil.sc.ie <- brazil.sc.ie[, c(1:3, 5, 4)]
write_xlsx(brazil.sc.ie, "/your/local/filepath/transfer.xlsx")

(sum(brazil.sc.ie$CSCM_cumsum) - sum(brazil.sc.ie$cases_cumsum))*100/sum(brazil.sc.ie$CSCM_cumsum)
(sum(brazil.sc.ie$SCM_cumsum) - sum(brazil.sc.ie$cases_cumsum))*100/sum(brazil.sc.ie$SCM_cumsum)
