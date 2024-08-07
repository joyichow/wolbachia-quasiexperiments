# Brazil data SCM
library(doParallel)
library(foreach)

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

brazil.interventions <- function(i, brazil.den, t_ints_real){
  print(paste("current site:", i))
  placebo_id <- c(33:51)
  
  dengue.main.res <- countSynth(data=brazil.den,                                 
                                predictors=NULL,
                                dependent="notified_norm",
                                unit.variable="site_ID",
                                time.variable = "time_ID",
                                treatment.identifier = i,
                                # controls.identifier = c(33:51),                   # for main analysis and in-time placebo 
                                controls.identifier = placebo_id[placebo_id != i],  # for in-space placebo
                                t_int = t_ints_real[i],
                                min_1se=F,
                                K=2)
  print("checkpoint A")
  
  synth.orig <- dengue.main.res$dataprep.main$Y0plot %*% as.numeric(dengue.main.res$unit.weight.SCM)
  cscm.cf <- dengue.main.res$dataprep.main$Y0plot %*% as.numeric(dengue.main.res$unit.weight.full.sample)
  year <- 1:53
  
  adf <- as.data.frame(cbind(year,synth.orig,rep("SCM",length(year))))
  bdf <- as.data.frame(cbind(year,cscm.cf,rep("CSCM",length(year))))
  cdf <- as.data.frame(cbind(year,dengue.main.res$dataprep.main$Y1plot,rep("Observed",length(year))))
  
  colnames(adf)=colnames(bdf)=colnames(cdf)=c("Quarter","Y","Method")
  
  print("checkpoint B")
  
  min_Y <- as.numeric(site_maxmins[i,2])
  max_Y <- as.numeric(site_maxmins[i,3])
  
  adf$Y <- as.numeric(adf$Y)*(max_Y - min_Y) + min_Y
  bdf$Y <- as.numeric(bdf$Y)*(max_Y - min_Y) + min_Y
  cdf$Y <- as.numeric(cdf$Y)*(max_Y - min_Y) + min_Y
  
  print("checkpoint C")
  
  scm_error <- (as.numeric(adf[1:t_ints_real[i], "Y"]) - as.numeric(cdf[1:t_ints_real[i], "Y"]))^2
  scm_rmse <- sqrt(sum(scm_error)/t_ints_real[i])
  
  cscm_error <- (as.numeric(bdf[1:t_ints_real[i], "Y"]) - as.numeric(cdf[1:t_ints_real[i], "Y"]))^2
  cscm_rmse <- sqrt(sum(cscm_error)/t_ints_real[i])
  
  error_df <- data.frame(treatment_id = i,
                         scm_rmse = scm_rmse,
                         cscm_rmse = cscm_rmse
                         )
  
  adf$Method <- "SCM"
  bdf$Method <- "CSCM"
  
  cfdf <- rbind(adf, bdf)
  cfdf <- cfdf %>% mutate(Y = as.numeric(as.character(Y)),
                          Quarter = as.numeric(as.character(Quarter)))
  cdf <- cdf %>% mutate(Y = as.numeric(as.character(Y)),
                        Quarter = as.numeric(as.character(Quarter)))
  
  cfdf_x = cfdf$Quarter-t_ints_real[i]
  cfdf_y = cfdf$Y
  cfdf_method = cfdf$Method
  cdf_y = cdf$Y
  
  case_data_list <- data.frame(
    cfdf_x = cfdf$Quarter-t_ints_real[i], 
    cfdf_y = cfdf$Y, 
    cfdf_method = cfdf$Method,
    cdf_y = cdf$Y
  )
  
  # Second method of IE calculation (sum)
  # df <- as.data.frame(subset(case_data_list, cfdf_x >= 0 & cfdf_x < 8))   # for in-time placebo
  df <- as.data.frame(subset(case_data_list, cfdf_x >= 0))                  # main analysis and in-space placebo 
  colnames(df) = c("time", "counterfactual", "method", "cases")
  df$time <- as.integer(df$time)
  
  df <- df %>%
    filter(time >= 0)
  
  df_new <- df %>%
    spread(key = method, value = counterfactual)
  
  df_new$cases_cumsum <- cumsum(df_new$cases)
  df_new$CSCM_cumsum <- cumsum(df_new$CSCM)
  df_new$SCM_cumsum <- cumsum(df_new$SCM)
  
  # Calculate the intervention efficacy
  df_new$IE_CSCM <- (df_new$CSCM_cumsum - df_new$cases_cumsum)*100 / df_new$CSCM_cumsum
  df_new$IE_SCM <- (df_new$SCM_cumsum - df_new$cases_cumsum)*100 / df_new$SCM_cumsum
  
  return(list(error_df, case_data_list, df_new))
}

## treatment_ids <- c(1:32)
t_ints_real <- c(rep(42, 3), rep(43, 11), rep(45, 13), rep(51, 5), rep(45, 19)) 

## for in-time placebo 
# treatment_ids <- c(1:32)
# t_ints_real <- c(rep(42-8, 3), rep(43-8, 11), rep(45-8, 13), rep(51-8, 5), rep(45-8, 19))

# for in-space placebo
treatment_ids <- c(33:51)

registerDoParallel(cores = 18)

#results_interventions 
brazil_results <- foreach(i = treatment_ids, .combine = "rbind", .export = c("countSynth", "csynth_inner", "factor.init", "factor.sim", 
                                                                             "holdout.csynth", "brazil.interventions"), 
                          .packages = c('dplyr', 'zoo', 'readxl', 'ggplot2', 'Synth', 'glmnet', 'dplyr', 
                                        'osqp', 'optimx', 'tidyr', 'lubridate', 'zoo', 'data.table')) %dopar% {
                                          
                                          source("/your/local/filepath/CSCM_helper_functions.R")
                                          tryCatch({
                                          sink("brazil.txt", append=TRUE)
                                          #change dataset
                                          test <- brazil.interventions(i, brazil.den, t_ints_real)
                                          },
                                          error=function(...) NULL)
                                          return(test)
                                        }

stopImplicitCluster()

brz_error_df <- do.call("rbind", brazil_results[,1])

write_xlsx(brz_error_df, "/your/local/filepath/wolb_interventions/transfer.xlsx")

# Get aggregated IE by event-time
brz.agg.dfs <- lapply(brazil_results[,3], function(df) {
  df$time <- df$time * 3
  breaks <- c(0, 6, 12, 18, 24, 30, 36)
  labels <- c("1-6 months", "7-12 months", "13-18 months", "19-24 months", ">24 months", ">24 months")
  df$event_time <- cut(df$time, breaks = breaks, labels = labels, include.lowest = TRUE)
  return(df)
})

brz.combined.ie <- do.call(rbind, brz.agg.dfs)
brz.grouped.ie <- split(brz.combined.ie, brz.combined.ie$event_time)

# Calculate the aggregate IE for each dataframe
brz.et <- lapply(brz.grouped.ie, function(group) {
  if (nrow(group) > 0) {
    cscm_sum <- sum(group$CSCM)
    scm_sum <- sum(group$SCM)
    cases_sum <- sum(group$cases)
    
    # Calculate IE
    IE_CSCM <- (cscm_sum - cases_sum)*100 / cscm_sum
    IE_SCM <- (scm_sum - cases_sum)*100 / scm_sum
    
    # Return a dataframe
    data.frame(event_time = unique(group$event_time), IE_CSCM = IE_CSCM, IE_SCM = IE_SCM)
  }
})

brz.et <- do.call(rbind, brz.et)

write_xlsx(brz.et, "/your/local/filepath/transfer.xlsx")

# New IE calculation extraction
brazil.sc.ie <- data.frame()

ie.extract <- lapply(brazil_results[,3], function(df) {
  df <- df[, 5:9]
  df <- df[nrow(df), ]
  return(df)
})

# Combine the rows into a new dataframe
brazil.sc.ie <- do.call(rbind, ie.extract)
brazil.sc.ie <- brazil.sc.ie[, c(1:3, 5, 4)]
write_xlsx(brazil.sc.ie, "/your/local/filepath/transfer.xlsx")

(sum(brazil.sc.ie$CSCM_cumsum) - sum(brazil.sc.ie$cases_cumsum))*100/sum(brazil.sc.ie$CSCM_cumsum)
(sum(brazil.sc.ie$SCM_cumsum) - sum(brazil.sc.ie$cases_cumsum))*100/sum(brazil.sc.ie$SCM_cumsum)

# Run confidence intervals first (brazil-confints.R)
library(dplyr)
agg.IE.df <- do.call("rbind", brazil_results[,3])

agg.IE.df <- agg.IE.df %>%
  group_by(IE_df_x, IE_df_method) %>%
  summarise(mean_value = mean(IE_df_y, na.rm = TRUE))

extract_row <- function(mat, row_number, max_cols) {
  row <- mat[row_number, ]
  length(row) <- max_cols  # Extend the length of row vector
  return(row)
}

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
