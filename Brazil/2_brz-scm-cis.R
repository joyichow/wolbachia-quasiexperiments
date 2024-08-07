source("/your/local/filepath/CSCM_helper_functions.R") # change to your local file path

libraries <- c(
  "dplyr", "lubridate", "meboot", "tidyverse", 
  "doParallel", "foreach", "Synth", "glmnet", "osqp"
)

for (lib in libraries) {
  if (!require(lib, character.only = TRUE)) {
    install.packages(lib, dependencies = TRUE)
    library(lib, character.only = TRUE)
  }
}

brazil.den <- read.csv("/your/local/filepath/Niteroi_NeighbourhoodData_ExclJurujuba_Apr2021.csv")
brazil.den <- brazil.den[c(1:2, 7)]

brazil.den <- brazil.den %>%
  pivot_wider(names_from = neighborhood, values_from = notified)

brazil.ref <- brazil.den

brazil.den <- as.ts(brazil.den)
brazil.den <- brazil.den[,-1]

set.seed(123)
(store <- meboot(brazil.den, trim = list(xmin = 0), reps = 1000))
check <- store$ensemble

bootstrap.dfs <- vector("list", 1000)

for (i in 1:ncol(check)) {
  time_series <- check[, i]
  reshaped_time_series <- matrix(time_series, nrow = 159, ncol = 51, byrow = FALSE)
  bootstrap.dfs[[i]] <- as.data.frame(reshaped_time_series)
  bootstrap.dfs[[i]] <- cbind(brazil.ref[,1], bootstrap.dfs[[i]])
  
  colnames(bootstrap.dfs[[i]]) <- colnames(brazil.ref)
}

process_data_frame <- function(df) {
  df <- df %>%
    pivot_longer(cols = -yearmonth, names_to = "neighborhood", values_to = "notified")
  
  df <- df %>%
    left_join(df %>%
                distinct(neighborhood) %>%
                mutate(site_ID = row_number()), 
              by = "neighborhood")
  
  df <- df %>%
    mutate(yearmonth = as.Date(paste0(yearmonth, "_01"), format="%Y_%m_%d"),
           quarter = paste(year(yearmonth), quarter(yearmonth), sep = "_"))
  
  df <- df %>%
    group_by(neighborhood, quarter, site_ID) %>%
    summarize(notified = sum(notified), .groups = 'drop')
  
  df <- df %>%
    mutate(time_ID = as.integer(factor(quarter, levels = unique(quarter))))
  
  df$site_ID <- as.numeric(df$site_ID)
  
  df <- as.data.frame(df)
  
  return(df)
}

bootstrap.dfs <- lapply(bootstrap.dfs, process_data_frame)

normalize <- function(df) {
  df <- df %>%
    group_by(neighborhood) %>%
    mutate(notified_norm = (notified - min(notified)) / (max(notified) - min(notified)))
  
  site_maxmins <- df %>%
    filter(site_ID %in% 1:51) %>%
    group_by(site_ID) %>%
    summarise(min_notified = min(notified), max_notified = max(notified), .groups = 'drop')
  
  list(normalized_df = as.data.frame(df), min_max_values = as.data.frame(site_maxmins))
}

normalized.res <- lapply(bootstrap.dfs, normalize)

bootstrap.dfs.normalized <- lapply(normalized.res, `[[`, "normalized_df")
min_max_values <- lapply(normalized.res, `[[`, "min_max_values")

brazil.interventions <- function(i, t_ints_real, bootstrap.dfs.normalized){
  
  col_names <- as.character(seq(from = 1, to = 53 - t_ints_real[i]))
  values <- lapply(col_names, function(x) numeric(0))
  
  IE_SCM <- data.frame(setNames(values, col_names))
  IE_CSCM <- data.frame(setNames(values, col_names))
  SCM_case <- data.frame(setNames(values, col_names))
  CSCM_case <- data.frame(setNames(values, col_names))
  obs_case <- data.frame(setNames(values, col_names))
  
  for (k in seq_along(bootstrap.dfs)){
    print(paste("current site:", i, "; bs:", k))
    current.df <- bootstrap.dfs.normalized[[k]]
    
    dengue.main.res <- countSynth(data=current.df,                                 
                                  predictors=NULL,
                                  dependent="notified_norm",
                                  unit.variable="site_ID",
                                  time.variable = "time_ID",
                                  treatment.identifier = i,
                                  controls.identifier = c(33:51),               
                                  t_int = t_ints_real[i],
                                  min_1se=F,
                                  K=2)
    
    synth.orig <- dengue.main.res$dataprep.main$Y0plot %*% as.numeric(dengue.main.res$unit.weight.SCM)
    cscm.cf <- dengue.main.res$dataprep.main$Y0plot %*% as.numeric(dengue.main.res$unit.weight.full.sample)
    year <- 1:53
    
    adf <- as.data.frame(cbind(year,synth.orig,rep("SCM",length(year))))
    bdf <- as.data.frame(cbind(year,cscm.cf,rep("CSCM",length(year))))
    cdf <- as.data.frame(cbind(year,dengue.main.res$dataprep.main$Y1plot,rep("Observed",length(year))))
    
    colnames(adf)=colnames(bdf)=colnames(cdf)=c("Year","Y","Method")
    
    max.norm <- min_max_values[[k]][i,3]
    min.norm <- min_max_values[[k]][i,2]
    
    adf$Y <- as.numeric(adf$Y)*(max.norm - min.norm) + min.norm
    bdf$Y <- as.numeric(bdf$Y)*(max.norm - min.norm) + min.norm
    cdf$Y <- as.numeric(cdf$Y)*(max.norm - min.norm) + min.norm
    
    SCM_case <- rbind(SCM_case, adf$Y)
    CSCM_case <- rbind(CSCM_case, bdf$Y)
    obs_case <- rbind(obs_case, cdf$Y)
    
    IE_scm_data <- data.frame((as.numeric(adf[t_ints_real[i]:53, "Y"]) - as.numeric(cdf[t_ints_real[i]:53, "Y"]))*100/ as.numeric(adf[t_ints_real[i]:53, "Y"]))
    
    IE_cscm_data <- data.frame((as.numeric(bdf[t_ints_real[i]:53, "Y"]) - as.numeric(cdf[t_ints_real[i]:53, "Y"]))*100/ as.numeric(bdf[t_ints_real[i]:53, "Y"]))
    
    IE_scm_data$Method <- 'SCM'
    IE_cscm_data$Method <- 'CSCM'
    
    IE_scm_data$Quarter <- c(t_ints_real[i]:53)
    IE_cscm_data$Quarter <- c(t_ints_real[i]:53)
    
    colnames(IE_cscm_data)=colnames(IE_scm_data)=c("intervention_efficacy", "Method", "Quarter")
    
    IE_SCM <- rbind(IE_SCM, (IE_scm_data$intervention_efficacy))
    IE_CSCM <- rbind(IE_CSCM, (IE_cscm_data$intervention_efficacy))
  }
  return(list(IE_SCM, IE_CSCM, SCM_case, CSCM_case, obs_case))
}

treatment_ids <- c(1:32)

# treatment_ids <- c(33:51) # in space placebo
t_ints_real <- c(rep(42, 3), rep(43, 11), rep(45, 13), rep(51, 5), rep(45, 19))       #treatment + in-space placebo

# t_ints_real <- c(rep(42-5, 3), rep(43-5, 11), rep(45-5, 13), rep(51-5, 5), rep(45, 19)) #in-time placebo

registerDoParallel(cores = 20)

#results_interventions 

brazil_confints.agg <- foreach(i = treatment_ids, .combine = "rbind", .export = c("countSynth", "csynth_inner", "factor.init", "factor.sim", 
                                                                                  "holdout.csynth", "brazil.interventions"), 
                               .packages = c('dplyr', 'zoo', 'readxl', 'ggplot2', 'Synth', 'glmnet', 'dplyr', 
                                             'osqp', 'optimx', 'tidyr', 'lubridate', 'zoo', 'data.table')) %dopar% {
                                               
                                               source("/Users/limlab/Desktop/JYC - wolbachia/CSCM_helper_functions.R")
                                               # tryCatch({
                                               sink("brazil.txt", append=TRUE)
                                               #change dataset
                                               test <- brazil.interventions(i, t_ints_real, bootstrap.dfs.normalized)
                                               # }, 
                                               # error=function(...) NULL)
                                               return(test)
                                             }

stopImplicitCluster()

# Get confidence intervals for each site 
IE.vals <- data.frame()
brazil.agg.scm <- data.frame(matrix(ncol = 0, nrow = 1000))
brazil.agg.cscm <- data.frame(matrix(ncol = 0, nrow = 1000))
brazil.agg.obs <- data.frame(matrix(ncol = 0, nrow = 1000))

for (i in 1:nrow(brazil_confints.agg)) {
  time <- t_ints_real[i]
  scm.df <- data.frame(brazil_confints.agg[i,3])
  cscm.df <- data.frame(brazil_confints.agg[i,4])
  obs.df <- data.frame(brazil_confints.agg[i,5])
  
  colnames(scm.df) = colnames(cscm.df) = colnames(obs.df) = c(1:53)

  scm.agg <- rowSums(scm.df[, as.numeric(colnames(scm.df)) >= time])
  cscm.agg <- rowSums(cscm.df[, as.numeric(colnames(cscm.df)) >= time])
  obs.agg <- rowSums(obs.df[, as.numeric(colnames(obs.df)) >= time])
  
  brazil.agg.scm <- cbind(brazil.agg.scm, scm.agg)
  brazil.agg.cscm <- cbind(brazil.agg.cscm, cscm.agg)
  brazil.agg.obs <- cbind(brazil.agg.obs, obs.agg)
  
  scm.ie <- (scm.agg - obs.agg)*100/scm.agg
  cscm.ie <- (cscm.agg - obs.agg)*100/cscm.agg

  scm.confint <- quantile(scm.ie, c(0.025, 0.975))
  cscm.confint <- quantile(cscm.ie, c(0.025, 0.975))
  
  IE.vals <- rbind(IE.vals, data.frame(scm.confint[[1]], scm.confint[[2]], cscm.confint[[1]], cscm.confint[[2]]))
}

colnames(IE.vals) = c("scm_lower", "scm_upper", "cscm_lower", "cscm_upper")

scm.site.ie <- sprintf("%.2f (%.2f-%.2f)", brazil.sc.ie[,4], IE.vals[,1], IE.vals[,2])
cscm.site.ie <- sprintf("%.2f (%.2f-%.2f)", brazil.sc.ie[,5], IE.vals[,3], IE.vals[,4])

site.ie <- data.frame(cbind(scm.site.ie, cscm.site.ie))
write_xlsx(site.ie, "/your/local/filepath/transfer.xlsx")

# Aggregate calculation
agg.scm.ie <- (rowSums(brazil.agg.scm) - rowSums(brazil.agg.obs))*100/rowSums(brazil.agg.scm)
agg.cscm.ie <- (rowSums(brazil.agg.cscm) - rowSums(brazil.agg.obs))*100/rowSums(brazil.agg.cscm)

agg.scm.confint <- quantile(agg.scm.ie, c(0.025, 0.975))
agg.cscm.confint <- quantile(agg.cscm.ie, c(0.025, 0.975))

# Get aggregate by time period
calculate_efficacy <- function(predicted_list, observed_list, start_times) {
  timeframes <- list(0:2, 3:4, 5:6, 7:8, 9:11)
  results <- list()
  
  for (i in seq_along(timeframes)) {
    timeframe <- timeframes[[i]]
    efficacy_samples <- list()
    
    for (sample in 1:nrow(predicted_list[[1]])) {
      sum_predicted <- 0
      sum_observed <- 0
      
      for (site in seq_along(predicted_list)) {
        start_quarter <- start_times[site]
        
        for (quarter in timeframe) {
          if ((start_quarter + quarter) <= ncol(predicted_list[[site]])) {
            sum_predicted <- sum_predicted + predicted_list[[site]][sample, start_quarter + quarter]
            sum_observed <- sum_observed + observed_list[[site]][sample, start_quarter + quarter]
          }
        }
      }
      
      ie <- (sum_predicted - sum_observed) * 100 / sum_predicted
      efficacy_samples[[sample]] <- ie
    }
    
    results[[i]] <- efficacy_samples
  }
  
  return(results)
}

scm.ie.et <- calculate_efficacy(predicted_list = brazil_confints.agg[,3], 
                                observed_list = brazil_confints.agg[,5], 
                                start_times = t_ints_real)

cscm.ie.et <- calculate_efficacy(predicted_list = brazil_confints.agg[,4], 
                                observed_list = brazil_confints.agg[,5], 
                                start_times = t_ints_real)

# Get aggregate confidence intervals 
sum_bootstrap_iterations <- function(df_list) {
  n_iterations <- nrow(df_list[[1]])
  
  sum_list <- vector("list", n_iterations)
  for (i in 1:n_iterations) {
    sum_list[[i]] <- lapply(df_list, function(df) df[i, 42:53]) %>%
      Reduce('+', .)
  }
  
  sum_matrix <- do.call(rbind, sum_list)
  
  return(sum_matrix)
}

aggregate_SCM <- sum_bootstrap_iterations(brazil_confints.agg[,3])
aggregate_CSCM <- sum_bootstrap_iterations(brazil_confints.agg[,4])
aggregate_obs <- sum_bootstrap_iterations(brazil_confints.agg[,5])

total.scm <- rowSums(aggregate_SCM)
total.cscm <- rowSums(aggregate_CSCM)
total.obs <- rowSums(aggregate_obs)

scm.ie <- (total.scm - total.obs)*100/total.scm
cscm.ie <- (total.cscm - total.obs)*100/total.cscm

SCM_timepoints <- timepoint.CI (brazil_confints.agg[,3])
CSCM_timepoints <- timepoint.CI (brazil_confints.agg[,4])
obs_timepoints <- timepoint.CI (brazil_confints.agg[,5])

# IE CIs for each intervention site 
calculate.CI <- function(data_list) {
  lapply(data_list, function(df) {
    apply(df, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
  })
}

confidence_intervals_SCM_normalized <- calculate.CI(brazil_confints[,1])
confidence_intervals_CSCM_normalized <- calculate.CI(brazil_confints[,2])

conf.int.avg.brz <- data.frame()

for (k in 1:32){
  current <- paste0("result.",k)
  scm_df <- rowMeans(confidence_intervals_SCM_normalized[[current]])
  upper_scm <- scm_df[2]
  lower_scm <- scm_df[1]
  
  cscm_df <- rowMeans(confidence_intervals_CSCM_normalized[[current]])
  upper_cscm <- cscm_df[2]
  lower_cscm <- cscm_df[1]
  
  conf.int.avg.brz <- rbind(conf.int.avg.brz, data.frame(lower_scm, upper_scm, lower_cscm, upper_cscm))
}


library(writexl)
write_xlsx(conf.int.avg.brz, "/your/local/filepath/confint_brz_scm.xlsx")
write_xlsx(as.data.frame(confidence_intervals_CSCM_normalized), "/your/local/filepath/confint_brz_scm.xlsx")