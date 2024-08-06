library(stats)
library(meboot)
library(readxl)
library(tidyr)
library(lubridate)
library(zoo)
library(data.table)
library(doParallel)
library(foreach)

mal <- read_xlsx(" Dengue Incidence.xlsx") # incidence = (total cfm dengue cases per total pop)*100000 (weekly incidence per 100000
mal <- head(mal, 332)
mal <- mal[-c(1:2)]

mal <- as.ts(mal)

mal_ref <- read_xlsx(" Dengue Incidence.xlsx")
mal_ref <- head(mal_ref , 332)

store <- meboot(mal, trim = list(xmin=0), reps = 1000)
check <- store$ensemble

bootstrap.dfs <- vector("list", 1000)

process_df <- function(df) {
  df$Dates <- seq(as.Date("2013-01-02"), as.Date("2019-05-12"), by = 7)
  df$Quarter <- as.yearqtr(df$Dates)
  df$Y <- as.numeric(df$Y)
  df_sum <- df %>%
    group_by(Quarter) %>%
    summarise(Y = sum(Y, na.rm = TRUE))
  df_sum$Quarter <- factor(df_sum$Quarter, levels = unique(df_sum$Quarter))
  df_sum$Quarter <- as.integer(df_sum$Quarter)
  return(df_sum)
}

for (i in 1:ncol(check)) {
  time_series <- check[, i]
  reshaped_time_series <- matrix(time_series, nrow = 332, ncol = 25, byrow = FALSE)
  bootstrap.dfs[[i]] <- as.data.frame(reshaped_time_series)
  bootstrap.dfs[[i]] <- cbind(mal_ref[,1:2], bootstrap.dfs[[i]])
  colnames(bootstrap.dfs[[i]]) <- colnames(mal_ref)
}

process_data_frame <- function(mal) {
  mal$Time <- seq(1900,length.out=nrow(mal))
  
  long_mal <- mal %>%
    pivot_longer(
      cols = -c('Year', 'Epid Week'), 
      names_to = "intervention_site", 
      values_to = "cases"
    )
  
  long_mal <- data.table(long_mal)
  
  mal_qtr <- long_mal[, time_ID := .GRP, by = .(Year, `Epid Week`)]
  
  mal_qtr$site_ID <- mal_qtr$intervention_site
  
  identifiers <- c('MH','SF','SL','SC','AL','AF','M1','M2','M3','M4','M5','M6','M7','M8','S1','S2','S3','S4', "C1", "C2" ,"C3" ,"C4", "A1" ,"A2" ,"A3")
  mal_qtr$site_ID <- as.numeric(factor(mal_qtr$site_ID, levels = identifiers))
  
  mal_qtr <- as.data.frame(mal_qtr)
  
  return(mal_qtr)
}

bootstrap.dfs <- lapply(bootstrap.dfs, process_data_frame)

bootstrapping.syn.controls <- function(i, t_ints_real, x_ints_real, bootstrap.dfs){
  col_names <- as.character(seq(from = 1, to = 26 - x_ints_real[1]))
  values <- lapply(col_names, function(x) numeric(0))
  
  IE_SCM <- data.frame(setNames(values, col_names))
  IE_CSCM <- data.frame(setNames(values, col_names))
  SCM_case <- data.frame(setNames(values, col_names))
  CSCM_case <- data.frame(setNames(values, col_names))
  obs_case <- data.frame(setNames(values, col_names))
  
  for (k in seq_along(bootstrap.dfs)){
    print(paste("current site:", i, "current bootstrap:", k))
    tryCatch({
      mal_qtr <- bootstrap.dfs[[k]]
      
      dengue.main.res <- countSynth(data=mal_qtr,
                                    predictors=NULL,
                                    dependent="cases",
                                    unit.variable="site_ID",
                                    time.variable = "time_ID",
                                    treatment.identifier = i,
                                    controls.identifier = c(7:25),
                                    t_int=t_ints_real[i],
                                    min_1se=F,
                                    K=2)
      
      synth.orig <- dengue.main.res$dataprep.main$Y0plot %*% as.numeric(dengue.main.res$unit.weight.SCM)
      cscm.cf <- dengue.main.res$dataprep.main$Y0plot %*% as.numeric(dengue.main.res$unit.weight.full.sample)
      year <- 1:332
      
      adf <- as.data.frame(cbind(year,synth.orig,rep("SCM",length(year)), mal_ref$Year, mal_ref$`Epid Week`))
      bdf <- as.data.frame(cbind(year,cscm.cf,rep("CSCM",length(year)), mal_ref$Year, mal_ref$`Epid Week`))
      cdf <- as.data.frame(cbind(year,dengue.main.res$dataprep.main$Y1plot,rep("Observed",length(year)), mal_ref$Year, mal_ref$`Epid Week`))
      
      colnames(adf)=colnames(bdf)=colnames(cdf)=c("Year","Y","Method", "Labelled_Year", 'Epid_Week')
      
      scm_error <- (as.numeric(adf[1:t_ints_real[i], "Y"]) - as.numeric(cdf[1:t_ints_real[i], "Y"]))^2
      scm_rmse <- sqrt(sum(scm_error)/t_ints_real[i])
      
      cscm_error <- (as.numeric(bdf[1:t_ints_real[i], "Y"]) - as.numeric(cdf[1:t_ints_real[i], "Y"]))^2
      cscm_rmse <- sqrt(sum(cscm_error)/t_ints_real[i])
      
      adf <- process_df(adf)
      bdf <- process_df(bdf)
      cdf <- process_df(cdf)
      
      adf <- as.data.frame(adf)
      bdf <- as.data.frame(bdf)
      cdf <- as.data.frame(cdf)
      
      SCM_case <- rbind(SCM_case, adf$Y)
      CSCM_case <- rbind(CSCM_case, bdf$Y)
      obs_case <- rbind(obs_case, cdf$Y)
    }, error=function(...) NULL)
  }
  return(list(SCM_case, CSCM_case, obs_case))
}

treatment_ids <- c(1:6)
t_ints_real <- c(250,227,227,255,221,245)
x_ints_real <- c(20,18,18,20,17,19)

# cl <- makeCluster(cores )
registerDoParallel(cores=8)

bootstrapped.CI.agg <- foreach(i = treatment_ids, .combine = "rbind", .export = c("bootstrapping.syn.controls", "process_df", "process_data_frame"), 
                            
                            .packages = c('dplyr', 'zoo', 'readxl', 'ggplot2', 'Synth', 'glmnet', 'dplyr', 
                                          'osqp', 'optimx', 'tidyr', 'lubridate', 'zoo', 'data.table', 'stats', 'meboot')) %dopar% {
                                                 source("/Users/joyichow/Desktop/HPC_Server/CSCM_helper_functions.R")
                                                 tryCatch({
                                                   sink("log.txt", append=TRUE)
                                                   CI_IE <- bootstrapping.syn.controls(i, t_ints_real, x_ints_real, bootstrap.dfs)
                                                   
                                                 }, 
                                                 error=function(...) NULL)
                                                 return(CI_IE)
                                          }

stopImplicitCluster()

# New way of calculating confidence intervals

# Get confidence intervals for each event time 
calculate_efficacy <- function(predicted_list, observed_list, start_times) {
  timeframes <- list(0:2, 3:4, 5:6, 7:8)
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
    lower <- quantile(unlist(efficacy_samples), 0.025)
    upper <- quantile(unlist(efficacy_samples), 0.975)
    results[[i]] <- list("Lower" = lower, "Upper" = upper)
  }
  return(results)
}

scm.mal <- do.call(rbind, calculate_efficacy(predicted_list = bootstrapped.CI.agg[,1],
                                             observed_list = bootstrapped.CI.agg[,3],
                                             start_times = x_ints_real))

cscm.mal <- do.call(rbind, calculate_efficacy(predicted_list = bootstrapped.CI.agg[,2],
                                              observed_list = bootstrapped.CI.agg[,3],
                                              start_times = x_ints_real))

# Calculate site specific CIs
IE.vals <- data.frame(matrix(ncol = 0, nrow = 1000))
mal.agg.scm <- data.frame(matrix(ncol = 0, nrow = 1000))
mal.agg.cscm <- data.frame(matrix(ncol = 0, nrow = 1000))
mal.agg.obs <- data.frame(matrix(ncol = 0, nrow = 1000))

for (i in 1:nrow(bootstrapped.CI.agg)){
  time <- x_ints_real[i]
  scm.df <- data.frame(bootstrapped.CI.agg[i,1])
  cscm.df <- data.frame(bootstrapped.CI.agg[i,2])
  obs.df <- data.frame(bootstrapped.CI.agg[i,3])
  
  colnames(scm.df) = colnames(cscm.df) = colnames(obs.df) = c(1:26)
  
  scm.agg <- rowSums(scm.df[, as.numeric(colnames(scm.df)) >= time])
  cscm.agg <- rowSums(cscm.df[, as.numeric(colnames(cscm.df)) >= time])
  obs.agg <- rowSums(obs.df[, as.numeric(colnames(obs.df)) >= time])
  
  mal.agg.scm <- cbind(mal.agg.scm, scm.agg)
  mal.agg.cscm <- cbind(mal.agg.cscm, cscm.agg)
  mal.agg.obs <- cbind(mal.agg.obs, obs.agg)
  
  scm.ie <- (scm.agg - obs.agg)*100/scm.agg
  cscm.ie <- (cscm.agg - obs.agg)*100/cscm.agg
  
  scm.confint <- quantile(scm.ie, c(0.025, 0.975))
  cscm.confint <- quantile(cscm.ie, c(0.025, 0.975))
  
  IE.vals <- rbind(IE.vals, data.frame(scm.confint[[1]], scm.confint[[2]], cscm.confint[[1]], cscm.confint[[2]]))
}

colnames(IE.vals) = c("scm_lower", "scm_upper", "cscm_lower", "cscm_upper")

#Aggregate calculation for confidence interval
agg.scm.ie <- (rowSums(mal.agg.scm) - rowSums(mal.agg.obs))*100/rowSums(mal.agg.scm)
agg.cscm.ie <- (rowSums(mal.agg.cscm) - rowSums(mal.agg.obs))*100/rowSums(mal.agg.cscm)

agg.scm.confint <- quantile(agg.scm.ie, c(0.025, 0.975))
agg.cscm.confint <- quantile(agg.cscm.ie, c(0.025, 0.975))