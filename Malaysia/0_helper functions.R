# helper functions 

# Aggregate weekly data to quarterly to reduce noise 
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

# Formatting for confidence intervals
format_interval <- function(central, lower, upper) {
  sprintf("%.2f (%.2f-%.2f)", central, lower, upper)
}

# Main function for synthetic control
syn.control.interventions <- function(i, mal_qtr, treatment_ids, x_ints_real, t_ints_real){
  # placebo.donour <- treatment_ids[treatment_ids != i]                         # for in-space placebo
  
  dengue.main.res <- countSynth(data=mal_qtr,
                                predictors=NULL,
                                dependent="cases",
                                unit.variable="site_ID",
                                time.variable = "time_ID",
                                treatment.identifier = i,
                                # controls.identifier = placebo.donour,         # for in-space placebo
                                controls.identifier = c(7:25),                  # for in-time placebo
                                t_int=t_ints_real[i],
                                min_1se=F,
                                K=2)
  
  synth.orig <- dengue.main.res$dataprep.main$Y0plot %*% as.numeric(dengue.main.res$unit.weight.SCM)
  cscm.cf <- dengue.main.res$dataprep.main$Y0plot %*% as.numeric(dengue.main.res$unit.weight.full.sample)
  year <- 1:332
  
  adf <- as.data.frame(cbind(year,synth.orig,rep("SCM",length(year)), mal$Year, mal$`Epid Week`))
  bdf <- as.data.frame(cbind(year,cscm.cf,rep("CSCM",length(year)), mal$Year, mal$`Epid Week`))
  cdf <- as.data.frame(cbind(year,dengue.main.res$dataprep.main$Y1plot,rep("Observed",length(year)), mal$Year, mal$`Epid Week`))
  
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
  
  error_df <- data.frame(treatment_id = i,
                         scm_rmse = scm_rmse,
                         cscm_rmse = cscm_rmse)
  
  adf$Method <- 'SCM'
  bdf$Method <- 'CSCM'
  
  cfdf <- rbind(adf, bdf)
  cfdf <- cfdf %>% mutate(Y = as.numeric(as.character(Y)),
                          Quarter = as.numeric(as.character(Quarter)))
  cdf <- cdf %>% mutate(Y = as.numeric(as.character(Y)),
                        Quarter = as.numeric(as.character(Quarter)))
  
  cfdf_x = cfdf$Quarter-x_ints_real[i]
  cfdf_y = cfdf$Y
  cfdf_method = cfdf$Method
  cdf_y = cdf$Y
  
  case_data_list <- data.frame(
    cfdf_x = cfdf$Quarter-x_ints_real[i], 
    cfdf_y = cfdf$Y, 
    cfdf_method = cfdf$Method,
    cdf_y = cdf$Y
  )
  
  # Second method of IE calculation (sum)
  #  df <- as.data.frame(subset(case_data_list, cfdf_x >= 0 & cfdf_x < 8))        # for in-time placebo
  df <- as.data.frame(subset(case_data_list, cfdf_x >= 0))                        # for main analysis and in-space placebo
  colnames(df) = c("time", "counterfactual", "method", "cases")
  df$time <- as.integer(df$time)
  
  df_new <- df %>% 
    spread(key = method, value = counterfactual)
  
  df_new$cases_cumsum <- cumsum(df_new$cases)
  df_new$CSCM_cumsum <- cumsum(df_new$CSCM)
  df_new$SCM_cumsum <- cumsum(df_new$SCM)
  
  df_new$IE_CSCM <- (df_new$CSCM_cumsum - df_new$cases_cumsum)*100 / df_new$CSCM_cumsum
  df_new$IE_SCM <- (df_new$SCM_cumsum - df_new$cases_cumsum)*100 / df_new$SCM_cumsum
  
  return(list(error_df, case_data_list, df_new))
}