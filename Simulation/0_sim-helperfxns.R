simulation.scm <- function(i, long.sim, treatment_ids, t_ints_real){
  
  dengue.main.res <- countSynth(data=long.sim,
                                predictors=NULL,
                                dependent="cases",
                                unit.variable="intervention_site",
                                time.variable = "time_ID",
                                treatment.identifier = i,
                                controls.identifier = c(3:size+2),                  
                                t_int=t_ints_real[i],
                                min_1se=F,
                                K=2)
  
  synth.orig <- dengue.main.res$dataprep.main$Y0plot %*% as.numeric(dengue.main.res$unit.weight.SCM)
  cscm.cf <- dengue.main.res$dataprep.main$Y0plot %*% as.numeric(dengue.main.res$unit.weight.full.sample)
  year <- 1:851
  
  adf <- as.data.frame(cbind(year,synth.orig,rep("SCM",length(year))))
  bdf <- as.data.frame(cbind(year,cscm.cf,rep("CSCM",length(year))))
  cdf <- as.data.frame(cbind(year,dengue.main.res$dataprep.main$Y1plot,rep("Observed",length(year))))
  
  colnames(adf)=colnames(bdf)=colnames(cdf)=c("Year","Y","Method")
  
  scm_error <- (as.numeric(adf[1:t_ints_real[i], "Y"]) - as.numeric(cdf[1:t_ints_real[i], "Y"]))^2
  scm_rmse <- sqrt(sum(scm_error)/t_ints_real[i])
  
  cscm_error <- (as.numeric(bdf[1:t_ints_real[i], "Y"]) - as.numeric(cdf[1:t_ints_real[i], "Y"]))^2
  cscm_rmse <- sqrt(sum(cscm_error)/t_ints_real[i])
  
  cfdf <- rbind(adf, bdf)
  cfdf <- cfdf %>% mutate(Y = as.numeric(as.character(Y)),
                          Year = as.numeric(as.numeric(Year)))
  cdf <- cdf %>% mutate(Y = as.numeric(as.character(Y)),
                        Year = as.numeric(as.numeric(Year)))
  
  cfdf_x = cfdf$Year-t_ints_real[i]
  cfdf_y = cfdf$Y
  cfdf_method = cfdf$Method
  cdf_y = cdf$Y
  
  case_data_list <- data.frame(
    cfdf_x = cfdf$Year-t_ints_real[i], 
    cfdf_y = cfdf$Y, 
    cfdf_method = cfdf$Method,
    cdf_y = cdf$Y
  )
  
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
  
  return(df_new)
}


# Formatting function for results dataframe 
format <- function(data) {
  formatted_string <- paste0(
    sprintf("%.2f", data[[1]])
  )
  data.frame(formatted_string)
}

aggregate_by_quarter <- function(agg_df, i, x_ints_real) {
  
  agg_df$time_ID <- c(1:851) - x_ints_real
  agg_df$site <- i
  
  agg_df <- subset(agg_df, time_ID >= 0)
  
  agg_df$cases.cs <- cumsum(agg_df$cases)
  agg_df$est.cs <- cumsum(agg_df$yhat)
  
  return(list(agg_df))
}

calculate_ie <- function(df, est_col, case_col) {
  IE <- (sum(df[[est_col]]) - sum(df[[case_col]])) * 100 / sum(df[[est_col]])
  
  return(IE)
}


