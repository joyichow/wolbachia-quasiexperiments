brazil.interventions <- function(i, brazil.den, t_ints_real){
  print(paste("current site:", i))
  placebo_id <- c(33:51)
  
  dengue.main.res <- countSynth(data=brazil.den,                                 
                                predictors=NULL,
                                dependent="notified_norm",
                                unit.variable="site_ID",
                                time.variable = "time_ID",
                                treatment.identifier = i,
                                controls.identifier = c(33:51),                      # for actual treatment
                                # controls.identifier = placebo_id[placebo_id != i],  # for in-space placebo
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
  
  IE_scm <- mean((as.numeric(adf[t_ints_real[i]:53, "Y"]) - as.numeric(cdf[t_ints_real[i]:53, "Y"]))*100/ as.numeric(adf[t_ints_real[i]:53, "Y"]), na.rm = TRUE)
  
  IE_cscm <- mean((as.numeric(bdf[t_ints_real[i]:53, "Y"]) - as.numeric(cdf[t_ints_real[i]:53, "Y"]))*100/ as.numeric(bdf[t_ints_real[i]:53, "Y"]), na.rm = TRUE)
  
  IE_scm_6m <- mean((as.numeric(adf[(t_ints_real[i]+2):53, "Y"]) - as.numeric(cdf[(t_ints_real[i]+2):53, "Y"]))*100/ as.numeric(adf[(t_ints_real[i]+2):53, "Y"]), na.rm = TRUE)
  
  IE_cscm_6m <- mean((as.numeric(bdf[(t_ints_real[i]+2):53, "Y"]) - as.numeric(cdf[(t_ints_real[i]+2):53, "Y"]))*100/ as.numeric(bdf[(t_ints_real[i]+2):53, "Y"]), na.rm = TRUE)
  
  # Uncomment for in-time placebo
  # scm_error <- (as.numeric(adf[1:t_ints_real[i], "Y"]) - as.numeric(cdf[1:t_ints_real[i], "Y"]))^2
  # scm_rmse <- sqrt(sum(scm_error)/t_ints_real[i])
  # 
  # cscm_error <- (as.numeric(bdf[1:t_ints_real[i], "Y"]) - as.numeric(cdf[1:t_ints_real[i], "Y"]))^2
  # cscm_rmse <- sqrt(sum(cscm_error)/t_ints_real[i])
  # 
  # IE_scm <- mean((as.numeric(adf[t_ints_real[i]:(t_ints_real[i]+4), "Y"]) - as.numeric(cdf[t_ints_real[i]:(t_ints_real[i]+4), "Y"]))*100/ as.numeric(adf[t_ints_real[i]:(t_ints_real[i]+4), "Y"]))
  # 
  # IE_cscm <- mean((as.numeric(bdf[t_ints_real[i]:(t_ints_real[i]+4), "Y"]) - as.numeric(cdf[t_ints_real[i]:(t_ints_real[i]+4), "Y"]))*100/ as.numeric(bdf[t_ints_real[i]:(t_ints_real[i]+4), "Y"]))
  # 
  # IE_scm_6m <- mean((as.numeric(adf[(t_ints_real[i]+2):(t_ints_real[i]+4), "Y"]) - as.numeric(cdf[(t_ints_real[i]+2):(t_ints_real[i]+4), "Y"]))*100/ as.numeric(adf[(t_ints_real[i]+2):(t_ints_real[i]+4), "Y"]))
  # 
  # IE_cscm_6m <- mean((as.numeric(bdf[(t_ints_real[i]+2):(t_ints_real[i]+4), "Y"]) - as.numeric(cdf[(t_ints_real[i]+2):(t_ints_real[i]+4), "Y"]))*100/ as.numeric(bdf[(t_ints_real[i]+2):(t_ints_real[i]+4), "Y"]))
  # 
  # print("checkpoint D")
  
  error_df <- data.frame(treatment_id = i,
                         scm_rmse = scm_rmse,
                         cscm_rmse = cscm_rmse,
                         IE_scm = IE_scm,
                         IE_cscm = IE_cscm,
                         IE_scm_6m = IE_scm_6m,
                         IE_cscm_6m = IE_cscm_6m)
  
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
  
  IE_scm_data <- data.frame((as.numeric(adf[t_ints_real[i]:53, "Y"]) - as.numeric(cdf[t_ints_real[i]:53, "Y"]))*100/ as.numeric(adf[t_ints_real[i]:53, "Y"]))
  
  IE_cscm_data <- data.frame((as.numeric(bdf[t_ints_real[i]:53, "Y"]) - as.numeric(cdf[t_ints_real[i]:53, "Y"]))*100/ as.numeric(bdf[t_ints_real[i]:53, "Y"]))
  
  IE_scm_data$Method <- 'SCM'
  IE_cscm_data$Method <- 'CSCM'
  
  IE_scm_data$Quarter <- c(t_ints_real[i]:53)
  IE_cscm_data$Quarter <- c(t_ints_real[i]:53)
  
  colnames(IE_cscm_data)=colnames(IE_scm_data)=c("intervention_efficacy", "Method", "Quarter")
  
  IE_df <- rbind(IE_scm_data, IE_cscm_data)
  IE_df <- IE_df %>% mutate(intervention_efficacy = as.numeric(as.character(intervention_efficacy)),
                            Quarter = as.numeric(as.character(Quarter)))
  
  IE_data_list <- data.frame(
    treatment_id = i,
    IE_df_x = IE_df$Quarter-t_ints_real[i],
    IE_df_y = IE_df$intervention_efficacy,
    IE_df_method = IE_df$Method
  )
  
  # Second method of IE calculation (sum)
  df <- as.data.frame(case_data_list)
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
  
  return(list(error_df, case_data_list, IE_data_list, df_new))
}

extract_row <- function(mat, row_number, max_cols) {
  row <- mat[row_number, ]
  length(row) <- max_cols  # Extend the length of row vector
  return(row)
}

# Formatting confidence intervals
format <- function(data) {
  formatted_string <- paste0(
    sprintf("%.2f", data[[1]]), 
    " (", 
    sprintf("%.2f", data[[2]]), 
    "-", 
    sprintf("%.2f", data[[3]]), 
    ")"
  )
  data.frame(formatted_string)
}

# Aggregate weekly data into quarterly data
aggregate_by_quarter <- function(df, i, x_ints_real) {
  agg_df <- aggregate(cbind(notified_norm*100, Est, pLL, pUL) ~ quarter, data = df, sum)
  
  colnames(agg_df) = c("quarter", "notified_norm", "Est", "pLL", "pUL")
  
  agg_df$time_ID <- c(1:53) - x_ints_real
  agg_df$site <- i
  
  agg_df <- subset(agg_df, time_ID >= 0)
  
  agg_df$cases.cs <- cumsum(agg_df$notified_norm)
  agg_df$est.cs <- cumsum(agg_df$Est)
  agg_df$pLL.cs <- cumsum(agg_df$pLL)
  agg_df$pUL.cs <- cumsum(agg_df$pUL)
  
  return(list(agg_df))
}

# Calculate the intervention efficacies 
calculate_ie <- function(df, est_col, case_col, pLL_col, pUL_col) {
  IE <- (sum(df[[est_col]]) - sum(df[[case_col]])) * 100 / sum(df[[est_col]])
  lIE <- (sum(df[[pLL_col]]) - sum(df[[case_col]])) * 100 / sum(df[[pLL_col]])
  uIE <- (sum(df[[pUL_col]]) - sum(df[[case_col]])) * 100 / sum(df[[pUL_col]])
  
  list(IE = IE, lIE = lIE, uIE = uIE)
}

# Calculate intervention efficacies by event time
calculate.et <- function(df) {
  df$time <- df$time * 3
  breaks <- c(0, 6, 12, 18, 24, 30)
  labels <- c("1-6 months", "7-12 months", "13-18 months", "19-24 months", ">24 months")
  df$event_time <- cut(df$time, breaks = breaks, labels = labels, include.lowest = TRUE)
  
  df_grouped <- split(df, df$event_time)
  
  agg_df <- lapply(df_grouped, function(group) {
    if (nrow(group) > 0) {
      est_sum <- sum(group$Est)
      pLL_sum <- sum(group$pLL)
      pUL_sum <- sum(group$pUL)
      cases_sum <- sum(group$notified_norm)
      
      IE <- (est_sum - cases_sum)*100 / est_sum
      l.IE <- (pLL_sum - cases_sum)*100 / pLL_sum
      u.IE <- (pUL_sum - cases_sum)*100 / pUL_sum
      
      data.frame(event_time = unique(group$event_time), IE = IE, lower = l.IE, upper = u.IE)
    }
  })
  
  agg_df <- do.call(rbind, agg_df)
  return(agg_df)
}