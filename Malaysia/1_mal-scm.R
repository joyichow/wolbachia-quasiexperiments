# Depending on whether you're running validity checks or the main analysis, 
# please make the changes in the code on Lines 68, 76, 132, 150, 218 by commenting or uncommenting the appropriate lines 

# Please change file locations for helper scripts on Lines 6 and 173

source("/path/to/your/local/CSCM_helper_functions.R") # Change this to your local file path (Bonander - CSCM_helper_functions.R)

# Load/install required packages
packages <- c("ggplot2", "Synth", "glmnet", "dplyr", "osqp", "optimx", "readxl", 
              "tidyr", "lubridate", "zoo", "data.table", "doParallel", "foreach")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
sapply(packages, require, character.only=TRUE)

# Version check 
ifelse(.Machine$sizeof.pointer == 4, warning("It appears you are running a 32-bit version of R. Please switch to 64-bit to ensure that the optimization does not produce different results due to rounding errors."),
       print("64-bit version check: OK."))

# Load and preprocess data
mal <- read_xlsx(" Dengue Incidence.xlsx") # incidence = (total cfm dengue cases per total pop)*100000 (weekly incidence per 100000
mal <- head(mal, 332)

mal$Time <- seq(1900,length.out=332)

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

# For main analysis 
treatment_ids <- c(1:6)
t_ints_real <- c(250,227,227,255,221,245)
x_ints_real <- c(20,18,18,20,17,19)

# Uncomment for in-time placebo
# treatment_ids <- c(1:6)
# t_ints_real <- c(250-104,227-104,227-104,255-104,221-104,245-104)
# x_ints_real <- c(20-8,18-8,18-8,20-8,17-8,19-8)

# Uncomment for in-space placebo
# treatment_ids <- c(7:25)
# t_ints_real <- c(rep(221, 19))
# x_ints_real <- c(rep(17,19))

# Parallel processing setup
cl <- makeCluster(detectCores())
registerDoParallel(cl)

results_interventions <- foreach(i = treatment_ids, .combine = "rbind", .export = c("countSynth", "csynth_inner", "factor.init", "factor.sim", "holdout.csynth", "syn.control.interventions", "process_df"), 
                                 .packages = c('dplyr', 'zoo', 'readxl', 'ggplot2', 'Synth', 'glmnet', 'osqp', 'optimx', 'tidyr', 'lubridate', 'data.table')) %dopar% {
                                   source("/path/to/your/local/CSCM_helper_functions.R") # Change this to your local file path
                                   tryCatch({
                                     test <- syn.control.interventions(i, mal_qtr, treatment_ids, x_ints_real, t_ints_real)
                                   }, error = function(...) NULL)
                                   return(test)
                                 }

stopCluster(cl)

########################
## Individual site IE ##
########################

# Get IEs for each site
mal.sc.ie <- data.frame()

ie.extract <- lapply(results_interventions[,3], function(df) {
  df <- df[, 5:9]
  df <- df[nrow(df), ]
  return(df)
})

mal.sc.ie <- do.call(rbind, ie.extract)
mal.sc.ie <-mal.sc.ie[, c(1:3, 5, 4)]

# Aggregate calculation for sites 
agg.SCM <- (sum(mal.sc.ie $SCM_cumsum) - sum(mal.sc.ie $cases_cumsum))*100/sum(mal.sc.ie $SCM_cumsum)
agg.CSCM <- (sum(mal.sc.ie$CSCM_cumsum) - sum(mal.sc.ie $cases_cumsum))*100/sum(mal.sc.ie $CSCM_cumsum)

# FOR MAIN ANALYSIS ONLY: Get the confidence intervals for each site's IE (run SCM_parallelised_CI.R first)
scm.site.ie <- format_interval(mal.sc.ie[,4], IE.vals[,1], IE.vals[,2])
cscm.site.ie <- format_interval(mal.sc.ie[,5], IE.vals[,3], IE.vals[,4])
site.ie <- data.frame(cbind(scm.site.ie, cscm.site.ie))

# FOR MAIN ANALYSIS ONLY : Aggregate IE confints
scm.agg.ie <- format_interval(agg.SCM, agg.scm.confint[,1], agg.scm.confint[,2])
cscm.agg.ie <- format_interval(agg.CSCM, agg.cscm.confint[,1], agg.cscm.confint[,2])

site.ie <- (rbind(site.ie, data.frame(scm.site.ie = scm.agg.ie, 
                                      cscm.site.ie = cscm.agg.ie)))

write_xlsx(site.ie, "/path/to/your/local/transfer.xlsx") # Change this to your local file path

# FOR IN-TIME AND IN-SPACE PLACEBO, exporting the IEs
scm.site.ie <- sprintf("%.2f", mal.sc.ie[,4])
cscm.site.ie <- sprintf("%.2f", mal.sc.ie[,5])
site.ie <- data.frame(cbind(scm.site.ie, cscm.site.ie))

scm.agg.ie <- sprintf("%.2f", agg.SCM)
cscm.agg.ie <- sprintf("%.2f", agg.CSCM)

site.ie <- (rbind(site.ie, data.frame(scm.site.ie = scm.agg.ie,
                                      cscm.site.ie = cscm.agg.ie)))

########################
### IE by event time ###
########################
mal.agg.dfs <- lapply(results_interventions[,3], function(df) {
  df$time <- df$time * 3
  breaks <- c(0, 6, 12, 18, 24, 30, 36)
  labels <- c("1-6 months", "7-12 months", "13-18 months", "19-24 months", ">24 months", ">24 months")
  df$event_time <- cut(df$time, breaks = breaks, labels = labels, include.lowest = TRUE)
  return(df)
})

mal.combined.ie <- do.call(rbind, mal.agg.dfs)
mal.grouped.ie <- split(mal.combined.ie, mal.combined.ie$event_time)

mal.et <- lapply(mal.grouped.ie, function(group) {
  if (nrow(group) > 0) {
    cscm_sum <- sum(group$CSCM)
    scm_sum <- sum(group$SCM)
    cases_sum <- sum(group$cases)
    
    # Calculate IE
    IE_CSCM <- (cscm_sum - cases_sum)*100 / cscm_sum
    IE_SCM <- (scm_sum - cases_sum)*100 / scm_sum
    data.frame(event_time = unique(group$event_time), IE_CSCM = IE_CSCM, IE_SCM = IE_SCM)
  }
})

mal.et <- do.call(rbind, mal.et)

# FOR MAIN-ANALYSIS ONLY: confidence intervals for event time (run SCM_parallelised_CI.R)
scm.site.ie <- format_interval(mal.et[,3], scm.mal[,1], scm.mal[,2])
cscm.site.ie <- format_interval(mal.et[,2], cscm.mal[,3], cscm.mal[,4])

scm.site.ie <- sprintf("%.2f", mal.et[,3])
mal.et[,2] <- scm.site.ie
mal.et[,3] <- cscm.site.ie

write_xlsx(mal.et, "/path/to/your/local/transfer.xlsx") # Change this to your local file path


