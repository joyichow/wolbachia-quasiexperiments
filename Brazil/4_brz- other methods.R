libraries <- c(
  "readxl", "tidyr", "ggplot2", "data.table", 
  "lubridate", "boot", "truncreg", "dplyr", "did"
)

# Load all libraries
for (lib in libraries) {
  if (!require(lib, character.only = TRUE)) {
    install.packages(lib, dependencies = TRUE)
    library(lib, character.only = TRUE)
  }
}

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

f <- function(model, data, i, newdata, population) {
  require(pscl)
  offset_vector <- rep(log(population), nrow(data))
  
  c.int <- model[["coefficients"]][["count"]][["(Intercept)"]]
  c.cut <- model[["coefficients"]][["count"]][["cutoff"]]
  z.int <- model[["coefficients"]][["zero"]][["(Intercept)"]]
  z.cut <- model[["coefficients"]][["zero"]][["cutoff"]]
  
  m <- zeroinfl(cases ~ cutoff + offset(population),
                data = data[i, ], dist = "negbin",
                start = list(count = c(c.int, c.cut), zero = c(z.int, z.cut)))
  mparams <- as.vector(t(do.call(rbind, coef(summary(m)))[, 1:2]))
  yhat <- predict(m, newdata, type = "response")
  return(c(mparams, yhat))
}

t_ints_real <- c(rep(42, 3), rep(43, 11), rep(45, 13), rep(51, 5), rep(45,19)) 

# for main analysis
treatment_ids <- c(1:32)

treatment_df <- data.frame(
  site_ID = treatment_ids,
  first.treat = c(rep(42, 3), rep(43, 11), rep(45, 13), rep(51, 5)),
  treat = c(rep(1,32))
)

# for in-space placebo 
treatment_ids <- c(33:51)

treatment_df <- data.frame(
  site_ID = treatment_ids,
  first.treat = c(rep(45,19)),
  treat = c(rep(1,19))
)

################################################################
##################### Pre-post analysis ########################
################################################################
interventions.df <- brazil.den[brazil.den$site_ID %in% treatment_ids,]

interventions.df <- left_join(interventions.df, treatment_df, by = "site_ID")

interventions.df  <- interventions.df  %>%
  mutate(rel.time = time_ID - first.treat)

interventions.df$cutoff <- ifelse(interventions.df$rel.time >= 0, 1, 0)

f <- function(model, data, i, newdata) {
  tryCatch({
    require(truncreg)
    
    c.int <- model[["coefficients"]][["(Intercept)"]]
    c.cut <- model[["coefficients"]][["cutoff"]]
    
    m <- truncreg(notified_norm*100 ~ cutoff,
                  data = data[i, ], point = 100, direction = "right", start = list(c.int, c.cut))
    
    mparams <- as.vector(t(do.call(rbind, as.data.frame(coef(summary(m))))[,1:2]))
    
    yhat <- predict(m, newdata, type = "response")
    return(c(mparams, yhat))
  }, error=function(...) NA)  # Returning NA in case of error
}

brz_pp <- data.frame()
pp.IE <- data.frame()
agg.pp <- data.frame()
count <- 1 

for (i in treatment_ids) {
  intervention.subs <- subset(interventions.df, interventions.df$site_ID == i)
  
  pp.model <- truncreg(notified_norm*100 ~ cutoff, data = intervention.subs,
                       point = 100, direction = "right")
  
  brz_pp_counterfactual <- intervention.subs
  brz_pp_counterfactual$cutoff <- 0
  
  set.seed(10)
  (res <- boot(intervention.subs, f, R = 1000, newdata = brz_pp_counterfactual, model = pp.model,
               parallel = "snow", ncpus = 8))
  
  yhat <- t(sapply((1:nrow(brz_pp_counterfactual)), function(i){
    out <- boot.ci(res, index = i, type =c("perc"))
    with(out, c(Est = (t0), pLL = (percent[4]), pUL = (percent[5])))
  } ))
  
  brz_pp_counterfactual <- cbind(brz_pp_counterfactual, yhat)
  brz_pp_counterfactual <- as.data.frame(brz_pp_counterfactual)
  brz_pp_agg <- as.data.frame(aggregate_by_quarter(brz_pp_counterfactual, i, t_ints_real[count]))
  
  last_row <- tail(brz_pp_agg, 1)
  
  IE <- calculate_ie(df = brz_pp_agg, est_col = 'Est', case_col = 'notified_norm', pLL_col = 'pLL', pUL_col = 'pUL')
  
  pp.IE <- rbind(pp.IE, data.frame(IE = IE[1],
                                   lower = IE[2],
                                   upper = IE[3]))
  
  brz_pp <- rbind(brz_pp, brz_pp_agg)
  agg.pp <- rbind(agg.pp, last_row)
  
  count = count + 1
}

# Aggregate IE calculation for pre-post
IE <- calculate_ie(df = agg.pp, est_col = 'est.cs', case_col = 'cases.cs', pLL_col = 'pLL.cs', pUL_col = 'pUL.cs')
pp.IE <- rbind(pp.IE, data.frame(IE = IE[[1]], lIE = IE[[2]], uIE = IE[[3]]))

pp.IE.final <- format(pp.IE)

# Aggregate IE calculation by event time for pre-post
agg_df <- calculate.et(brz_pp)
pp.agg_df <- format(agg_df[2:4])

##############################################################
################ Difference in differences ###################
##############################################################
# 2x2 DiD
intervention_times <- setNames(t_ints_real, treatment_ids)

brz_did <- brazil.den %>%
  mutate(time = ifelse(time_ID >= intervention_times[treatment_ids], 1, 0))

brz_did <- left_join(brazil.den, treatment_df, by = "site_ID")
brz_did$first.treat[is.na(brz_did$first.treat)] <- 0
brz_did$treat[is.na(brz_did$treat)] <- 0

f.did <- function(model, data, i, newdata) {
  tryCatch({
    require(truncreg)
  
  c.int <- model[["coefficients"]][["(Intercept)"]]
  c.treat <- model[["coefficients"]][["treat"]]
  c.time <- model[["coefficients"]][["time"]]
  c.did <- model[["coefficients"]][["did"]]
  
  m <- truncreg(notified_norm*100 ~ treat + time + did, data = data[i,],
                point = 100, direction = "right",
                start = list(c.int, c.treat, c.time, c.did))
  
  mparams <- as.vector(t(do.call(rbind, as.data.frame(coef(summary(m))))[,1:2]))
  
  yhat <- predict(m, newdata, type = "response")
  return(c(mparams, yhat))
  }, error=function(...) NA)
}

did.IE <- data.frame()
agg.did <- data.frame()
brz_did.total <- data.frame()
count <- 1

t_ints_real <- c(rep(42, 3), rep(43, 11), rep(45, 13), rep(51, 5), rep(45,19)) 

for (i in treatment_ids) {
  exclude_sites <- setdiff(treatment_ids, i)
  
  # For main analysis 
  # brz_did_subset <- brz_did[!brz_did$site_ID %in% exclude_sites, ]
  # 
  # brz_did_subset$time <- ifelse(brz_did_subset$time_ID >=t_ints_real[count], 1, 0)
  # brz_did_subset$did <- brz_did_subset$time * brz_did_subset$treat
  
  # For in-space placebo 
  brz_did_subset <- brz_did[brz_did$site_ID %in% exclude_sites, ]
  brz_did_subset$treat <- 0
  brz_did_subset <- rbind(brz_did_subset, data.frame(subset(brazil.den, site_ID == i), time = 0, treat = 1))
  brz_did_subset$time <- ifelse(brz_did_subset$time_ID >=t_ints_real[counter], 1, 0)
  brz_did_subset$did <- brz_did_subset$time * brz_did_subset$treat
  
  didreg <- truncreg(notified_norm*100 ~ treat + time + did, data = brz_did_subset,
                     point = 100, direction = "right")
  
  brz_did_counterfactual <- brz_did_subset
  brz_did_counterfactual$did <- 0
  
  set.seed(10)
  (res <- boot(brz_did_subset, f.did, R = 1000, newdata = brz_did_counterfactual, 
               model = didreg, parallel = "snow", ncpus = 8))
  
  yhat <- t(sapply((1:nrow(brz_did_counterfactual)), function(i){
    out <- boot.ci(res, index = i, type =c("perc"))
    with(out, c(Est = t0, pLL = percent[4], pUL = percent[5]))
  } ))
  
  brz_did_counterfactual <- cbind(brz_did_counterfactual, yhat)
  brz_did_counterfactual <- as.data.frame(subset(brz_did_counterfactual, site_ID == i))
  brz_did_agg <- as.data.frame(aggregate_by_quarter(brz_did_counterfactual, i, t_ints_real[[count]]))
  
  last_row <- tail(brz_did_agg, 1)
  
  IE <- calculate_ie(df = brz_did_agg, est_col = 'Est', case_col = 'notified_norm', pLL_col = 'pLL', pUL_col = 'pUL')
  
  did.IE <- rbind(did.IE, data.frame(IE = IE[1],
                                     lower = IE[2],
                                     upper = IE[3]))
  
  brz_did.total <- rbind(brz_did.total, brz_did_agg)
  agg.did <- rbind(agg.did, last_row)
  
  count = count + 1
}

# Aggregate IE calculation for 2x2 DiD
IE <- calculate_ie(df = agg.did, est_col = 'est.cs', case_col = 'cases.cs', pLL_col = 'pLL.cs', pUL_col = 'pUL.cs')
did.IE <- rbind(did.IE, data.frame(IE = IE[[1]], lIE = IE[[2]], uIE = IE[[3]]))

did.IE.final <- format(did.IE)

# Aggregate IE calculation by event time for 2x2 DiD
did.agg_df <- calculate.et(brz_did.total)
did.agg_df <- format(did.agg_df[2:4])

##############################################################################
##################### Regression Discontinuity Design ########################
##############################################################################
panel.IE <- data.frame()
counter <- 1 

f.rdd <- function(model, data, i, newdata) {
  tryCatch({
  require(truncreg)
  
  c.int <- model[["coefficients"]][["(Intercept)"]]
  c.cut <- model[["coefficients"]][["cutoff"]]
  c.time <- model[["coefficients"]][["rel.time"]]
  c.time2 <- model[["coefficients"]][["I((rel.time)^2)"]]
  c.ctime <- model[["coefficients"]][["cutoff:rel.time"]]
  c.ctime2 <- model[["coefficients"]][["cutoff:I((rel.time)^2)"]]
  
  m <- truncreg(notified_norm*100 ~ cutoff + rel.time + I((rel.time)^2) + cutoff:rel.time + cutoff:I((rel.time)^2),
                data = data[i, ], point = 100, direction = "right",
                start = list(count = c(c.int, c.cut, c.time, c.time2, c.ctime, c.ctime2)))
  
  mparams <- as.vector(t(do.call(rbind, as.data.frame(coef(summary(m))))[,1:2]))
  yhat <- predict(m, newdata, type = "response")
  return(c(mparams, yhat))
  }, error=function(...) NA)
}

rdd.agg <- data.frame()
rdd.IE <- data.frame()
brz_rdd <- data.frame()
count <- 1 

for (i in treatment_ids) {
  intervention.subs <- subset(interventions.df, interventions.df$site_ID == i)
  
  panel.rdd.model <- truncreg(notified_norm*100 ~ cutoff + rel.time + I((rel.time)^2) + cutoff:rel.time + cutoff:I((rel.time)^2),
                              data = intervention.subs, point = 100, direction = "right")
  
  brz_rdd_counterfactual <- intervention.subs
  brz_rdd_counterfactual$cutoff <- 0
  
  set.seed(10)
  (res <- boot(intervention.subs, f.rdd, R = 1000, newdata = brz_rdd_counterfactual, model = panel.rdd.model,
               parallel = "snow", ncpus = 8))
  
  yhat <- t(sapply((1:nrow(brz_rdd_counterfactual)), function(i){
    out <- boot.ci(res, index = i, type =c("perc"))
    with(out, c(Est = t0, pLL = percent[4], pUL = percent[5]))
  } ))
  
  # placeholder <- matrix(NA, nrow = t_ints_real[[count]]-1, ncol = ncol(yhat))
  # colnames(placeholder) <- colnames(yhat)
  # yhat_extended <- rbind(placeholder, yhat)
  
  brz_rdd_counterfactual <- cbind(brz_rdd_counterfactual, yhat)
  
  colnames(mal_rdd_counterfactual) = c("Year", "Epid Week", "intervention_site", "cases", "time_ID", "site_ID", "first.treat",
                                       "treat", "rel.time", "cutoff", "Est", "pLL", "pUL")
  
  brz_rdd_counterfactual <- as.data.frame(subset(brz_rdd_counterfactual))
  brz_rdd_agg <- as.data.frame(aggregate_by_quarter(brz_rdd_counterfactual, i, t_ints_real[[count]]))
  
  last_row <- tail(brz_rdd_agg, 1)
  
  IE <- calculate_ie(df = brz_rdd_agg, est_col = 'Est', case_col = 'notified_norm', pLL_col = 'pLL', pUL_col = 'pUL')
  
  rdd.IE <- rbind(rdd.IE, data.frame(IE = IE[1],
                                     lower = IE[2],
                                     upper = IE[3]))
  
  brz_rdd <- rbind(brz_rdd, brz_rdd_agg)
  rdd.agg <- rbind(rdd.agg, last_row)
  
  count = count + 1
}

# Aggregate IE calculation for rdd
IE <- calculate_ie(df = rdd.agg, est_col = 'est.cs', case_col = 'cases.cs', pLL_col = 'pLL.cs', pUL_col = 'pUL.cs')
rdd.IE <- rbind(rdd.IE, data.frame(IE = IE[[1]], lIE = IE[[3]], uIE = IE[[2]]))

rdd.IE.final <- format(rdd.IE)

# Aggregate IE calculation by event time for rdd
rdd.agg_df <- calculate.et(brz_rdd)
rdd.agg_df <- format(rdd.agg_df[2:4])

# Final table
i.sites.df <- c(1:32, "Agg.")

final_df <- cbind(data.frame(c(1:32, "Agg.")), pp.IE.final, rdd.IE.final, did.IE.final)
colnames(final_df) = c("Site", "Pre-post", "RDD", "DiD (2x2)")

library(writexl)
write_xlsx(rdd.IE.final, "/your/local/filepath/transfer.xlsx") # Change file location

# Final table by event time 
final_et <- cbind(pp.agg_df,did.agg_df,rdd.agg_df)
write_xlsx(final_et, "/your/local/filepath/transfer.xlsx")

################################################################################
# Final table for in-space placebo
i.sites.df <- c(33:51, "Agg.")

final_df <- cbind(data.frame(c(33:51, "Agg.")), pp.IE.final, rdd.IE.final, did.IE.final)
colnames(final_df) = c("Site", "Pre-post", "RDD", "DiD (2x2)")
write_xlsx(final_df, "/your/local/filepath/transfer.xlsx")

final_et <- cbind(pp.agg_df,rdd.agg_df, did.agg_df, placebo.et.final[1:4,])
write_xlsx(final_et, "/your/local/filepath/transfer.xlsx")
