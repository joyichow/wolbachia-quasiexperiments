lapply(c('readxl', 'tidyr', 'ggplot2', 'data.table', 'lubridate', 'boot', 'pscl', 
         'dplyr', 'doParallel', 'foreach', 'did'), require, character.only = TRUE)

# Constants
h.size <- 3.7 # Avg household size in Malaysia
i.sites <- c('MH','SF','SL','SC','AL','AF') # Intervention sites
t_ints_real <- c(250-104,227-104,227-104,255-104,221-104,245-104) # Fake treatment times 

identifiers <- c('MH','SF','SL','SC','AL','AF',
                 'M1','M2','M3','M4','M5','M6','M7','M8',
                 'S1','S2','S3','S4', 
                 "C1", "C2" ,"C3" ,"C4", 
                 "A1" ,"A2" ,"A3")

# Load and preprocess data
mal <- read_xlsx(" Dengue Incidence.xlsx") %>%
  head(332) %>%
  mutate(Time = seq(1900, length.out=332))

pop_data <- read_xlsx("site-population.xlsx") %>%
  mutate(site.pop = round(households * h.size))

pop_vector <- setNames(object = pop_data$site.pop, nm = pop_data$site)
rownames(pop_data) <- pop_data$site

mal <- mal %>%
  dplyr::mutate(across(.cols = intersect(names(mal), names(pop_vector)),
                       .fns = ~ ceiling( .x * pop_vector[cur_column()] / 100000)))
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

################################################################
##################### Pre-post analysis ########################
################################################################
# Pre-post analysis
treatment_df <- data.frame(
  intervention_site = i.sites,
  first.treat = t_ints_real, # adjusting intervention time
  treat = c(rep(1,6))
)

interventions.df <- mal_qtr[mal_qtr$intervention_site %in% i.sites,]

intervention_times <- setNames(t_ints_real, i.sites)

mal_did <- mal_qtr %>% 
  mutate(time = ifelse(time_ID >= intervention_times[intervention_site], 1, 0))

interventions.df <- left_join(interventions.df, treatment_df, by = "intervention_site")
interventions.df  <- interventions.df  %>%
  mutate(rel.time = time_ID - first.treat)

interventions.df$cutoff <- ifelse(interventions.df$rel.time >= 0, 1, 0)

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

mal_pp <- data.frame()
pp.IE <- data.frame()
agg.pp <- data.frame()

count <- 1 
x_ints_real <- c(20-8, 18-8, 18-8, 20-8, 17-8, 19-8)

for (i in i.sites) {
  intervention.subs <- subset(interventions.df, interventions.df$intervention_site == i & rel.time <= t_ints_real[[count]])
  
  offset_vector <- as.numeric(rep(log(pop_data[as.character(i), 4]), nrow(intervention.subs)))
  pp.model <- zeroinfl(cases ~ cutoff + offset(offset_vector),data = intervention.subs, dist = "negbin")
  
  mal_pp_counterfactual <- intervention.subs
  mal_pp_counterfactual$cutoff <- 0
  
  # Confidence interval calculation 
  set.seed(10) 
  (res <- boot(intervention.subs, f, R = 1000, newdata = mal_pp_counterfactual, population = offset_vector, model = pp.model,
               parallel = "snow", ncpus = 8))
  
  yhat <- t(sapply((1:nrow(mal_pp_counterfactual)), function(i){
    out <- boot.ci(res, index = i, type =c("perc"))
    with(out, c(Est = (t0), pLL = (percent[4]), pUL = (percent[5])))
  } ))
  
  mal_pp_counterfactual <- cbind(mal_pp_counterfactual, yhat)
  mal_pp_counterfactual <- as.data.frame(mal_pp_counterfactual)
  mal_pp_agg <- as.data.frame(aggregate_by_quarter(mal_pp_counterfactual, i, x_ints_real[[count]]))
  
  last_row <- tail(mal_pp_agg, 1)
  
  IE <- calculate_ie(df = mal_pp_agg, est_col = 'Est', case_col = 'cases', pLL_col = 'pLL', pUL_col = 'pUL')
  
  pp.IE <- rbind(pp.IE, data.frame(IE = IE[1],
                                   lower = IE[2],
                                   upper = IE[3]))
  
  mal_pp <- rbind(mal_pp, mal_pp_agg)
  agg.pp <- rbind(agg.pp, last_row)
  
  count = count + 1
}

# Aggregate IE calculation for pre-post
IE <- calculate_ie(df = agg.pp, est_col = 'est.cs', case_col = 'cases.cs', pLL_col = 'pLL.cs', pUL_col = 'pUL.cs')
pp.IE <- rbind(pp.IE, data.frame(IE = IE[[1]], lIE = IE[[2]], uIE = IE[[3]]))

pp.IE.final <- format(pp.IE)

# Aggregate IE calculation by event time for pre-post
pp.agg_df <- calculate.et(mal_pp)
pp.agg_df <- format(pp.agg_df[2:4])

##############################################################
############## 2x2 Difference in differences #################
##############################################################
did.IE <- data.frame()
agg.did <- data.frame()
mal_did.total <- data.frame()
count <- 1

f.did <- function(model, data, i, newdata, population) {
  require(pscl)
  offset_vector <- rep(log(population), nrow(data))
  
  c.int <- model[["coefficients"]][["count"]][["(Intercept)"]]
  c.treat <- model[["coefficients"]][["count"]][["treat"]]
  c.time <- model[["coefficients"]][["count"]][["time"]]
  c.did <- model[["coefficients"]][["count"]][["did"]]
  
  z.int <- model[["coefficients"]][["zero"]][["(Intercept)"]]
  z.treat <- model[["coefficients"]][["zero"]][["treat"]]
  z.time <- model[["coefficients"]][["zero"]][["time"]]
  z.did <- model[["coefficients"]][["zero"]][["did"]]
  
  m <- zeroinfl(cases ~ treat + time + did + offset(offset_vector),
                data = data[i, ], dist = "negbin",
                start = list(count = c(c.int, c.treat, c.time, c.did)), zero = c(list(z.int, z.treat, z.time, z.did)))
  mparams <- as.vector(t(do.call(rbind, coef(summary(m)))[, 1:2]))
  yhat <- predict(m, newdata, type = "response")
  return(c(mparams, yhat))
}

for (i in i.sites) {
  exclude_sites <- setdiff(i.sites, i)
  
  mal_did_subset <- mal_did[!mal_did$intervention_site %in% exclude_sites, ]
  mal_did_subset$time <- ifelse(mal_did_subset$time_ID >=t_ints_real[count], 1, 0)
  mal_did_subset$did <- mal_did_subset$time * mal_did_subset$treat
  mal_did_subset <- subset(mal_did_subset, time_ID <=(t_ints_real[count] + 103))
  
  offset_vector <- as.numeric(rep(log(pop_data[as.character(i), 4]), nrow(mal_did_subset)))
  didreg <- zeroinfl(cases ~ treat + time + did + offset(offset_vector), data = mal_did_subset, dist = "negbin")
  
  mal_did_counterfactual <- mal_did_subset
  mal_did_counterfactual$did <- 0
  
  set.seed(10)
  (res <- boot(mal_did_subset, f.did, R = 1000, newdata = mal_did_counterfactual, 
               population = as.numeric(pop_data[as.character(i), 4]), 
               model = didreg, parallel = "snow", ncpus = 18))
  
  yhat <- t(sapply((1:nrow(mal_did_counterfactual)), function(i){
    out <- boot.ci(res, index = i, type =c("perc"))
    with(out, c(Est = (t0), pLL = (percent[4]), pUL = (percent[5])))
  } ))
  
  mal_did_counterfactual <- cbind(mal_did_counterfactual, yhat)
  mal_did_counterfactual <- as.data.frame(mal_did_counterfactual)
  mal_did_agg <- as.data.frame(aggregate_by_quarter(mal_did_counterfactual, i, x_ints_real[[count]]))
  
  last_row <- tail(mal_did_agg, 1)
  
  IE <- calculate_ie(df = mal_did_agg, est_col = 'Est', case_col = 'cases', pLL_col = 'pLL', pUL_col = 'pUL')
  
  did.IE <- rbind(did.IE, data.frame(IE = IE[1],
                                     lower = IE[2],
                                     upper = IE[3]))
  
  mal_did.total <- rbind(mal_did.total, mal_did_agg)
  agg.did <- rbind(agg.did, last_row)
  
  count = count + 1
}

# Aggregate IE calculation for 2x2 DiD
IE <- calculate_ie(df = agg.did, est_col = 'est.cs', case_col = 'cases.cs', pLL_col = 'pLL.cs', pUL_col = 'pUL.cs')
did.IE <- rbind(did.IE, data.frame(IE = IE[[1]], lIE = IE[[2]], uIE = IE[[3]]))

did.IE.final <- format(did.IE)

# Aggregate IE calculation by event time for 2x2 DiD
did.agg_df <- calculate.et(mal_did.total)
did.agg_df <- format(did.agg_df[2:4])

##############################################################################
##################### Regression Discontinuity Design ########################
##############################################################################
interventions.df <- mal_qtr[mal_qtr$intervention_site %in% i.sites,]

interventions.df <- left_join(interventions.df, treatment_df, by = "intervention_site")
interventions.df  <- interventions.df  %>%
  mutate(rel.time = time_ID - first.treat)

interventions.df$cutoff <- ifelse(interventions.df$rel.time >= 0, 1, 0)

f.rdd <- function(model, data, i, newdata, population) {
  require(pscl)
  offset_vector <- rep(log(population), nrow(data))
  
  c.int <- model[["coefficients"]][["count"]][["(Intercept)"]]
  c.cut <- model[["coefficients"]][["count"]][["cutoff"]]
  c.time <- model[["coefficients"]][["count"]][["rel.time"]]
  c.time2 <- model[["coefficients"]][["count"]][["I((rel.time)^2)"]]
  c.ctime <- model[["coefficients"]][["count"]][["cutoff:rel.time"]]
  c.ctime2 <- model[["coefficients"]][["count"]][["cutoff:I((rel.time)^2)"]]
  
  z.int <- model[["coefficients"]][["zero"]][["(Intercept)"]]
  z.cut <- model[["coefficients"]][["zero"]][["cutoff"]]
  z.time <- model[["coefficients"]][["zero"]][["rel.time"]]
  z.time2 <- model[["coefficients"]][["zero"]][["I((rel.time)^2)"]]
  z.ctime <- model[["coefficients"]][["zero"]][["cutoff:rel.time"]]
  z.ctime2 <- model[["coefficients"]][["zero"]][["cutoff:I((rel.time)^2)"]]
  
  m <- zeroinfl(cases ~ cutoff + rel.time + I((rel.time)^2) + cutoff:rel.time + cutoff:I((rel.time)^2) + offset(offset_vector),
                data = data[i, ], dist = "negbin",
                start = list(count = c(c.int, c.cut, c.time, c.time2, c.ctime, c.ctime2), zero = c(z.int, z.cut, z.time, z.time2, z.ctime, z.ctime2)))
  mparams <- as.vector(t(do.call(rbind, coef(summary(m)))[, 1:2]))
  yhat <- predict(m, newdata, type = "response")
  return(c(mparams, yhat))
}

rdd.agg <- data.frame()
rdd.IE <- data.frame()
mal_rdd <- data.frame()
count <- 1 

for (i in i.sites) {
  intervention.subs <- subset(interventions.df, interventions.df$intervention_site == i 
                              & time_ID <=(t_ints_real[count] + 103))
  
  offset_vector <- as.numeric(rep(log(pop_data[as.character(i), 4]), nrow(intervention.subs)))
  
  panel.rdd.model <- zeroinfl(cases ~ cutoff + rel.time + I((rel.time)^2) + cutoff:rel.time + cutoff:I((rel.time)^2) + offset(offset_vector),
                              data = intervention.subs, dist = "negbin")
  
  summary(panel.rdd.model)
  
  mal_rdd_counterfactual <- intervention.subs
  mal_rdd_counterfactual$cutoff <- 0
  
  set.seed(10)
  (res <- boot(intervention.subs, f.rdd, R = 1000, newdata = mal_rdd_counterfactual, population = as.numeric(pop_data[as.character(i), 4]), model = panel.rdd.model,
               parallel = "snow", ncpus = 18))
  
  yhat <- t(sapply(((t_ints_real[[count]]):nrow(mal_rdd_counterfactual)), function(i){
    out <- boot.ci(res, index = i, type =c("perc"))
    with(out, c(Est = (t0), pLL = (percent[4]), pUL = (percent[5])))
  } ))
  
  placeholder <- matrix(0, nrow = t_ints_real[[count]]-1, ncol = ncol(yhat))
  colnames(placeholder) <- colnames(yhat)
  yhat_extended <- rbind(placeholder, yhat)
  
  mal_rdd_counterfactual <- cbind(mal_rdd_counterfactual, yhat_extended)
  
  colnames(mal_rdd_counterfactual) = c("Year", "Epid Week", "intervention_site", "cases", "time_ID", "site_ID", "first.treat",
                                       "treat", "rel.time", "cutoff", "Est", "pLL", "pUL")
  
  mal_rdd_counterfactual <- as.data.frame(subset(mal_rdd_counterfactual))
  mal_rdd_agg <- as.data.frame(aggregate_by_quarter(mal_rdd_counterfactual, i, x_ints_real[[count]]))
  
  last_row <- tail(mal_rdd_agg, 1)
  
  IE <- calculate_ie(df = mal_rdd_agg, est_col = 'Est', case_col = 'cases', pLL_col = 'pLL', pUL_col = 'pUL')
  
  rdd.IE <- rbind(rdd.IE, data.frame(IE = IE[1],
                                     lower = IE[2],
                                     upper = IE[3]))
  
  mal_rdd <- rbind(mal_rdd, mal_rdd_agg)
  rdd.agg <- rbind(rdd.agg, last_row)
  
  count = count + 1
}

# Aggregate IE calculation for rdd
IE <- calculate_ie(df = rdd.agg, est_col = 'est.cs', case_col = 'cases.cs', pLL_col = 'pLL.cs', pUL_col = 'pUL.cs')
rdd.IE <- rbind(rdd.IE, data.frame(IE = IE[[1]], lIE = IE[[2]], uIE = IE[[3]]))

rdd.IE.final <- format(rdd.IE)

# Aggregate IE calculation by event time for rdd
rdd.agg_df <- calculate.et(mal_rdd)
rdd.agg_df <- format(rdd.agg_df[2:4])

##############################################################################
############################# Final Tables ###################################
##############################################################################
final_df <- cbind(data.frame(c(1:6, "Agg.")), pp.IE.final, rdd.IE.final, did.IE.final)
final_et <- cbind(pp.agg_df,rdd.agg_df,did.agg_df, format(panel.agg[6:8]))

colnames(final_et) = colnames(final_df) = c("Site", "Pre-post", "RDD", "DiD (2x2)")
final_sheets <- list("Sheet1" = final_df, "Sheet2" = final_et)
write_xlsx(final_sheets, "/Users/joyichow/Desktop/wolb_interventions/transfer.xlsx") # Change file location
