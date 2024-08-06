##############################################################
############# Panel Difference in differences ################
##############################################################

############################## For main analysis ################################
mal_did <- left_join(mal_qtr, treatment_df, by = "intervention_site")
mal_did$first.treat[is.na(mal_did$first.treat)] <- 0
mal_did$treat[is.na(mal_did$treat)] <- 0
mal_did <- subset(mal_did, intervention_site != "Time")

set.seed(123)
dengue.attgt <- att_gt(yname = "cases",
                       gname = "first.treat",
                       idname = "site_ID",
                       tname = "time_ID",
                       data = mal_did,
                       bstrap = TRUE,
                       biters = 10,
                       pl = TRUE,
                       cores = 10)

dengue.gs <- aggte(dengue.attgt, type = "group")
summary(dengue.gs)
ggdid(dengue.gs)

count <- 1

pack.did <- data.frame()

for (i in i.sites){
  post.treat <- mean(subset(mal_did, intervention_site == i & time_ID >= t_ints_real[[count]])$cases)
  mal_did$time <- ifelse(mal_did$time_ID >=t_ints_real[count], 1, 0)

  counterfactual <- post.treat - (-0.7118) # Values obtained from dengue.gs
  up.c <- post.treat - (-1.5002)
  lo.c <- post.treat - (0.0767)

  IE <- ((counterfactual - post.treat)*100)/counterfactual
  up.IE <- ((up.c - post.treat)*100)/up.c
  lo.IE <- ((lo.c - post.treat)*100)/lo.c

  count <- count + 1
  pack.did <- rbind(pack.did, data.frame(mean_IE = IE, lower_IE = lo.IE, upper_IE = up.IE))
}

post.treat <- mean(subset(mal_did, time*treat == 1)$cases)

counterfactual <- post.treat - (-0.7118)
up.c <- post.treat - (-1.5002)
lo.c <- post.treat - (0.0767)

IE <- ((counterfactual - post.treat)*100)/counterfactual
up.IE <- ((up.c - post.treat)*100)/up.c
lo.IE <- ((lo.c - post.treat)*100)/lo.c

conf.ints <- (c(up.IE, lo.IE))

pack.did <- rbind(pack.did, data.frame(mean_IE = IE, lower_IE = min(conf.ints), upper_IE = max(conf.ints)))

pack.did.final <- format(pack.did)
################################################################################

# For in-space placebo
set.seed(123)
placebo_results <- data.frame()
avg.panel.ie <- data.frame()
placebo.et.df <- data.frame()
exclude_sites <- c('MH','SF','SL','SC','AL','AF')
mal_did_placebo <- mal_did[!mal_did$intervention_site %in% exclude_sites,]

for (i in i.sites){
  pseudo_treatment_sites <- sample(mal_did_placebo$intervention_site, 6)
  
  if (!(i %in% pseudo_treatment_sites)){
    pseudo_treatment_sites <- c(sample(pseudo_treatment_sites, 5), i)
  }
  
  t_ints_real <- c(250,227,227,255,221,245)
  t_ints_real <- sample(t_ints_real, length(pseudo_treatment_sites), replace = TRUE)
  
  treatment_df <- data.frame(
    intervention_site = pseudo_treatment_sites,
    first.treat = t_ints_real,
    treat = c(rep(1, 6)))
  
  interventions.df <- mal_qtr[mal_qtr$intervention_site %in% i.sites,]
  intervention_times <- setNames(t_ints_real, pseudo_treatment_sites)
  interventions.df <- interventions.df %>%
    mutate(time = ifelse(time_ID >= intervention_times[intervention_site], 1, 0))
  
  interventions.df <- left_join(interventions.df, treatment_df, by = "intervention_site")
  
  interventions.df$first.treat[is.na(interventions.df$first.treat)] <- 0
  interventions.df$treat[is.na(interventions.df$treat)] <- 0
  interventions.df$time[is.na(interventions.df$time)] <- 0
  
  dengue.attgt_placebo <- att_gt(yname = "cases",
                                 gname = "first.treat",
                                 idname = "site_ID",
                                 tname = "time_ID",
                                 data = interventions.df,
                                 pl = TRUE,
                                 cores = 8)
  
  # For individual sites 
  dengue.gs_placebo <- aggte(dengue.attgt_placebo, type = "group", na.rm = TRUE)
  
  post.treat <- mean(subset(interventions.df, intervention_site == i & time_ID >= treatment_df[treatment_df$intervention_site == i, "first.treat"])$cases) 
  counterfactual <-  post.treat - dengue.gs_placebo[["overall.att"]] 
  up.c <- post.treat - (counterfactual + (dengue.gs_placebo[["crit.val.egt"]]*dengue.gs_placebo[["overall.se"]]))
  lo.c <- post.treat - (counterfactual - (dengue.gs_placebo[["crit.val.egt"]]*dengue.gs_placebo[["overall.se"]]))
  
  IE <- ((counterfactual - post.treat) * 100) / counterfactual
  up.IE <- ((up.c - post.treat) * 100) / up.c
  lo.IE <- ((lo.c - post.treat) * 100) / lo.c
  
  placebo_results <- rbind(placebo_results, data.frame(Site = i, post.treat = post.treat,
                                                       counterfactual = counterfactual, up.c = up.c, lo.c = lo.c, 
                                                       mean_IE = IE, lower_IE = lo.IE, upper_IE = up.IE))
  
  # For event time
  dengue.dy <- aggte(dengue.attgt_placebo, type = "dynamic", na.rm = TRUE)
  
  dengue.dy.placebo <- data.frame(cbind(dengue.dy[["egt"]], dengue.dy[["att.egt"]], dengue.dy[["se.egt"]]))
  colnames(dengue.dy.placebo) = c("time", "att", "se")
  dengue.dy.placebo$upper.ci <- dengue.dy.placebo$att + dengue.dy[["crit.val.egt"]]*dengue.dy.placebo$se
  dengue.dy.placebo$lower.ci <- dengue.dy.placebo$att - dengue.dy[["crit.val.egt"]]*dengue.dy.placebo$se 
  
  current <- subset(interventions.df, intervention_site == i)
  current$time <- ifelse(current$time_ID >= treatment_df[treatment_df$intervention_site == i, "first.treat"], 1, 0)
  current$rel.time <- current$time_ID - current$first.treat
  
  current <- subset(current, rel.time >= 0)
  dengue.dy.placebo <- subset(dengue.dy.placebo, time >= 0)
  
  merged_df <- merge(current, dengue.dy.placebo, by.x = "rel.time", by.y = "time")
  
  merged_df <- subset(merged_df, intervention_site == i)
  
  merged_df$counterfactual <- merged_df$cases - merged_df$att
  merged_df$u.counterfactual <- merged_df$cases - merged_df$upper.ci
  merged_df$l.counterfactual <- merged_df$cases - merged_df$lower.ci
  
  placebo.et.df <- rbind(placebo.et.df, merged_df)
}

# Aggregate for site for in-space placebo
agg.IE.placebo <- (sum(placebo_results$counterfactual) - sum(placebo_results$post.treat))*100/sum(placebo_results$counterfactual)
agg.uIE.placebo <-(sum(placebo_results$up.c) - sum(placebo_results$post.treat))*100/sum(placebo_results$up.c)
agg.lIE.placebo <-(sum(placebo_results$lo.c) - sum(placebo_results$post.treat))*100/sum(placebo_results$lo.c)

placebo.df <- placebo_results[,c(1,6:8)]
placebo.df <- rbind(placebo.df, data.frame(Site = "Agg.", mean_IE = agg.IE.placebo, lower_IE = agg.uIE.placebo, upper_IE = agg.lIE.placebo))

# By event time for in-space placebo
placebo.et.df$Quarter <- ifelse(placebo.et.df$`Epid Week` <= 13, 1,
                                ifelse(placebo.et.df$`Epid Week` <= 26, 2,
                                       ifelse(placebo.et.df$`Epid Week` <= 39, 3,
                                              ifelse(placebo.et.df$`Epid Week` <= 52, 4, 4))))

placebo.et.df$YearQuarter <- paste(placebo.et.df$Year, 'Q', placebo.et.df$Quarter, sep='')

placebo.et.df <- placebo.et.df %>%
  group_by(YearQuarter, site_ID, intervention_site) %>%
  summarise(
    u.counterfactual = sum(u.counterfactual, na.rm = TRUE),
    l.counterfactual = sum(l.counterfactual, na.rm = TRUE),
    cases = sum(cases, na.rm = TRUE),
    counterfactual = sum(counterfactual, na.rm = TRUE)
  )

placebo.et.df <- placebo.et.df %>%
  arrange(intervention_site, YearQuarter) %>%
  group_by(intervention_site) %>%
  mutate(time = row_number() - 1)

placebo.et.df$time <- placebo.et.df$time * 3
breaks <- c(0, 6, 12, 18, 24, 30, 36)
labels <- c("1-6 months", "7-12 months", "13-18 months", "19-24 months", ">24 months", ">24 months")
placebo.et.df$event_time <- cut(placebo.et.df$time, breaks = breaks, labels = labels, include.lowest = TRUE)
average_cases_p <- aggregate(cbind(cases, counterfactual, u.counterfactual, l.counterfactual) ~ event_time, placebo.et.df, sum)

average_cases_p$IE <- (average_cases_p$counterfactual - average_cases_p$cases)*100/average_cases_p$counterfactual
average_cases_p$u.IE <- (average_cases_p$u.counterfactual - average_cases_p$cases)*100/average_cases_p$u.counterfactual
average_cases_p$l.IE <- (average_cases_p$l.counterfactual - average_cases_p$cases)*100/average_cases_p$l.counterfactual

placebo.et.final <- format(average_cases_p[6:8])

##############################################################
# for in time placebo
mal_did <- left_join(mal_qtr, treatment_df, by = "intervention_site")
mal_did$first.treat[is.na(mal_did$first.treat)] <- 0
mal_did$treat[is.na(mal_did$treat)] <- 0
mal_did <- subset(mal_did, intervention_site != "Time")

count = 1
mal.did.placebo <- data.frame()

for (i in i.sites){
  mal.did.new <- subset(mal_did, time_ID <= (t_ints_real[[count]] + 103) & intervention_site == i)
  
  mal.did.placebo <- rbind(mal.did.placebo, mal.did.new)
  count = count + 1
}

mal.did.placebo <- rbind(mal.did.placebo, subset(mal_did, site_ID >= 7))

set.seed(123)
dengue.attgt <- att_gt(yname = "cases",
                       gname = "first.treat",
                       idname = "site_ID",
                       tname = "time_ID",
                       data = mal.did.placebo,
                       panel = FALSE,
                       pl = TRUE,
                       cores = 8)

dengue.gs <- aggte(dengue.attgt, type = "group", na.rm = TRUE)
dengue.dy <- aggte(dengue.attgt, type = "dynamic", na.rm = TRUE)

# Extracting data from the function 
dengue.dy.df <- data.frame(cbind(dengue.dy[["egt"]], dengue.dy[["att.egt"]], dengue.dy[["se.egt"]]))
colnames(dengue.dy.df) = c("time", "att", "se")
header <- subset(mal_did, intervention_site == "AL" & time_ID <= 104)
dengue.dy.df <- cbind(subset(dengue.dy.df, time >= 0), header$Year, header$`Epid Week`)
colnames(dengue.dy.df) = c("time", "att", "se", "Year", "Epid Week")

dengue.dy.df$upper.ci <- dengue.dy.df$att + dengue.dy[["crit.val.egt"]]*dengue.dy.df$se
dengue.dy.df$lower.ci <- dengue.dy.df$att - dengue.dy[["crit.val.egt"]]*dengue.dy.df$se 

dengue.dy.df$Quarter <- ifelse(dengue.dy.df$`Epid Week` <= 13, 1,
                               ifelse(dengue.dy.df$`Epid Week` <= 26, 2,
                                      ifelse(dengue.dy.df$`Epid Week` <= 39, 3,
                                             ifelse(dengue.dy.df$`Epid Week` <= 52, 4, 4))))

dengue.dy.df$YearQuarter <- paste(dengue.dy.df$Year, 'Q', dengue.dy.df$Quarter, sep='')

# Getting case data 
panel.did.et <- data.frame() 
count <- 1

for (i in i.sites){
  current <- subset(mal_did, intervention_site == i)
  current$time <- ifelse(current$time_ID >=t_ints_real[count], 1, 0)
  current$rel.time <- current$time_ID - current$first.treat
  
  current$Quarter <- ifelse(current$`Epid Week` <= 13, 1,
                            ifelse(current$`Epid Week` <= 26, 2,
                                   ifelse(current$`Epid Week` <= 39, 3,
                                          ifelse(current$`Epid Week` <= 52, 4, 4))))
  
  current$YearQuarter <- paste(current$Year, 'Q', current$Quarter, sep='')
  current$site <- i
  current$site_ID <- count
  
  panel.did.et <- rbind(panel.did.et, current)
  count = count + 1
}

panel.did.et <- subset(panel.did.et, rel.time >= 0)
panel_et <- merge(panel.did.et, dengue.dy.df, by.x = "rel.time", by.y = "time")
panel_et$counterfactual <- panel_et$cases - panel_et$att
panel_et$u.counterfactual <- panel_et$cases - panel_et$upper.ci
panel_et$l.counterfactual <- panel_et$cases - panel_et$lower.ci

panel.agg <- panel_et %>%
  group_by(YearQuarter.x, site_ID, intervention_site) %>%
  summarise(
    u.counterfactual = sum(u.counterfactual, na.rm = TRUE),
    l.counterfactual = sum(l.counterfactual, na.rm = TRUE),
    cases = sum(cases, na.rm = TRUE),
    counterfactual = sum(counterfactual, na.rm = TRUE)
  )

panel.agg <- panel.agg %>%
  arrange(intervention_site, YearQuarter.x) %>%
  group_by(intervention_site) %>%
  mutate(time = row_number() - 1)

panel.agg$time <- panel.agg$time*3

breaks <- c(0, 6, 12, 18, 24, 30)
labels <- c("1-6 months", "7-12 months", "13-18 months", "19-24 months", ">24 months")

panel.agg$event_time <- cut(panel.agg$time, breaks = breaks, labels = labels, include.lowest = TRUE)

panel.agg <- aggregate(cbind(cases, counterfactual, u.counterfactual, l.counterfactual) ~ event_time, panel.agg, sum)
panel.agg$IE <- (panel.agg$counterfactual - panel.agg$cases)*100/panel.agg$counterfactual
panel.agg$u.IE <- (panel.agg$u.counterfactual - panel.agg$cases)*100/panel.agg$u.counterfactual
panel.agg$l.IE <- (panel.agg$l.counterfactual - panel.agg$cases)*100/panel.agg$l.counterfactual

# Individual site IE 
count <- 1
pack.did <- data.frame()

for (i in i.sites){
  post.treat <- mean(subset(mal_did, intervention_site == i & time_ID >= t_ints_real[[count]])$cases) 
  mal_did$time <- ifelse(mal_did$time_ID >=t_ints_real[count], 1, 0)
  
  counterfactual <- post.treat - (1.1598) # Values obtained from dengue.gs
  up.c <- post.treat - (1.5866)
  lo.c <- post.treat - (0.733)
  
  IE <- ((counterfactual - post.treat)*100)/counterfactual
  up.IE <- ((up.c - post.treat)*100)/up.c
  lo.IE <- ((lo.c - post.treat)*100)/lo.c
  
  count <- count + 1
  pack.did <- rbind(pack.did, data.frame(mean_IE = IE, lower_IE = lo.IE, upper_IE = up.IE)) 
}

post.treat <- mean(subset(mal_did, time*treat == 1)$cases) 

counterfactual <- post.treat - (-0.7118)
up.c <- post.treat - (-1.5002)
lo.c <- post.treat - (0.0767)

IE <- ((counterfactual - post.treat)*100)/counterfactual
up.IE <- ((up.c - post.treat)*100)/up.c
lo.IE <- ((lo.c - post.treat)*100)/lo.c

conf.ints <- (c(up.IE, lo.IE))

pack.did <- rbind(pack.did, data.frame(mean_IE = IE, lower_IE = min(conf.ints), upper_IE = max(conf.ints))) 

pack.did.final <- format(pack.did)