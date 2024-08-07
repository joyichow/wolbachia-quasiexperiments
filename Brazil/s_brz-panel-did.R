# Panel DiD
intervention_times <- setNames(t_ints_real, treatment_ids)

brz_did <- brazil.den %>%
  mutate(time = ifelse(time_ID >= intervention_times[treatment_ids], 1, 0))

brz_did <- left_join(brazil.den, treatment_df, by = "site_ID")
brz_did$first.treat[is.na(brz_did$first.treat)] <- 0
brz_did$treat[is.na(brz_did$treat)] <- 0

set.seed(123)
dengue.attgt <- att_gt(yname = "notified_norm",
                       gname = "first.treat",
                       idname = "site_ID",
                       tname = "time_ID",
                       data = brz_did,
                       pl = TRUE,
                       cores = 8)

pre.test <- conditional_did_pretest(yname = "notified_norm",
                                    tname = "time_ID",
                                    idname = "site_ID",
                                    gname = "first.treat",
                                    xformla=~1,
                                    data = brz_did)

summary(pre.test)

dengue.gs <- aggte(dengue.attgt, type = "group")
dengue.dy <- aggte(dengue.attgt, type = "dynamic")

# Extracting data from the function 
dengue.dy.df <- data.frame(cbind(dengue.dy[["egt"]], dengue.dy[["att.egt"]], dengue.dy[["se.egt"]]))
colnames(dengue.dy.df) = c("time", "att", "se")
dengue.dy.df$upper.ci <- dengue.dy.df$att + dengue.dy[["crit.val.egt"]]*dengue.dy.df$se
dengue.dy.df$lower.ci <- dengue.dy.df$att - dengue.dy[["crit.val.egt"]]*dengue.dy.df$se 

panel.did.et <- data.frame() 
count <- 1

# IE by time 
for (i in treatment_ids){
  current <- subset(brz_did, site_ID == i)
  current$time <- ifelse(current$time_ID >=t_ints_real[count], 1, 0)
  current$rel.time <- current$time_ID - current$first.treat
  
  panel.did.et <- rbind(panel.did.et, current)
  count = count + 1
}

panel.did.et <- subset(panel.did.et, rel.time >= 0)
dengue.dy.df <- subset(dengue.dy.df, time >= 0)

merged_df <- merge(panel.did.et, dengue.dy.df, by.x = "rel.time", by.y = "time")
merged_df$time <- merged_df$rel.time * 3

merged_df$counterfactual <- merged_df$notified_norm - merged_df$att
merged_df$u.counterfactual <- merged_df$notified_norm - merged_df$upper.ci
merged_df$l.counterfactual <- merged_df$notified_norm - merged_df$lower.ci

breaks <- c(0, 6, 12, 18, 24, 30)
labels <- c("1-6 months", "7-12 months", "13-18 months", "19-24 months", ">24 months")
merged_df$event_time <- cut(merged_df$time, breaks = breaks, labels = labels, include.lowest = TRUE)
c <- aggregate(cbind(notified_norm, counterfactual, u.counterfactual, l.counterfactual) ~ event_time, merged_df, sum)

average_cases$IE <- (average_cases$counterfactual - average_cases$notified_norm)*100/average_cases$counterfactual
average_cases$u.IE <- (average_cases$u.counterfactual - average_cases$notified_norm)*100/average_cases$u.counterfactual
average_cases$l.IE <- (average_cases$l.counterfactual - average_cases$notified_norm)*100/average_cases$l.counterfactual

count <- 1

pack.did <- data.frame()

# Individual site calculation
for (i in treatment_ids){
  post.treat <- mean(subset(brz_did, site_ID == i & time_ID >= t_ints_real[[count]])$notified_norm)
  brz_did$time <- ifelse(brz_did$time_ID >=t_ints_real[count], 1, 0)
  
  counterfactual <- post.treat - (-0.013)
  up.c <- post.treat - (-0.0414)
  lo.c <- post.treat - (0.0154)
  
  IE <- ((counterfactual - post.treat)*100)/counterfactual
  up.IE <- ((up.c - post.treat)*100)/up.c
  lo.IE <- ((lo.c - post.treat)*100)/lo.c
  
  conf.ints <- (c(up.IE, lo.IE))
  
  count <- count + 1
  pack.did <- rbind(pack.did, data.frame(mean_IE = IE, lower_IE = min(conf.ints), upper_IE = max(conf.ints))) 
}

# Aggregate calculation
post.treat <- mean(subset(brz_did, time*treat == 1)$notified_norm) 

counterfactual <- post.treat - (-0.7118)
up.c <- post.treat - (-1.5002)
lo.c <- post.treat - (0.0767)

IE <- ((counterfactual - post.treat)*100)/counterfactual
up.IE <- ((up.c - post.treat)*100)/up.c
lo.IE <- ((lo.c - post.treat)*100)/lo.c

pack.did.final <- format(pack.did)
write_xlsx(pack.did.final, "/your/local/filepath/transfer.xlsx")

# In-space placebo calculation
set.seed(123)
placebo_results <- data.frame()
avg.panel.ie <- data.frame()
placebo.et.df <- data.frame()
exclude_sites <- c(1:32)
brz_did_placebo <- brz_did[!brz_did$site_ID %in% exclude_sites,]

for (i in treatment_ids){
  pseudo_treatment_sites <- sample(brz_did_placebo$site_ID, 4)
  
  if (!(i %in% pseudo_treatment_sites)){
    pseudo_treatment_sites <- c(sample(pseudo_treatment_sites, 3), i)
  }
  
  t_ints_real <- c(42, 43, 44, 51)
  t_ints_real <- sample(t_ints_real, length(pseudo_treatment_sites), replace = TRUE)
  
  treatment_df <- data.frame(
    site_ID = pseudo_treatment_sites,
    first.treat = t_ints_real,
    treat = c(rep(1, 4)))
  
  interventions.df <- brazil.den[brazil.den$site_ID %in% treatment_ids,]
  intervention_times <- setNames(t_ints_real, pseudo_treatment_sites)
  interventions.df <- interventions.df %>%
    mutate(time = ifelse(time_ID >= intervention_times[site_ID], 1, 0))
  
  interventions.df <- left_join(interventions.df, treatment_df, by = "site_ID")
  
  interventions.df$first.treat[is.na(interventions.df$first.treat)] <- 0
  interventions.df$treat[is.na(interventions.df$treat)] <- 0
  interventions.df$time[is.na(interventions.df$time)] <- 0
  
  dengue.attgt_placebo <- att_gt(yname = "notified_norm",
                                 gname = "first.treat",
                                 idname = "site_ID",
                                 tname = "time_ID",
                                 data = interventions.df,
                                 pl = TRUE,
                                 cores = 8)
  
  print(i)
  
  pre.test <- conditional_did_pretest(yname = "notified_norm",
                                      tname = "time_ID",
                                      idname = "site_ID",
                                      gname = "first.treat",
                                      xformla=~1,
                                      data = brz_did)
  
  summary(pre.test)
  
  # For individual sites 
  dengue.gs_placebo <- aggte(dengue.attgt_placebo, type = "group", na.rm = TRUE)
  
  post.treat <- mean(subset(interventions.df, site_ID == i & time_ID >= treatment_df[treatment_df$site_ID == i, "first.treat"])$notified_norm) 
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
  
  current <- subset(interventions.df, site_ID == i)
  current$time <- ifelse(current$time_ID >= treatment_df[treatment_df$site_ID == i, "first.treat"], 1, 0)
  current$rel.time <- current$time_ID - current$first.treat
  
  current <- subset(current, rel.time >= 0)
  dengue.dy.placebo <- subset(dengue.dy.placebo, time >= 0)
  
  merged_df <- merge(current, dengue.dy.placebo, by.x = "rel.time", by.y = "time")
  
  merged_df <- subset(merged_df, site_ID == i)
  
  merged_df$counterfactual <- merged_df$notified_norm - merged_df$att
  merged_df$u.counterfactual <- merged_df$notified_norm - merged_df$upper.ci
  merged_df$l.counterfactual <- merged_df$notified_norm - merged_df$lower.ci
  
  placebo.et.df <- rbind(placebo.et.df, merged_df)
}

# Aggregate for site for in-space placebo
agg.IE.placebo <- (sum(placebo_results$counterfactual) - sum(placebo_results$post.treat))*100/sum(placebo_results$counterfactual)
agg.uIE.placebo <-(sum(placebo_results$up.c) - sum(placebo_results$post.treat))*100/sum(placebo_results$up.c)
agg.lIE.placebo <-(sum(placebo_results$lo.c) - sum(placebo_results$post.treat))*100/sum(placebo_results$lo.c)

placebo.df <- placebo_results[,c(1,6:8)]
placebo.df <- rbind(placebo.df, data.frame(Site = "Agg.", mean_IE = agg.IE.placebo, lower_IE = agg.uIE.placebo, upper_IE = agg.lIE.placebo))

# By event time for in-space placebo
placebo.et.df$time <- placebo.et.df$rel.time * 3
breaks <- c(0, 6, 12, 18, 24, 30, 36)
labels <- c("1-6 months", "7-12 months", "13-18 months", "19-24 months", ">24 months", ">24 months")
placebo.et.df$event_time <- cut(placebo.et.df$time, breaks = breaks, labels = labels, include.lowest = TRUE)
average_cases_p <- aggregate(cbind(notified_norm, counterfactual, u.counterfactual, l.counterfactual) ~ event_time, placebo.et.df, sum)

average_cases_p$IE <- (average_cases_p$counterfactual - average_cases_p$notified_norm)*100/average_cases_p$counterfactual
average_cases_p$u.IE <- (average_cases_p$u.counterfactual - average_cases_p$notified_norm)*100/average_cases_p$u.counterfactual
average_cases_p$l.IE <- (average_cases_p$l.counterfactual - average_cases_p$notified_norm)*100/average_cases_p$l.counterfactual

placebo.et.final <- format(average_cases_p[6:8])

# Panel DiD
intervention_times <- setNames(t_ints_real, treatment_ids)
brz_did <- brazil.den %>%
  mutate(time = ifelse(time_ID >= intervention_times[treatment_ids], 1, 0))

brz_did <- left_join(brazil.den, treatment_df, by = "site_ID")
brz_did$first.treat[is.na(brz_did$first.treat)] <- 0
brz_did$treat[is.na(brz_did$treat)] <- 0

count = 1
brz.did.placebo <- data.frame()

for (i in treatment_ids){
  brz.did.new <- subset(brz_did, time_ID <= (t_ints_real[[count]] + 7) & site_ID == i)
  
  brz.did.placebo <- rbind(brz.did.placebo, brz.did.new)
  count = count + 1
}

brz.did.placebo <- rbind(brz.did.placebo, subset(brz_did, site_ID >= 33))

set.seed(123)
dengue.attgt <- att_gt(yname = "notified_norm",
                       gname = "first.treat",
                       idname = "site_ID",
                       tname = "time_ID",
                       data = brz.did.placebo,
                       panel = FALSE,
                       pl = TRUE,
                       cores = 8)

dengue.gs <- aggte(dengue.attgt, type = "group", na.rm = TRUE)
dengue.dy <- aggte(dengue.attgt, type = "dynamic", na.rm = TRUE)

# Extracting data from the function 
dengue.dy.df <- data.frame(cbind(dengue.dy[["egt"]], dengue.dy[["att.egt"]], dengue.dy[["se.egt"]]))
colnames(dengue.dy.df) = c("time", "att", "se")
dengue.dy.df$upper.ci <- dengue.dy.df$att + dengue.dy[["crit.val.egt"]]*dengue.dy.df$se
dengue.dy.df$lower.ci <- dengue.dy.df$att - dengue.dy[["crit.val.egt"]]*dengue.dy.df$se 

panel.did.et <- data.frame() 
count <- 1

# IE by time 
for (i in treatment_ids){
  current <- subset(brz_did, site_ID == i)
  current$time <- ifelse(current$time_ID >=t_ints_real[count], 1, 0)
  current$rel.time <- current$time_ID - current$first.treat
  current$site <- i
  
  panel.did.et <- rbind(panel.did.et, current)
  count = count + 1
}

panel.did.et <- subset(panel.did.et, rel.time >= 0 & rel.time <= 7)
dengue.dy.df <- subset(dengue.dy.df, time >= 0 & time <= 7)

panel_et <- merge(panel.did.et, dengue.dy.df, by.x = "rel.time", by.y = "time")
panel_et$counterfactual <- panel_et$notified_norm - panel_et$att
panel_et$u.counterfactual <- panel_et$notified_norm - panel_et$upper.ci
panel_et$l.counterfactual <- panel_et$notified_norm - panel_et$lower.ci

panel_et$time <- panel_et$rel.time * 3
breaks <- c(0, 6, 12, 18, 24, 30)
labels <- c("1-6 months", "7-12 months", "13-18 months", "19-24 months", ">24 months")
panel_et$event_time <- cut(panel_et$time, breaks = breaks, labels = labels, include.lowest = TRUE)
average_cases <- aggregate(cbind(notified_norm, counterfactual, u.counterfactual, l.counterfactual) ~ event_time, panel_et, sum)

average_cases$IE <- (average_cases$counterfactual - average_cases$notified_norm)*100/average_cases$counterfactual
average_cases$u.IE <- (average_cases$u.counterfactual - average_cases$notified_norm)*100/average_cases$u.counterfactual
average_cases$l.IE <- (average_cases$l.counterfactual - average_cases$notified_norm)*100/average_cases$l.counterfactual

count <- 1

pack.did <- data.frame()

for (i in treatment_ids){
  post.treat <- mean(subset(brz_did, site_ID == i & time_ID >= t_ints_real[[count]])$notified_norm)
  brz_did$time <- ifelse(brz_did$time_ID >=t_ints_real[count], 1, 0)
  
  counterfactual <- post.treat - (-0.0231)
  up.c <- post.treat - (0.0067)
  lo.c <- post.treat - (-0.0528)
  
  IE <- ((counterfactual - post.treat)*100)/counterfactual
  up.IE <- ((up.c - post.treat)*100)/up.c
  lo.IE <- ((lo.c - post.treat)*100)/lo.c
  
  conf.ints <- (c(up.IE, lo.IE))
  
  count <- count + 1
  pack.did <- rbind(pack.did, data.frame(mean_IE = IE, lower_IE = min(conf.ints), upper_IE = max(conf.ints))) 
}

post.treat <- mean(subset(brz_did, time*treat == 1)$notified_norm) 

counterfactual <- post.treat - (-0.7118)
up.c <- post.treat - (-1.5002)
lo.c <- post.treat - (0.0767)

IE <- ((counterfactual - post.treat)*100)/counterfactual
up.IE <- ((up.c - post.treat)*100)/up.c
lo.IE <- ((lo.c - post.treat)*100)/lo.c

pack.did.final <- format(pack.did)
write_xlsx(pack.did.final, "/your/local/filepath/transfer.xlsx")

count <- 1

for (i in treatment_ids){
  brz_did$time <- ifelse(brz_did$time_ID >= t_ints_real[[count]] & brz_did$site_ID <= 32, 1, 0)
}

post.treat <- mean(subset(brz_did, time == 1)$notified_norm) 

counterfactual <- post.treat - (-0.0231)
up.c <- post.treat - (0.0067)
lo.c <- post.treat - (-0.0528)

IE <- ((counterfactual - post.treat)*100)/counterfactual
up.IE <- ((up.c - post.treat)*100)/up.c
lo.IE <- ((lo.c - post.treat)*100)/lo.c

conf.ints <- (c(up.IE, lo.IE))

pack.did <- rbind(pack.did, data.frame(mean_IE = IE, lower_IE = min(conf.ints), upper_IE = max(conf.ints))) 

pack.did.final <- format(pack.did)
write_xlsx(pack.did.final, "/your/local/filepath/transfer.xlsx")