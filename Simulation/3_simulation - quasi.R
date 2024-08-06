libraries <- c(
  "readxl", "tidyr", "lubridate", "zoo", "data.table", 
  "doParallel", "foreach", "dplyr", "tidyverse", "did"
)

for (lib in libraries) {
  if (!require(lib, character.only = TRUE)) {
    install.packages(lib, dependencies = TRUE)
    library(lib, character.only = TRUE)
  }
}

# Set up the cluster using all available cores
numCores <- 6
cl <- makeCluster(numCores)
registerDoParallel(cl)

control_sizes <- c(15, 20, 25, 30, 35, 40) # Control sizes
iterations <- as.numeric(1000)

results_dir <- "simulation_results"
dir.create(results_dir, showWarnings = FALSE)

deviation_results <- list()
intermediate_filepath <- "simulation_results_size_%d_iter_%d.Rdata"

for (size in control_sizes) {
  foreach(iter = 1:iterations, .combine = c, .packages = c("readxl", "tidyr", "lubridate", "zoo", "data.table", 
                                                           "dplyr", "tidyverse", "did", "foreach", "doParallel")) %dopar% {
    
    
    source("/Users/joyichow/Desktop/wolb_interventions/dengue simulator/utility.R")
    source("/Users/joyichow/Desktop/CSCM/CSCM_helper_functions.R")
    source("/Users/joyichow/Desktop/wolb_interventions/dengue simulator/dengue simulation.R")
    source("/Users/joyichow/Desktop/wolb_interventions/simulation_helperfxns.R")
                                                             
    tryCatch({
      
      set.seed(123 + iter)
      size <- as.numeric(size)
      print(paste("Size:", size, "Iter:", iter))
      
      # 1. Create intervention sites and also store intervention site data without 
      sim.sites <- data.frame()      # for all intervention sites without intervention
      sim.stats <- data.frame()      # to store timing 
      sim.isites.non <- data.frame() # for intervention sites without intervention
      
      deviation_df <- data.frame() # Placeholder for the simulation results
      
      for (i in 1:2){
        # generate random numbers for human population, start, end, wolb population 
        human_pop <- floor(runif(1, min=1000, max=15000))
        ihuman_pop <- floor(runif(1, min=10, max=500))
        wt_mos_pop <- floor(runif(1, min=1000, max=2000))
        start <- floor(runif(1, min=300, max=500))
        
        # set up initials for population
        {
          Sh[1]<- human_pop
          Eh[1]<- 0
          Ih[1]<- ihuman_pop
          Rh[1]<- 0
          
          An[1]<- wt_mos_pop
          Sn[1]<- 0
          En[1]<- 0
          In[1]<- 20
          
          Aw[1]<- 0
          Sw[1]<- 0
          Ew[1]<- 0
          Iw[1]<- 0
        }
        
        # intervention site with no intervention
        intervention <- generate_intervention(c_intervention, start, 1000, 0, μw)
        result <- run_iit(iteration, mosquitoes)
        result$cumulative_intervention <- intervention$cumulative
        
        # get poissoned report
        poissioned_Ih<-add_poisson_noise(result$Ih)
        result$Ih <- poissioned_Ih
        
        sim.isites.non <- rbind(sim.isites.non, poissioned_Ih)
        
        # simulate male-wolbachia mosquitoes release
        min_wolb_mos <- as.numeric(sum(result$An[start:1001]))
        
        mos_pop <- floor(runif(1, min=min_wolb_mos*1.2, max=min_wolb_mos*2))
        
        sim.stats <- rbind(sim.stats, data.frame(start = start,
                                                 hpop = human_pop, 
                                                 Ih_intial = ihuman_pop, 
                                                 mpop = mos_pop, 
                                                 wt_mpop = wt_mos_pop,
                                                 total_wt_mos = min_wolb_mos))
        
        intervention<- generate_intervention(c_intervention, start, 1000, mos_pop, μw)
        
        # run compartments model with population initials and released mosquito track
        result <- run_iit(iteration, mosquitoes)
        result$cumulative_intervention <- intervention$cumulative
        
        # get poissoned report
        poissioned_Ih<-add_poisson_noise(result$Ih)
        result$Ih <- poissioned_Ih
        
        sim.sites <- rbind(sim.sites, poissioned_Ih)
      }
      
      sim.isites.non <- t(sim.isites.non)
      sim.sites <- t(sim.sites)
      
      # 2. Create control sites 
      set.seed(143)
      for (i in 1:size){
        # generate random numbers for human population, start, end, wolb population 
        human_pop <- floor(runif(1, min=1000, max=15000))
        start <- floor(runif(1, min=300, max=500))
        ihuman_pop <- floor(runif(1, min=10, max=500))
        wt_mos_pop <- floor(runif(1, min=1000, max=10000))
        
        # set up initials for population
        {
          Sh[1]<- human_pop
          Eh[1]<- 0
          Ih[1]<- ihuman_pop
          Rh[1]<- 0
          
          An[1]<- wt_mos_pop 
          Sn[1]<- 0
          En[1]<- 0
          In[1]<- 20
          
          Aw[1]<- 0
          Sw[1]<- 0
          Ew[1]<- 0
          Iw[1]<- 0
        }
        
        # simulate male-wolbachia mosquitoes release
        intervention<- generate_intervention(c_intervention, start, 1000, 0, μw)
        
        # run compartments model with population initials and released mosquito track
        result <- run_iit(iteration, mosquitoes)
        
        min_wolb_mos <- as.numeric(sum(result$An[start:1001]))
        
        sim.stats <- rbind(sim.stats, data.frame(start = start,
                                                 hpop = human_pop, 
                                                 Ih_intial = ihuman_pop, 
                                                 mpop = 0, 
                                                 wt_mpop = wt_mos_pop,
                                                 total_wt_mos = min_wolb_mos))
        
        # get poissoned report
        poissioned_Ih<-add_poisson_noise(result$Ih)
        sim.sites <- cbind(sim.sites, poissioned_Ih)
      }
      
      rm(result)
      
      ############ Bit of data cleaning#############
      sim.sites <- sim.sites[-(1:150), , drop = FALSE]
      sim.isites.non <- sim.isites.non[-(1:150), , drop = FALSE]
      rownames(sim.isites.non) = rownames(sim.sites) = as.character(1:851)
      colnames(sim.sites) <- as.character(1:(size+2))
      
      # changing to rate incidence: (cases/pop)*100000
      population <- sim.stats$hpop
      rate.sim.sites <- sweep(sim.sites, 2, population, FUN = "/") * 100000
      rate.sim.sites.non <- sweep(sim.isites.non, 2, c(population[1:2]), FUN = "/") * 100000
      
      # converting to long format 
      long.sim <- as.data.frame(rate.sim.sites) %>% 
        mutate(time_ID = row_number()) %>% 
        pivot_longer(
          cols = -time_ID,  
          names_to = "intervention_site", 
          values_to = "cases"
        )
      
      long.sim$site_ID <- as.numeric(long.sim$intervention_site)
      long.sim <- as.data.frame(long.sim)
      
      treatment_ids <- c(1:2)
      t_ints_real <- c(sim.stats$start[1:2]-150)
      
      ###############################################################################
      ########################## Calculating the true IE  ###########################
      ###############################################################################
      true_ie <- data.frame()
      agg_true <- data.frame()
      
      for (i in treatment_ids){
        new_df <- cbind(rate.sim.sites[,i], rate.sim.sites.non[,i], day = c(1:851))
        new_df <- as.data.frame(new_df)
        new_df$day <- as.numeric(new_df$day)
        
        new_df <- subset(as.data.frame(new_df), day >= as.numeric(t_ints_real[i]))
        
        new_df$cases_cumsum <- cumsum(new_df$V2)
        new_df$intervention_cumsum <- cumsum(new_df$V1)
        
        new_df$IE <- (new_df$cases_cumsum - new_df$intervention_cumsum)*100 / new_df$cases_cumsum
        
        ie <- (tail(new_df, 1))$IE
        
        true_ie <- rbind(true_ie, data.frame(site = i, IE = ie))
        agg_true <- rbind(agg_true, tail(new_df, 1))
      }
      
      agg.true <- (sum(agg_true$cases_cumsum) - sum(agg_true$intervention_cumsum))*100/sum(agg_true$cases_cumsum)
      
      true_ie <- rbind(true_ie, agg.true)
      true_ie <- format(true_ie[2])
      
      rm(agg_true)
      
      ###############################################################################
      ####################### Synthetic control - SCM & cSCM ########################
      ###############################################################################
      long.sim$intervention_site <- as.numeric(long.sim$intervention_site)
      simulated_scm_results <- list()
      
      for(i in 1:2){
        simulated_scm_results[[i]] <- simulation.scm(i, long.sim, treatment_ids, t_ints_real)
      }
      
      # Individual site IEs
      sim.sc.ie <- data.frame()
      
      ie.extract <- lapply(simulated_scm_results, function(df) {
        df <- df[, 5:9]
        df <- df[nrow(df), ]
        return(df)
      })
      
      sim.sc.ie <- do.call(rbind, ie.extract)
      sim.sc.ie <- sim.sc.ie[, c(1:3, 5, 4)]
      
      agg.SCM <- (sum(sim.sc.ie $SCM_cumsum) - sum(sim.sc.ie $cases_cumsum))*100/sum(sim.sc.ie $SCM_cumsum)
      agg.CSCM <- (sum(sim.sc.ie$CSCM_cumsum) - sum(sim.sc.ie $cases_cumsum))*100/sum(sim.sc.ie $CSCM_cumsum)
      
      sim.sc.ie <- rbind(sim.sc.ie[,4:5], data.frame(IE_SCM = agg.SCM, IE_CSCM = agg.CSCM))
      
      scm.ie <- format(sim.sc.ie[1])
      cscm.ie <- format(sim.sc.ie[2])
      
      rm(simulated_scm_results)
      
      ###############################################################################
      ################### Data cleaning for other quasi-methods #####################
      ###############################################################################
      sim.cases <- as.data.frame(sim.sites) %>% 
        mutate(time_ID = row_number()) %>% 
        pivot_longer(
          cols = -time_ID,  
          names_to = "intervention_site", 
          values_to = "cases"
        )
      
      sim.cases <- as.data.frame(sim.cases)
      
      pop_vector <- setNames(object = sim.stats$hpop, nm = c(1:(size+2)))
      
      ###############################################################################
      ############################### Pre-post ######################################
      ###############################################################################
      i.sites <- c(1:2)
      
      x_ints_real <- t_ints_real
      
      treatment_df <- data.frame(
        intervention_site = i.sites,
        first.treat = t_ints_real,
        treat = c(rep(1,2))
      )
      
      interventions.df <- sim.cases[sim.cases$intervention_site %in% i.sites,]
      interventions.df$intervention_site <- as.integer(interventions.df$intervention_site)
      intervention_times <- setNames(t_ints_real, i.sites)
      
      interventions.df <- left_join(interventions.df, treatment_df, by = "intervention_site")
      
      interventions.df  <- interventions.df  %>%
        mutate(rel.time = time_ID - first.treat)
      
      interventions.df$cutoff <- ifelse(interventions.df$rel.time >= 0, 1, 0)
      
      sim_pp <- data.frame()
      pp.IE <- data.frame()
      agg.pp <- data.frame()
      
      for (i in i.sites) {
        intervention.subs <- subset(interventions.df, interventions.df$intervention_site == i)
        
        offset_vector <- as.numeric(rep(log(pop_vector[i]), nrow(intervention.subs)))
        pp.model <- glm(cases ~ cutoff + offset(offset_vector), family = quasipoisson(link = "log"), data = intervention.subs)
        
        pp_counterfactual <- intervention.subs
        pp_counterfactual$cutoff <- 0
        
        yhat <- predict(pp.model, newdata = pp_counterfactual, type = "response")
        
        pp_counterfactual <- cbind(pp_counterfactual, yhat)
        pp_counterfactual <- as.data.frame(pp_counterfactual)
        pp_agg <- as.data.frame(aggregate_by_quarter(pp_counterfactual, i, x_ints_real[[i]]))
        
        last_row <- tail(pp_agg, 1)
        
        IE <- calculate_ie(df = pp_agg, est_col = 'yhat', case_col = 'cases')
        
        pp.IE <- rbind(pp.IE, data.frame(IE = IE[1]))
        
        sim_pp <- rbind(sim_pp, pp_agg)
        agg.pp <- rbind(agg.pp, last_row)
      }
      
      # Aggregate IE calculation for pre-post
      IE <- calculate_ie(df = agg.pp, est_col = 'est.cs', case_col = "cases.cs")
      pp.IE <- rbind(pp.IE, data.frame(IE = IE[[1]]))
      
      pp.IE.final <- format(pp.IE)
      
      ###############################################################################
      ############################### Panel DiD #####################################
      ###############################################################################
      sim.cases$intervention_site <- as.numeric(sim.cases$intervention_site)
      
      sim_did <- left_join(sim.cases, treatment_df, by = "intervention_site")
      sim_did$first.treat[is.na(sim_did$first.treat)] <- 0
      sim_did$treat[is.na(sim_did$treat)] <- 0
      
      set.seed(123)
      
      dengue.attgt <- att_gt(yname = "cases",
                             gname = "first.treat",
                             idname = "intervention_site",
                             tname = "time_ID",
                             data = sim_did,
                             bstrap = TRUE,
                             biters = 1000,
                             pl = TRUE,
                             cores = 10)
      
      dengue.gs <- aggte(dengue.attgt, type = "group")
      
      pack.did <- data.frame()
      
      for (i in i.sites){
        post.treat <- mean(subset(sim_did, intervention_site == i & time_ID >= t_ints_real[[i]])$cases)
        sim_did$time <- ifelse(sim_did$time_ID >=t_ints_real[[i]], 1, 0)
        
        counterfactual <- post.treat - (dengue.gs[["overall.att"]]) # Values obtained from dengue.gs
        
        IE <- ((counterfactual - post.treat)*100)/counterfactual
        
        pack.did <- rbind(pack.did, data.frame(mean_IE = IE))
      }
      
      # Aggregate calculation 
      post.treat <- mean(subset(sim_did, treat*first.treat >= 1)$cases) 
      counterfactual <- post.treat - (dengue.gs[["overall.att"]]) # Values obtained from dengue.gs
      
      IE <- ((counterfactual - post.treat)*100)/counterfactual
      
      pack.did <- rbind(pack.did, data.frame(mean_IE = IE)) 
      pack.did.final <- format(pack.did)
      
      ###############################################################################
      ########################## Difference-in-differences ##########################
      ###############################################################################
      
      did.IE <- data.frame()
      agg.did <- data.frame()
      did.total <- data.frame()
      
      for (i in i.sites) {
        exclude_sites <- setdiff(i.sites, i)
        
        did_subset <- sim_did[!sim_did$intervention_site %in% exclude_sites, ]
        did_subset$time <- ifelse(did_subset$time_ID >=t_ints_real[i], 1, 0)
        did_subset$did <- did_subset$time * did_subset$treat
        
        offset_vector <- as.numeric(rep(log(pop_vector[i]), nrow(did_subset)))
        
        didreg <- glm(cases ~ treat + time + did + offset(offset_vector), 
                      data = did_subset, family = poisson(link = "log"))
        
        did_counterfactual <- did_subset
        did_counterfactual$did <- 0
        
        yhat <- predict(didreg, newdata = did_counterfactual, type = "response")
        
        did_counterfactual <- cbind(did_counterfactual, yhat)
        did_counterfactual <- as.data.frame(subset(did_counterfactual, intervention_site == i))
        did_agg <- as.data.frame(aggregate_by_quarter(did_counterfactual, i, x_ints_real[[i]]))
        
        last_row <- tail(did_agg, 1)
        
        IE <- calculate_ie(df = did_agg, est_col = 'yhat', case_col = 'cases')
        
        did.IE <- rbind(did.IE, data.frame(IE = IE[1]))
        
        did.total <- rbind(did.total, did_agg)
        agg.did <- rbind(agg.did, last_row)
      }
      
      # Aggregate IE calculation for 2x2 DiD
      IE <- calculate_ie(df = agg.did, est_col = 'est.cs', case_col = 'cases.cs')
      did.IE <- rbind(did.IE, data.frame(IE = IE[[1]]))
      
      did.IE.final <- format(did.IE)
      
      ###############################################################################
      ###################### Regression Discontinuity Design ########################
      ###############################################################################
      rdd.agg <- data.frame()
      rdd.IE <- data.frame()
      sim_rdd <- data.frame()
      
      for (i in i.sites) {
        intervention.subs <- subset(interventions.df, interventions.df$intervention_site == i)
        
        offset_vector <- as.numeric(rep(log(pop_vector[i]), nrow(intervention.subs)))
        
        panel.rdd.model <- glm(cases ~ cutoff + rel.time + I((rel.time)^2) + cutoff:rel.time + cutoff:I((rel.time)^2) + offset(offset_vector), 
                               family = quasipoisson(link = "log"), data = intervention.subs)
        
        rdd_counterfactual <- intervention.subs
        rdd_counterfactual$cutoff <- 0
        
        yhat <- predict(panel.rdd.model, newdata = rdd_counterfactual, type = "response")
        
        rdd_counterfactual <- cbind(rdd_counterfactual, yhat)
        
        sim_rdd_counterfactual <- as.data.frame(subset(rdd_counterfactual))
        sim_rdd_agg <- as.data.frame(aggregate_by_quarter(sim_rdd_counterfactual, i, x_ints_real[[i]]))
        
        last_row <- tail(sim_rdd_agg, 1)
        
        IE <- calculate_ie(df = sim_rdd_agg, est_col = 'yhat', case_col = 'cases')
        
        rdd.IE <- rbind(rdd.IE, data.frame(IE = IE[1]))
        
        sim_rdd <- rbind(sim_rdd, sim_rdd_agg)
        rdd.agg <- rbind(rdd.agg, last_row)
      }
      
      # Aggregate IE calculation for rdd
      IE <- calculate_ie(df = rdd.agg, est_col = 'est.cs', case_col = 'cases.cs')
      rdd.IE <- rbind(rdd.IE, data.frame(IE = IE[[1]]))
      
      rdd.IE.final <- format(rdd.IE)
      
      ##############################################################
      ###################### Final Tables ##########################
      ##############################################################
      i.sites.df <- c(i.sites, "Agg.")
      
      final_df <- cbind(data.frame(c(i.sites, "Agg.")), pp.IE.final, rdd.IE.final, did.IE.final, pack.did.final, scm.ie, cscm.ie, true_ie)
      colnames(final_df) = c("Site", "Pre-post", "RDD", "DiD (2x2)", "Panel DiD", "SCM", "cSCM", "True")
      
      final_df <- final_df[3,]
      final_df_t <- t(final_df[-1])
      final_df_numeric <- as.numeric(final_df_t)
      
      true_ie <- final_df_numeric[[7]]
      percentage_differences <- ((true_ie - final_df_numeric) / true_ie) * 100
      
      percentage_differences <- percentage_differences[-7]
      deviation_df <- rbind(deviation_df, percentage_differences)
      
      colnames(deviation_df) = c("pre-post", "rdd", "2x2", "panel", "scm", "cscm")
      
      saveRDS(list(deviation_df = deviation_df,
                   final_df = final_df), 
              file = file.path(results_dir, sprintf("results_size_%d_iter_%d.rds", size, iter)))
      
      rm(offset_vector, final_df_t, dengue_attgt, didreg, sim.stats, sim.sites, deviation_df, long.sim, pack.did, pack.did.final, 
         panel.rdd.model,dengue.gs, pp.model, interventions.df, rdd.agg, rdd.IE, rate.sim.sites, sim.isites.non, poissioned_Ih,
         sim_rdd, did.IE, did.IE.final, agg.did, did.total, sim_pp, pp.IE, agg.pp, sim.sc.ie)
      gc()
    },
    error=function(...) NULL)
                                                           }
}

stopCluster(cl)