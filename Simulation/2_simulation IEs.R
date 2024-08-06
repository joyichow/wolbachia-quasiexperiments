library(data.table)

setwd("/Users/joyichow/Desktop/wolb_interventions/simulation_results")

results_dir <- "simulation_results"
control_sizes <- c(10, 15, 20, 25, 30, 35, 40) # Control sizes
iterations <- 1000

combined_final_dfs <- list()
combined_deviation_dfs <- list()

for (size in control_sizes) {
  iteration_final_dfs <- list()
  iteration_deviation_dfs <- list()
  
  for (iter in 1:iterations) {
    
    file_name <- file.path(sprintf("results_size_%d_iter_%d.rds", size, iter))
    
    tryCatch({
      results <- readRDS(file_name)
      
      iteration_final_dfs[[iter]] <- results$final_df
      iteration_deviation_dfs[[iter]] <- results$deviation_df
    }, error = function(e) {
      cat("Error reading file:", file_name, "\n", sep=" ")

      iteration_final_dfs[[iter]] <- NULL
      iteration_deviation_dfs[[iter]] <- NULL
    })
  }

  iteration_final_dfs <- Filter(Negate(is.null), iteration_final_dfs)
  iteration_deviation_dfs <- Filter(Negate(is.null), iteration_deviation_dfs)

  if (length(iteration_final_dfs) > 0) {
    combined_final_dfs[[as.character(size)]] <- do.call(rbind, iteration_final_dfs)
  }
  
  if (length(iteration_deviation_dfs) > 0) {
    combined_deviation_dfs[[as.character(size)]] <- do.call(rbind, iteration_deviation_dfs)
  }
}

means_list <- list()
sds_list <- list()

for (size in names(combined_deviation_dfs)) {
  means_list[[size]] <- colMeans(combined_deviation_dfs[[size]], na.rm = TRUE)
  sds_list[[size]] <- apply(combined_deviation_dfs[[size]], 2, sd, na.rm = TRUE)
}

formatted_results <- data.frame()

for (size in names(means_list)) {
  means <- means_list[[size]]
  sds <- sds_list[[size]]
  
  formatted_values <- mapply(function(mean, sd) paste0(round(mean, 2), " (", round(sd, 2), ")"), mean = means, sd = sds)
  formatted_vector <- setNames(as.vector(formatted_values), names(means))
  formatted_results <- rbind(formatted_results, formatted_vector)
}

rownames(formatted_results) <- names(means_list)
colnames(formatted_results) <- c("pre-post", "rdd", "2x2", "panel", "scm", "cscm")

library(writexl)
write_xlsx(formatted_results, "/Users/joyichow/Desktop/wolb_interventions/transfer.xlsx")
