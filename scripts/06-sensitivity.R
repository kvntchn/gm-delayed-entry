# 06-sensitivity.R
# Sensitivity analyses
# Kevin Chen
# December 18, 2021

library(here)
source(here("scripts", "02-analysis.R"))
# source(here("scripts", "04-application.R"))


# Use observed distribution of time to registry #####
# time_to_registry <- gm.dt[canc != 2,.(R = c(1973, 1973, 1985)[plant[1]] - yin.gm[1]), studyno]$R
# saveRDS(time_to_registry, here("resources", "time_to_registry.rds"))
time_to_registry <- readRDS(here("resources", "time_to_registry.rds"))
time_to_registry[time_to_registry <= 0] <- 0
# hist(time_to_registry, breaks = seq(0, 35, 1))
# summary(gm.dt[,.(yin.gm[1]), studyno]$V1)

# # Get truth
# truth_a0 <- lapply(1:5, get_truth, a = 0, n = 5e4,
# 									 R.rdist = function(x) {sample(time_to_registry, x, T)})
# saveRDS(truth_a0, here("resources/sensitivity", "truth_a0.rds"))
# truth_a1 <- lapply(1:5, get_truth, a = 1, n = 5e4,
# 									 R.rdist = function(x) {sample(time_to_registry, x, T)})
# saveRDS(truth_a1, here("resources/sensitivity", "truth_a1.rds"))
# 
# # Get simulation results
# M <- 250
# pb <- txtProgressBar(min = 0, max = M, style = 3)
# message("\nScenario ", scenario, ", a = 0")
# a0 <- lapply(1:M, function(i) {
# 	setTxtProgressBar(pb, i)
# 	return(get_estimates(scenario,
# 									 R.rdist = function(x) {sample(time_to_registry, x, T)}))})
# saveRDS(a0, here("resources/sensitivity", paste0("a0_scenario", scenario, ".rds")))
# 
# pb <- txtProgressBar(min = 0, max = M, style = 3)
# message("\nScenario ", scenario, ", a = 1")
# a1 <- lapply(1:M, function(i) {
# 	setTxtProgressBar(pb, i)
# 	return(get_estimates(scenario, a = 1,
# 									 R.rdist = function(x) {sample(time_to_registry, x, T)}))})
# saveRDS(a1, here("resources/sensitivity", paste0("a1_scenario", scenario, ".rds")))