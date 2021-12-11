# 02-analysis.R
# Analyzing simulated data
# Kevin Chen
# December 6, 2021

library(here)
library(lubridate)
library(tidyverse)
library(magrittr)
library(data.table)
library(survival)

# Get code for simulation
source(here("scripts", "01-simulation.R"))

# Get package for running estimators
# devtools::install_github('osofr/stremr')
library(stremr)

get_estimates <- function(
	scenario = 1,
	a = 0,
	g_a = "A ~ W + shiftA",
	g_d = "D ~ W + barA",
	stratify = list(D = paste0("t == ", 1:20),
									A = paste0("t == ", 1:20))
) {
	
	# message("Simulating...")
	sim.dt <- sim_data(scenario, a = a)
	sim.dt[,`:=`(
		shiftA = shift(A, 1, 0),
		barA = shift(as.integer(cumsum(A) > 0), 1, 0)
		), id]
	
	# Get rules
	sim.dt[,`:=`(A_0 = 0, A_1 = N)]
	
	sim.dt <- sim.dt[t <= last_t_obs]
	
	# Make stremr object
	sim.stremr <- importData(
		sim.dt,
		ID = "id",
		t_name = "t",
		covars = c("R", "W", "S", "H", "N", "barA"),
		CENS = "D",
		TRT = "A",
		OUTCOME = "Y_obs")
	
	# Fit PS
	# message("Fitting PS")
	sim.stremr <- fitPropensity(
		sim.stremr,
		gform_TRT = g_a,
		gform_CENS = g_d,
		stratify_TRT = stratify,
		stratify_CENS = stratify,
	)
	
	w <- getIPWeights(
		sim.stremr,
		intervened_TRT = paste0("A_", a))
	
	wkm <- w %>% survNPMSM(sim.stremr) %$% estimates
	
	# nrow(w)
	# nrow(sim.dt)
	# Get rid of rows according to `obs` for Aalen filter
	awkm <- w[sim.dt$obs == 1,] %>%
		survNPMSM(sim.stremr) %$%
		estimates
	
	return(
		data.table(
			t = 1:20,
			# KM = wkm$St.KM,
			WKM = wkm$St.NPMSM,
			AWKM = awkm$St.NPMSM)
	)
}


# km %>% melt(
# 	measure.vars = names(.)[-1],
# 	variable.name = "estimator",
# 	value = "survival"
# ) %>% ggplot(aes(x = t, y = survival, lty = estimator)) +
# 	geom_step() +
# 	theme_bw()

get_truth <- function(scenario, a, n = 5e5) {
	dt <- sim_data(scenario, n = n, a = a, d = 0, long = F)
	km <- survfit(Surv(last_t, status) ~ 1, data = dt[,.(last_t, status)])
	return(data.table(
		t = 1:20,
		Truth = c(rep(1, 20 - length(km$surv)), km$surv)
		))}

# truth_a0 <- lapply(1:5, get_truth, a = 0, n = 5e4)
# saveRDS(truth_a0, here("resources", "truth_a0.rds"))
# truth_a1 <- lapply(1:5, get_truth, a = 1, n = 5e4)
# saveRDS(truth_a1, here("resources", "truth_a1.rds"))

M <- 250
pb <- txtProgressBar(min = 0, max = M, style = 3)
message("\nScenario ", scenario, ", a = 0")
a0 <- lapply(1:M, function(i) {
	setTxtProgressBar(pb, i)
	return(get_estimates(scenario))})
saveRDS(a0, here("resources", paste0("a0_scenario", scenario, ".rds")))

pb <- txtProgressBar(min = 0, max = M, style = 3)
message("\nScenario ", scenario, ", a = 1")
a1 <- lapply(1:M, function(i) {
	setTxtProgressBar(pb, i)
	return(get_estimates(scenario, a = 1))})
saveRDS(a1, here("resources", paste0("a1_scenario", scenario, ".rds")))
