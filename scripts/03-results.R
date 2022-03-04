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

# Get script from github
library(RCurl)
eval(parse(text = getURL(
	"https://raw.githubusercontent.com/kvntchn/HeadRs/master/00-my-theme.R",
	ssl.verifypeer = F
)))

# Get package for running estimators
# devtools::install_github('osofr/stremr')
library(stremr)

results.dir <- "resources"

truth_a0 <- readRDS(here(results.dir, "truth_a0.rds"))
truth_a1 <- readRDS(here(results.dir, "truth_a1.rds"))

data.frame(truth_a0[[5]][,2],
					 truth_a1[[5]][,2])

a0 <- rbindlist(
	lapply(1:5, function(scenario = 1) {
		melt(rbindlist(readRDS(here(results.dir, paste0("a0_scenario", scenario, ".rds")))
		)[,.(WKM = mean(WKM), AWKM = mean(AWKM)), t],
		measure.vars = c("WKM", "AWKM"),
		variable.name = "Estimator",
		value.name = "Survival")
	}),
	idcol = "scenario"
)
a1 <- rbindlist(
	lapply(1:5, function(scenario) {
		melt(rbindlist(readRDS(here(results.dir, paste0("a1_scenario", scenario, ".rds")))
		)[,.(WKM = mean(WKM), AWKM = mean(AWKM)), t],
		measure.vars = c("WKM", "AWKM"),
		variable.name = "Estimator",
		value.name = "Survival")
	}),
	idcol = "scenario"
)


sim_results <- rbindlist(list(
	rbindlist(lapply(list(rbindlist(truth_a0, idcol = "scenario"),
												rbindlist(truth_a1, idcol = "scenario")),
									 melt,
									 measure.vars = "Truth",
									 variable.name = "Estimator",
									 value.name = "Survival"),
						idcol = "rule"),
	rbindlist(list(a0, a1),
						idcol = "rule")
))

sim_results[,`:=`(
	scenario = paste0("Scenario ", scenario),
	rule = paste0(c("Never exposed", "Always exposed while at work")[rule])
		# paste0("Rule ", c("$a_0$", "$a_1$")[rule])
)]

sim_results %>% ggplot(aes(
	x = t, y = Survival, lty = Estimator,
	color = rule,
)) + geom_step() + facet_grid(scenario ~ .) +
	# coord_cartesian(ylim = c(0.8, 1)) +
	labs(x = "Years since hire", y = "Cancer-free survival") +
	mytheme + theme(
		legend.title = element_blank(),
		legend.position = "bottom",
		legend.box = "vertical",
		legend.margin = margin()
	) -> results.ggplot

# results.ggplot
# 
# # saveRDS(sim_results, here(paste0("reports/", results.dir), "sim_results.rds"))
# tikz(here(paste0("reports/", results.dir), paste0("results", ".tex")),
# 		 width = 4.5,
# 		 height = 7.75,
# 		 # height = 7.5/4,
# 		 standAlone = T)
# print(results.ggplot)
# dev.off()
# 
# lualatex("results\\.tex", here(paste0("reports/", results.dir)))
