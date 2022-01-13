# 00-hello.R
# Loading data and other preliminaries
# Kevin Chen
# November 2, 2021

library(here)
library(boxr); box_auth()
library(lubridate)
library(data.table)

# Load helper functions, themes, etc
source("~/headRs/00-my-theme.R")

exposure.lag <- 15

# Load full data
if (!"cohort_analytic" %in% ls()) {
	cohort_analytic <- box_read(880479268725)
	cohort_analytic <- cohort_analytic[,-"exposure.lag"]
	cohort_analytic[,(c("straight", "soluble", "synthetic")) := lapply(
		list(straight, soluble, synthetic), function(x) {
			return(shift(x, get("exposure.lag", envir = .GlobalEnv), 0))
		}),
		studyno]
	cohort_analytic[,(paste0("cum_", c("straight", "soluble", "synthetic"))) := lapply(
		list(straight, soluble, synthetic), function(x) {
			return(cumsum(x))
		}),
		studyno]
}

cohort_analytic[,`:=`(Plant = factor(plant))]
cohort_analytic[,`:=`(plant = c(1:3)[which.max(table(Plant))]), studyno]
cohort_analytic[, plant := as.numeric(plant)]
			
cohort_analytic <- cohort_analytic[
	# plant %in% c(1,2) &
	wh == 1 & nohist == 0 & yin.gm >= 1938]

