# 04-application.R
# Analyzing simulated data
# Kevin Chen
# December 6, 2021

library(here)
library(lubridate)
library(tidyverse)
library(magrittr)
library(data.table)
library(survival)

library(boxr); box_auth()
if (!"cohort_analytic" %in% ls()) {
	source(here::here("scripts", "00-hello.R"))
}

# Get package for running estimators
# devtools::install_github('osofr/stremr')
library(stremr)

years_fu <- 40
cohort_analytic <- cohort_analytic[
	year < year(yin) + years_fu &
		nohist == 0 &
		wh == 1]

# Covariate engineering
cohort_analytic[,`:=`(
	Straight = as.integer(straight > 0),
	Soluble = as.integer(soluble > 0),
	Synthetic = as.integer(synthetic > 0),
	R = as.integer(year >= 1985),
	Plant = factor(plant),
	Race = factor(race),
	Sex = factor(sex),
	Age = cut(age.year2/365, c(
		min(age.year2/365),
		quantile(age.year2[canc_co == 1]/365, seq(0.2, 0.8, 0.2)),
		max(age.year2/365)),
		include.lowest = T, dig.lab = 4),
	Hire = cut(yin.gm, c(
		min(yin.gm),
		quantile(yin.gm[canc_co == 1], seq(0.2, 0.8, 0.2)),
		max(yin.gm)),
		include.lowest = T, dig.lab = 4),
	N = as.integer(year <= year(yout)),
	off = {
		off <- off.gan + off.han + off.san
		off[year > year(yout)] <- 0
		off[year >= 1995] <- 0
		off[is.na(off)] <- 0},
	duration = time_length(difftime(year2, yin), 'year')
)]
# Lag some covariates
cohort_analytic[,`:=`(
	off = shift(cumsum(off), 15, 1),
	duration = shift(duration, 15, 1)
), studyno]
# Cut by distribution of cases after lag
cohort_analytic[, `:=`(
	Off = cut(off, unique(c(
		-Inf,
		0,
		quantile(na.omit(off[off != 0 & canc_co == 1]), seq(0.2, 0.8, 0.2)),
		max(off))),
		include.lowest = T),
	Duration = cut(duration, c(
		min(duration),
		quantile(duration[canc_co == 1], seq(0.2, 0.8, 0.2)),
		max(duration)
	),
	include.lowest = T)
)]
# sapply(cohort_analytic[canc_co == 1,.(Age, Hire, Off, Duration)], table, useNA = "always")
# Lag exposure; employment status
cohort_analytic[,`:=`(
	Straight = shift(Straight, 15, 0),
	Soluble = shift(Soluble, 15, 0),
	Synthetic = shift(Synthetic, 15, 0),
	N = shift(N, 15, 0)
), studyno]
# Get exposure history
cohort_analytic[,`:=`(
	barStraight = shift(Straight, 16, 0),
	barSoluble = shift(Soluble, 16, 0),
	barSynthetic = shift(Synthetic, 16, 0),
	N = shift(N, 15, 0)
), studyno]

get_estimates_gm <- function(
	a = 0,
	mwf = "Straight",
	Y = "canc_co",
	dt = copy(cohort_analytic),
	g_a = "A ~ Plant + Race + Sex + Hire + Age + Off + Straight + Soluble + Synthetic + barStraight + barSoluble + barSynthetic",
	g_d = "D ~ Plant + Race + Sex + Hire + Age + Off + Straight + Soluble + Synthetic + barStraight + barSoluble + barSynthetic",
	stratify = list(D = paste0("t == ", 1:years_fu),
									A = paste0("t == ", 1:years_fu))
) {
	
	# message("Simulating...")
	
	g_a <- gsub(paste0(" \\+ ", mwf), "", g_a)
	g_d <- gsub(paste0(" \\+ ", mwf), "", g_d)
	
	# Variable engineering
	dt[, `:=`(
		A = get(mwf),
		Y = get(Y),
		D = `All causes`
	)]
	dt[, `:=`(
		# barA = shift(as.integer(cumsum(A) > 0), 1, fill = 0),
		t = 1:.N,
		n_rows = .N
	), studyno]
	
	# Get rules
	dt[,`:=`(A_0 = 0, A_1 = N)]
	
	# Get rid of rows after event
	dt[,`:=`(
		last_row = min(t[Y == 1], max(t), na.rm = T)
	), studyno]
	
	dt <- dt[t <= last_row]
	
	# Make stremr object
	dt.stremr <- importData(
		dt,
		ID = "studyno",
		t_name = "t",
		covars = c("R", "Plant", "Race", "Sex", "Hire", "Age", "Off", "Straight", "Soluble", "Synthetic", "barStraight", "barSoluble", "barSynthetic"),
		CENS = "D",
		TRT = "A",
		OUTCOME = "Y")
	
	# Fit PS
	# message("Fitting PS")
	dt.stremr <- fitPropensity(
		dt.stremr,
		gform_TRT = g_a,
		gform_CENS = g_d,
		stratify_TRT = stratify,
		stratify_CENS = stratify,
	)
	
	w <- getIPWeights(
		dt.stremr,
		intervened_TRT = paste0("A_", a))
	
	wkm <- w %>% survNPMSM(dt.stremr) %$% estimates
	
	# Get rid of rows according to `obs` for Aalen filter
	awkm <- w[dt$R,] %>%
		survNPMSM(dt.stremr) %$%
		estimates
	
	return(
		data.table(
			t = 1:max(dt$t),
			# KM = wkm$St.KM,
			WKM = wkm$St.NPMSM,
			AWKM = awkm$St.NPMSM)
	)
}

# data.table(
# 			t = 1:max(dt$t),
# 			# KM = wkm$St.KM,
# 			WKM = wkm$St.NPMSM,
# 			AWKM = awkm$St.NPMSM) %>% melt(
# 	measure.vars = names(.)[-1],
# 	variable.name = "estimator",
# 	value = "survival"
# ) %>% ggplot(aes(x = t, y = survival, lty = estimator)) +
# 	geom_step() +
# 	theme_bw()

if (!is.null(mwf)) {
	message("\n", mwf, ", a = 0")
	a0 <- get_estimates_gm(a = 0,
												 mwf = mwf)
	saveRDS(a0, here("resources", paste0(tolower(mwf), "0_gm", ".rds")))
	
	message("\n", mwf, ", a = 1")
	a1 <- get_estimates_gm(a = 1,
												 mwf = mwf)
	saveRDS(a1, here("resources", paste0(tolower(mwf), "1_gm", ".rds")))
}