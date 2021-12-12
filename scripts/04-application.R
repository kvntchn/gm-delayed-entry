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

years_fu <- 35
gm.dt <- cohort_analytic[year >= year(yin) & year < year(yin) + years_fu & (year(yin) <= 1995 - years_fu | year <= 1994)]
gm.dt <- gm.dt[plant %in% 1:2]
gm.dt[, canc := ifelse(canc_corec == 1 | canc_es == 1 | canc_st == 1, 1, 
										ifelse(canc_corec == 2 | canc_es == 2 | canc_st == 2, 2, 0))]

# Covariate engineering
gm.dt[,`:=`(
	Straight = as.integer(straight > median(straight[year <= year(yout)])),
	Soluble = as.integer(soluble > median(soluble[year <= year(yout)])),
	Synthetic = as.integer(synthetic > median(synthetic[year <= year(yout)])),
	R = as.integer(year >= 1985 | (year >= 1973 & plant %in% 1:2)),
	Plant = factor(plant),
	Race = factor(race),
	Sex = factor(sex),
	Age = cut(age.year2/365, c(
		min(age.year2/365),
		quantile(age.year2[canc == 1]/365, seq(0.2, 0.8, 0.2)),
		max(age.year2/365)),
		include.lowest = T, dig.lab = 4),
	Hire = cut(yin.gm, c(
		min(yin.gm),
		quantile(yin.gm[canc == 1], seq(0.2, 0.8, 0.2)),
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
gm.dt[,`:=`(
	# Straight = shift(Straight, 10, 0),
	# Soluble = shift(Soluble, 10, 0),
	# Synthetic = shift(Synthetic, 10, 0),
	bar1Straight = shift(Straight, 1, 0),
	bar1Soluble = shift(Soluble, 1, 0),
	bar1Synthetic = shift(Synthetic, 1, 0),
	bar2Straight = shift(Straight, 2, 0),
	bar2Soluble = shift(Soluble, 2, 0),
	bar2Synthetic = shift(Synthetic, 2, 0),
	bar3Straight = shift(Straight, 3, 0),
	bar3Soluble = shift(Soluble, 3, 0),
	bar3Synthetic = shift(Synthetic, 3, 0),
	bar4Straight = shift(Straight, 4, 0),
	bar4Soluble = shift(Soluble, 4, 0),
	bar4Synthetic = shift(Synthetic, 4, 0),
	bar5Straight = shift(Straight, 5, 0),
	bar5Soluble = shift(Soluble, 5, 0),
	bar5Synthetic = shift(Synthetic, 5, 0),
	bar6Straight = shift(Straight, 6, 0),
	bar6Soluble = shift(Soluble, 6, 0),
	bar6Synthetic = shift(Synthetic, 6, 0)
	# N = shift(N, 10, 0),
	# off = shift(cumsum(off), 10, 1),
	# duration = shift(duration, 10, 1)
), studyno]

# Cut by distribution of cases after lag
gm.dt[, `:=`(
	Off = cut(off, unique(c(
		-Inf,
		0,
		quantile(na.omit(off[off != 0 & canc == 1]), seq(0.2, 0.8, 0.2)),
		max(off))),
		include.lowest = T),
	Duration = cut(duration, c(
		min(duration),
		quantile(duration[canc == 1], seq(0.2, 0.8, 0.2)),
		max(duration)
	),
	include.lowest = T)
)]

setorder(gm.dt, studyno, year)
gm.dt[, `:=`(
	# barA = shift(as.integer(cumsum(A) > 0), 1, fill = 0),
	t = 1:.N,
	n_rows = .N
), studyno]

# sapply(gm.dt[canc == 1,.(Age, Hire, Off, Duration)], table, useNA = "always")
# table(gm.dt[t >= 10, .(year, yin, yout, N, Straight), studyno]$Straight)
# table(gm.dt[N == 1,.(comply = sum(Soluble[-(1:4)] == N[-(1:4)]) == .N - 4), studyno]$comply)

t.quant <- quantile(gm.dt[canc == 1, t], seq(0, 1, 0.1))
gm.dt[year > year(yout) & Straight > 0, .(studyno, year, yout15, straight)]$su
# seq(1, years_fu, 5)

get_estimates_gm <- function(
	a = 1,
	mwf = "Straight",
	Y = "canc",
	dt = copy(gm.dt),
	g_a = paste0("A ~ Plant + Race + Sex + Hire + Age + Off + N + Straight + Soluble + Synthetic + ",
							 paste(as.vector(sapply(1:3, function(x) {
							 		paste0("bar", x, c("Straight", "Soluble", "Synthetic"))})),
							 			collapse = " + ")),
	g_d = paste0("A ~ Plant + Race + Sex + Hire + Age + Off + N + Straight + Soluble + Synthetic + ",
							 paste(as.vector(sapply(1:6, function(x) {
							 		paste0("bar", x, c("Straight", "Soluble", "Synthetic"))})),
							 			collapse = " + ")),
	stratify =
		# list(D = paste("t ==", 1:years_fu),
		# 		 A = paste("t ==", 1:years_fu))
		list(D = paste0("t >= ", t.quant[-length(t.quant)], " & t < ", na.omit(shift(t.quant, -1))),
				 A = paste0("t >= ", t.quant[-length(t.quant)], " & t < ", na.omit(shift(t.quant, -1)))),
	trunc_weights = 10^3
) {
	
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
	dt[,`:=`(A_0 = A * 0.75, A_1 = A)]
	# mean(dt[,.(followed = A == N), studyno]$followed)
	
	# Get rid of rows after event
	dt[,`:=`(
		last_row = min(t[Y == 1], max(t), na.rm = T)
	), studyno]
	
	dt <- dt[t <= last_row]
	# dt[t == 1 & R == 1,.(studyno, t, R)]
	
	# Make stremr object
	dt.stremr <- importData(
		dt,
		ID = "studyno",
		t_name = "t",
		covars = c("R", "Plant", "Race", "Sex", "Hire", "Age", "Off", "Straight", "Soluble", "Synthetic", "N",
							 as.vector(sapply(1:6, function(x) {
							 		paste0("bar", x, c("Straight", "Soluble", "Synthetic"))}))),
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
	
	wkm <- w %>% survNPMSM(dt.stremr, trunc_weights = trunc_weights) %$% estimates
	
	# Get rid of rows according to `obs` for Aalen filter
	awkm <- w[as.logical(dt$R),] %>%
		survNPMSM(dt.stremr, trunc_weights = trunc_weights) %$%
		estimates
	
	# format(object.size(w), 'MB')
	# format(object.size(w[,.(studyno, t, cum.IPAW, cum.stab.P, wt.by.t)]), "MB")
	
	return(
		list(
			w[,.(studyno, t, cum.IPAW, cum.stab.P, wt.by.t, R = dt$R)],
			data.table(
				t = 1:max(dt$t),
				# KM = wkm$St.KM,
				WKM = wkm$St.NPMSM,
				AWKM = awkm$St.NPMSM)
		))
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
