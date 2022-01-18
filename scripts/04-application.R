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

if (!"cohort_analytic" %in% ls()) {
	library(boxr); box_auth()
	source(here::here("scripts", "00-hello.R"))
	message("Note. Exposure lagged ", exposure.lag, " years")
}

# Get package for running estimators
# devtools::install_github('osofr/stremr')
library(stremr)

years_fu <- 55
gm.dt <- cohort_analytic[year >= year(yin) & year < (year(yin) + years_fu) & (
	year(yin) <= (1995 + exposure.lag - years_fu) | year <= (1994 + exposure.lag))]
gm.dt <- gm.dt[plant %in% 1:2]
# gm.dt <- gm.dt[sex == "M"]
gm.dt[, `:=`(
	canc = ifelse(
		canc_corec == 1 | canc_es == 1 | canc_st == 1, 1,
		ifelse(canc_corec == 2 | canc_es == 2 | canc_st == 2, 2, 0))
)]
# gm.dt[, `:=`(canc = canc_nhl)]
# gm.dt[, canc := canc_first]

message(nrow(gm.dt[canc == 1]), " cases in ", n_distinct(gm.dt$studyno), " subjects")

gm.dt[, `:=`(
	other_canc = ifelse(canc_first == 1 | canc_first == 2 & canc == 0, 1, 0)
)]

# How much are people exposed?
gm.dt[year >= year(yin) + 15 & year <= year(yout) + 15,.(
	studyno,
	straight = straight > 0.05,
	soluble = soluble > 0.05,
	synthetic = synthetic > 0.05
)][,.(straight = mean(straight),
			soluble = mean(soluble),
			synthetic = mean(synthetic)), studyno][,-1] %>% sapply(mean)

mwf.quantiles <- gm.dt[canc == 1,.(
	straight = quantile(cum_straight[cum_straight > 0], seq(0.25, 1 - 0.25, 0.25)),
	soluble  = quantile(cum_soluble[cum_soluble > 0.05], seq(0.25, 1 - 0.25, 0.25)),
	synthetic = quantile(cum_synthetic[cum_synthetic > 0], seq(1/3, 1 - 1/3, 1/3))
)]
mwf.quantiles

# Covariate engineering
gm.dt[,`:=`(
	Straight = as.integer(cum_straight > 0),
	Soluble = as.integer(cum_soluble > 0.05),
	Synthetic = as.integer(cum_synthetic > 0),
	cum_Straight = cut(
		cum_straight,
		unique(c(-Inf, 0, mwf.quantiles$straight, max(cum_straight))),
		include.lowest = T),
	cum_Soluble = cut(
		cum_soluble,
		unique(c(-Inf, 0.05, mwf.quantiles$soluble, max(cum_soluble))),
		include.lowest = T),
	cum_Synthetic = cut(
		cum_synthetic,
		unique(c(-Inf, 0, mwf.quantiles$synthetic, max(cum_synthetic))),
		include.lowest = T),
	R = as.integer(year >= 1985 | (year >= 1973 & plant %in% 1:2)),
	# R = as.integer(year >= 1985),
	Plant = factor(plant),
	Race = factor(race),
	Sex = factor(sex),
	Age = cut(age.year2/365, c(
		min(age.year2/365),
		quantile(age.year2[canc == 1]/365, seq(1/5, 1 - 1/5, 1/5)),
		max(age.year2/365)),
		include.lowest = T, dig.lab = 4),
	Hire = cut(yin.gm, c(
		min(yin.gm),
		quantile(yin.gm[canc == 1], seq(1/5, 1 - 1/5, 1/5)),
		max(yin.gm)),
		include.lowest = T, dig.lab = 4),
	# N = as.integer(year <= year(yout) + exposure.lag & year >= year(yin) + exposure.lag),
	N = as.integer(year <= year(yout) + exposure.lag),
	Off = {
		off <- off.gan + off.han + off.san
		off[year > year(yout)] <- 0
		off[year >= 1995] <- 0
		off[is.na(off)] <- 0
		cut(off, unique(c(
			-Inf, 0,
			quantile(off[canc == 1], seq(1/5, 1 - 1/5, 1/5)),
			max(off))))
	},
	Duration = {
		duration <- time_length(difftime(year2, yin), 'year')
		cut(duration, c(
			min(duration),
			quantile(duration[canc == 1], seq(0.5, 1 - 0.5, 0.5)),
			max(duration)
		),
		include.lowest = T)
	}
)]

if (class(gm.dt$Straight) == "logical" | class(gm.dt$Straight) == "integer") {
	gm.dt[, Any:=as.integer(Straight | Soluble | Synthetic)]
} else {
	gm.dt[, Any := Straight + Soluble + Synthetic]
} 

# Lag time off and some cumulative exposure summaries for (PS modeling)
gm.dt[, Off := shift(Off, get("exposure.lag", envir = .GlobalEnv), 1)]
if (class(gm.dt$cum_Straight) != "factor") {stop("Check class of cum_Straight")}
gm.dt[,(paste0("cum_", c("Straight", "Soluble", "Synthetic"))) := lapply(
	paste0("cum_", c("Straight", "Soluble", "Synthetic")), function(x) {
		shift(get(x),1, 1)
	}
), studyno]

setorder(gm.dt, studyno, year)
# gm.dt <- gm.dt[canc != 2]
gm.dt[, `:=`(
	t = 1:.N,
	n_rows = .N
), studyno]

# sapply(gm.dt[canc == 1,.(Age, Hire, Off, Duration)], table, useNA = "always")
# table(gm.dt[t >= 10, .(year, yin, yout, N, Straight), studyno]$Straight)
# table(gm.dt[N == 1,.(comply = sum(Soluble[-(1:4)] == N[-(1:4)]) == .N - 4), studyno]$comply)

t.quant <- quantile(gm.dt[canc == 1, t], seq(0, 1, 0.2))
t.quant <- unique(c(1, t.quant, max(gm.dt$t) + 1))
seq(1, years_fu, 5)

get_estimates_gm <- function(
	a = 1,
	mwf = "Straight",
	Y = "canc",
	dt = copy(gm.dt),
	g_a = paste0("A ~ Plant + Race + Sex + Hire + Age + Off + N + cum_Straight + cum_Soluble + cum_Synthetic + Straight + Soluble + Synthetic"),
	g_d = paste0("D ~ Plant + Race + Sex + Hire + Age + Off + N + cum_Straight + cum_Soluble + cum_Synthetic + Straight + Soluble + Synthetic"),
	stratify =
		list(D = paste("t ==", 1:years_fu),
				 A = paste("t ==", 1:years_fu)),
	# list(D = paste0("t >= ", t.quant[-length(t.quant)], " & t < ", na.omit(shift(t.quant, -1))),
	# 		 A = paste0("t >= ", t.quant[-length(t.quant)], " & t < ", na.omit(shift(t.quant, -1)))),
	trunc_weights = 5 * 10^2
) {
	
	g_a <- gsub(paste0(" \\+ ", mwf), "", g_a, ignore.case = T)
	
	if (mwf == "Any") {
		g_a <- gsub(paste(paste0(" \\+ ", c("Straight", "Soluble", "Synthetic")), collapse = "|"), "", g_a, ignore.case = T)
	}
	
	# Variable engineering
	dt[, `:=`(
		A = get(mwf),
		Y = get(Y),
		D = {
			D = `All causes`
			D[canc == 1 | canc == 2] <- 0
			D}
	)]
	dt[, `:=`(
		t = 1:.N,
		n_rows = .N
	), studyno]
	
	# Get rules
	dt[,`:=`(
		A_0 = A,
		A_1 = 0)]
	# mean(dt[,.(followed = sum(A == N)), studyno]$followed)
	# mean(dt[,.(followed = sum(A == 0)), studyno]$followed)
	
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
		covars = c("R", "Plant", "Race", "Sex", "Hire", "Age", "Off",
							 "Straight", "Soluble", "Synthetic",
							 "cum_Straight", "cum_Soluble", "cum_Synthetic", "Any", "N"),
		CENS = "D",
		TRT = c("A"),
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
	
	# Get post-intervention cumulative exposure history
	if (a == 1) {
		get_data(dt.stremr)[,(c(paste0("cum_", mwf), mwf)):=list(
			1, 0)]
		
		if (mwf == "Any") {
			get_data(dt.stremr)[,(paste0(rep(c("cum_", ""), each = 3), c("Straight", "Soluble", "Synthetic"))):=list(
				1, 1, 1, 0, 0, 0)]
		}
	}
	
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

# mwf <- "Synthetic"
# Point estimate ###
if (!is.null(mwf)) {
	message("\n", mwf, ", a = 0")
	a0 <- get_estimates_gm(a = 0,
												 dt = copy(gm.dt),
												 mwf = mwf)
	saveRDS(a0, here("resources", paste0(tolower(mwf), "0_gm", ".rds")))
	
	message("\n", mwf, ", a = 1")
	a1 <- get_estimates_gm(a = 1,
												 dt = copy(gm.dt),
												 mwf = mwf)
	saveRDS(a1, here("resources", paste0(tolower(mwf), "1_gm", ".rds")))
}

# Bootstrap ###
studyno <- unique(gm.dt$studyno)
B <- 500
# set.seed(1240)
# who.mcmc <- lapply(1:B, function(b = 1) {
# 	who.mcmc <- data.table(studyno = sample(studyno, replace = T))
# 	who.mcmc$id <- 1:nrow(who.mcmc)
# 	return(who.mcmc)
# 	})
# saveRDS(who.mcmc, here("resources", "who_mcmc.rds"))
who.mcmc <- readRDS(here("resources", "who_mcmc.rds"))
if (!is.null(mwf)) {
	start <- Sys.time()
	message("BS started at ", start, "\n")
	# Progress bar
	pb <- txtProgressBar(min = 0, max = B, style = 3)
	bs <- list()
	for (b in 1:B) {
		bs.dt <- merge(
			who.mcmc[[b]],
			gm.dt[studyno %in% who.mcmc[[b]]$studyno,
						c("studyno", "year", "canc", "other_canc", "All causes",
							"R", "Plant", "Race", "Sex", "Hire", "Age", "Off",
							"Straight", "Soluble", "Synthetic",
							"cum_Straight", "cum_Soluble", "cum_Synthetic", "Any", "N"),
						with = F],
			by = "studyno",
			all.x = T, allow.cartesian = T
		)
		bs.dt[,`:=`(studyno = id,
								id = studyno)]
		
		a0 <- get_estimates_gm(a = 0,
													 dt = copy(bs.dt),
													 mwf = mwf)
		
		a1 <- get_estimates_gm(a = 1,
													 dt = copy(bs.dt),
													 mwf = mwf)
		# Set progressbar
		setTxtProgressBar(pb, b)
		
		bs[[b]] <- list(a0 = a0, a1 = a1)
	}
	saveRDS(bs, here("resources", paste0(tolower(mwf), "_bs_gm", ".rds")))
}

