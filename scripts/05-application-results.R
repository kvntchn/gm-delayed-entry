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

# Get script from github
library(RCurl)
eval(parse(text = getURL(
	"https://raw.githubusercontent.com/kvntchn/HeadRs/master/00-my-theme.R",
	ssl.verifypeer = F
)))

str_a0 <- readRDS(here("resources", paste0("straight0_gm", ".rds")))
str_a1 <- readRDS(here("resources", paste0("straight1_gm", ".rds")))
sol_a0 <- readRDS(here("resources", paste0("soluble0_gm", ".rds")))
sol_a1 <- readRDS(here("resources", paste0("soluble1_gm", ".rds")))
syn_a0 <- readRDS(here("resources", paste0("synthetic0_gm", ".rds")))
syn_a1 <- readRDS(here("resources", paste0("synthetic1_gm", ".rds")))

gm.surv <- rbindlist(list(
	rbindlist(list(str_a0[[2]], str_a1[[2]]), idcol = "rule"),
	rbindlist(list(sol_a0[[2]], sol_a1[[2]]), idcol = "rule"),
	rbindlist(list(syn_a0[[2]], syn_a1[[2]]), idcol = "rule")
), idcol = "mwf")

gm.w <- rbindlist(list(
	rbindlist(list(str_a0[[1]], str_a1[[1]]), idcol = "rule"),
	rbindlist(list(sol_a0[[1]], sol_a1[[1]]), idcol = "rule"),
	rbindlist(list(syn_a0[[1]], syn_a1[[1]]), idcol = "rule")
), idcol = "mwf")


gm.surv[,`:=`(
	mwf = {
		mwf.labs <- c("Straight", "Soluble", "Synthetic")
		factor(mwf.labs[mwf], mwf.labs)
		},
	rule = c("$a_0$", "$a_1$")[rule]
)]

gm.w[,`:=`(
	mwf = {
		mwf.labs <- c("Straight", "Soluble", "Synthetic")
		factor(mwf.labs[mwf], mwf.labs)
		},
	rule = c("$a_0$", "$a_1$")[rule]
)]

gm.w <- rbindlist(list(
	gm.w[,.(mwf, rule, t, weight = cum.IPAW, estimator = "WKM")],
	gm.w[R == 1,.(mwf, rule, t, weight = cum.IPAW, estimator = "AWKM")]))


gm_w.tab <- dcast(gm.w[,.(
	summary = as.vector(summary(weight)),
	stat = {
		stat <- c("Minimum", "Q1", "Median", "Mean", "Q3", "Maximum")
		factor(stat, stat)
	}), .(t, mwf, rule, estimator)],
	estimator + mwf + rule + t ~ stat,
	value.var = "summary"
)

# gm.surv %>% melt(
# 	measure.vars = c("WKM", "AWKM"),
# 	variable.name = "Estimator",
# 	value.name = "Survival"
# ) %>% ggplot(aes(
# 	x = t, y = Survival,
# 	lty = Estimator,
# 	color = rule
# )) + geom_step() + facet_grid(mwf ~ .) +
# 	labs(
# 		y = "Cancer-free survival",
# 		x = "Years since hire") +
# 	theme_bw() + theme(
# 		legend.title = element_blank(),
# 		legend.position = 'bottom'
# 	) -> gm.surv.ggplot
# 
# gm.surv.ggplot
# 
# gm_w.tab[estimator == "WKM"] %>% ggplot(
# 	aes(x = t,
# 			y = Median,
# 			ymin = Minimum,
# 			ymax = Maximum,
# 			color = rule,
# 			fill = rule)
# ) + geom_line() +
# 	geom_ribbon(alpha = 0.2, color = NA) +
# 	facet_grid(mwf ~ .) +
# 	coord_cartesian(ylim = c(0, 50)) +
# 	theme_bw() +
# 	labs(y = "Median cumulative weight",
# 			 x = "Years since hire") +
# 	theme(legend.position = "bottom",
# 				legend.title = element_blank()) -> gm.w.ggplot
# 
# gm.w.ggplot
# 
# tikz(here("reports/resources", "gm_results.tex"),
# 		 width = 4.5, height = 5, standAlone = T)
# gm.surv.ggplot
# dev.off()
# lualatex("gm_results\\.tex", here("reports/resources"))
# 
# 
# tikz(here("reports/resources", "gm_weights.tex"),
# 		 width = 4.5, height = 5, standAlone = T)
# gm.w.ggplot
# dev.off()
# lualatex("gm_weights\\.tex", here("reports/resources"))
