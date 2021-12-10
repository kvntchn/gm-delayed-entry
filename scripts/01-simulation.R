# 01-simulation.R
# Simulating data
# Kevin Chen
# November 3, 2021

library(here)
library(lubridate)
library(data.table)

# Parameters
sim_param <- readRDS(here::here("resources", "sim_param.rds"))
rownames(sim_param) <- sim_param[,1]
sim_param <- sim_param[,-1]

sim_data <- function(
	scenario = 1,
	n = 50000,
	full = T,
	K = 20,
	d = NULL,
	a = NULL,
	long = T,
	param = sim_param) {
	
	env.current <- environment()
	
	param <- as.data.frame(t(param))[scenario,]
	
	U <- as.data.frame(
		sapply(paste(
			rep(c("U_H", "U_N", "U_A", "U_D", "U_Y"), each = K),
			1:K, sep = "."), function(x) {
				runif(n)
			}))
	
	R <- runif(n, 0, 20)
	W <- rbinom(n, 1, param$p_W)
	S <- rbinom(n, 1, param$p_S)
	H.k <- H.1 <- rbinom(n, 1, param$p_H)
	N.k <- N.1 <- rbinom(n, 1, plogis(
		param$beta_0N +
			param$beta_WN * W +
			param$beta_HN * H.k +
			U[,paste0("U_N.", 1)]
	))
	if (is.null(d)) {
		D.k <- D.1 <- rbinom(n, 1, plogis(
			param$beta_0D +
				param$beta_WD * W +
				U[, paste0("U_D.", 1)]
		))
	} else {
		D.k <- D.1 <- rep(d, n)
	}
	if (is.null(a)) {
		A.1 <- rbinom(n, 1, plogis(
			param$beta_0A + param$beta_WA * W + U$U_A.1
		))
	} else {
		A.1 <- rep(a, n)
	}
	A.1[N.k == 0] <- 0
	A.k <- A.1
	Y.k <- Y.1 <- rbinom(n, 1, plogis(
		param$beta_0Y +
			param$beta_WY * W +
			param$beta_AY * A.1 +
			param$beta_SY * S * A.1 +
			param$beta_HY * H.k +
			U$U_Y.1
	))
	
	# k <- 2
	# k <- k + 1
	for (k in 2:K) {
		H.k[H.k != 1] <- rbinom(sum(H.k != 1), 1, param$p_H)
		assign(paste0("H.", k), H.k, envir = env.current)
		N.k[N.k != 0] <- rbinom(sum(N.k != 0), 1, plogis(
			param$beta_0N +
				param$beta_WN * W[N.k != 0] +
				param$beta_HN * H.k[N.k != 0] +
				param$beta_AN * A.k[N.k != 0] +
				U[N.k != 0, paste0("U_N.", k)]
		))
		assign(paste0("N.", k), N.k, envir = env.current)
		ybar <- apply(sapply(1:(k - 1), function(k) {get(paste0("Y.", k))}), 1, sum)
		if (is.null(d)) {
			D.k[D.k != 1] <- rbinom(sum(D.k != 1), 1, plogis(
				param$beta_0D +
					param$beta_WD * W[D.k != 1] +
					param$beta_barAD * apply(
						sapply(1:(k - 1), function(k) {get(paste0("A.", k))}), 1, function(x) {sum(x) > 0})[D.k != 1] +
					# sapply(1:(k - 1), function(k) {get(paste0("A.", k))}), 1, function(x) {sum(x)})[D.k != 1] +
					param$beta_barYD * ybar[D.k != 1] +
					U[D.k != 1, paste0("U_D.", k)]
			))
		} else {
			D.k <- rep(d, n)
		}
		assign(paste0("D.", k), D.k, envir = env.current)
		if (is.null(a)) {
			A.k <- rbinom(n, 1, plogis(
				param$beta_AA * A.k + U[, paste0("U_A.", k)]
			))
		} else {
			A.k <- rep(a, n)
		}
		A.k[N.k == 0] <- 0
		assign(paste0("A.", k), A.k, envir = env.current)
		barA <- apply(sapply(1:(k - 1), function(k) {get(paste0("A.", k))}), 1, function(x) {sum(x) > 0})
		# barA <- apply(sapply(1:(k - 1), function(k) {get(paste0("A.", k))}), 1, function(x) {sum(x)})
		Y.k[Y.k != 1] <- rbinom(sum(Y.k != 1), 1, plogis(
			param$beta_0Y +
				param$beta_WY * W[Y.k != 1] +
				param$beta_AY * A.k[Y.k != 1] +
				param$beta_barAY * barA[Y.k != 1] +
				param$beta_SY * S[Y.k != 1] * apply(
					cbind(A.k, barA), 1, function(x) {sum(x) > 0})[Y.k != 1] +
				param$beta_HY * H.k[Y.k != 1] +
				U[Y.k != 1, paste0("U_Y.", k)]
		))
		assign(paste0("Y.", k), Y.k, envir = env.current)
	}
	
	if (!long) {
		dt <- as.data.frame(sapply(c("R", "W", "S", paste(
			c("H", "N", "D", "A", "Y"),
			rep(1:K, each = 6), sep = ".")),
			get, envir = env.current))
		
		y.which <- grep("Y\\.", names(dt))
		d.which <- grep("D\\.", names(dt))
		
		dt$last_t <- apply(dt, 1, function(x) {
			min(K, which(x[c(y.which, d.which)] == 1))
		})
		
		dt$status <- apply(dt, 1, function(x) {
			max(x[y.which])
		})
		
		return(as.data.table(dt))
	} else {
		dt <- as.data.frame(sapply(c("R", "W", "S", paste(
			c("H", "N", "D", "A", "Y"),
			rep(1:K, each = 6), sep = ".")),
			get, envir = env.current))
		dt$id <- 1:nrow(dt)
		
		ldt <- melt(
			as.data.table(dt),
			measure.vars = lapply(
				c("H", "N", "D", "A", "Y"), function(x) {
					paste(x, 1:K, sep = ".")
				}),
			value.name = c("H", "N", "D", "A", "Y"),
			variable.name = "t",
			variable.factor = F
		)
		
		ldt[,t := as.numeric(t)]
		
		setorder(ldt, id, t)
		
		ldt[, obs := as.integer(t + 1 >= R), id]
		ldt[,Y_obs := {
			if (sum(Y[1:floor(R[1])]) > 0) {
				Y_obs <- rep(0, .N)
			} else {
				Y_obs <- Y
				Y_obs[t < R] <- 0}
			Y_obs
		},id]
		
		ldt[,`:=`(
			last_t = min(t[Y == 1], t[D == 1], K),
			last_t_obs = min(t[Y_obs == 1], t[D == 1], K)
		), id]
		
		return(ldt[,.(id, t, R, W, S, H, N, D, A, Y, Y_obs, obs, last_t, last_t_obs)])
	}
}
