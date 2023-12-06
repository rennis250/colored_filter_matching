source("../aux/plotFuncs.R")

data <- read.csv("../data/data_table.csv")
data_dark <- read.csv("../data/data_table_dark_exc.csv")

expn <- "patch_app"

do <- data[data$mask_name == "obj_mask" & data$exp_name == expn, ]
dr <- data[data$mask_name == "refl_mask" & data$exp_name == expn, ]
dowoh <- data[data$mask_name == "obj_wout_high_mask" & data$exp_name == expn, ]
ddark <- data_dark[data_dark$mask_name == "obj_mask" & data_dark$exp_name == expn, ]

agg_data <- list()
for (cc in 1:4) {
	if (cc == 1) {
		d <- do
	} else if (cc == 2) {
		d <- dr
	} else if (cc == 3) {
		d <- dowoh
	} else {
		d <- ddark
	}

	obs_mean_data <- aggregate(d, list(body = d$body_glaven,
										illum = d$illum_glaven,
										obs = d$obs_name,
										hilo = d$hilo_glaven),
								mean)
	obs_mean_data <- subset(obs_mean_data, select = -c(obs,
													obs_name,
													body_glaven,
													illum_glaven,
													hilo_glaven,
													exp_name,
													mask_name))
	pop_mean_data <- aggregate(obs_mean_data, list(body_glaven = obs_mean_data$body,
												illum_glaven = obs_mean_data$illum,
												hilo_glaven = obs_mean_data$hilo),
								function(x) {
									c(MN = mean(x), SD = sd(x))
								})
	pop_mean_data <- subset(pop_mean_data, select = -c(body, illum, hilo))

	agg_data[[cc]] <- pop_mean_data
}

# let's look at how RMC, RSD, and tau_general change between wout_high and normal mask
rmc_change_high <- 100*c(mean(agg_data[[3]]$rmc_l_glaven[, 1]/agg_data[[1]]$rmc_l_glaven[, 1]),
			mean(agg_data[[3]]$rmc_m_glaven[, 1]/agg_data[[1]]$rmc_m_glaven[, 1]),
			mean(agg_data[[3]]$rmc_s_glaven[, 1]/agg_data[[1]]$rmc_s_glaven[, 1]))

rsd_change_high <- 100*c(mean(agg_data[[3]]$rsd_l_glaven[, 1]/agg_data[[1]]$rsd_l_glaven[, 1]),
			mean(agg_data[[3]]$rsd_m_glaven[, 1]/agg_data[[1]]$rsd_m_glaven[, 1]),
			mean(agg_data[[3]]$rsd_s_glaven[, 1]/agg_data[[1]]$rsd_s_glaven[, 1]))

tau_change_high <- 100*c(mean(agg_data[[3]]$tau_general_l_glaven[, 1]/agg_data[[1]]$tau_general_l_glaven[, 1]),
			mean(agg_data[[3]]$tau_general_m_glaven[, 1]/agg_data[[1]]$tau_general_m_glaven[, 1]),
			mean(agg_data[[3]]$tau_general_s_glaven[, 1]/agg_data[[1]]$tau_general_s_glaven[, 1]))

# let's look at how RMC, RSD, and tau_general change between refl and normal mask
rmc_change_refl <- 100*c(mean(agg_data[[2]]$rmc_l_glaven[, 1]/agg_data[[1]]$rmc_l_glaven[, 1]),
			mean(agg_data[[2]]$rmc_m_glaven[, 1]/agg_data[[1]]$rmc_m_glaven[, 1]),
			mean(agg_data[[2]]$rmc_s_glaven[, 1]/agg_data[[1]]$rmc_s_glaven[, 1]))

rsd_change_refl <- 100*c(mean(agg_data[[2]]$rsd_l_glaven[, 1]/agg_data[[1]]$rsd_l_glaven[, 1]),
			mean(agg_data[[2]]$rsd_m_glaven[, 1]/agg_data[[1]]$rsd_m_glaven[, 1]),
			mean(agg_data[[2]]$rsd_s_glaven[, 1]/agg_data[[1]]$rsd_s_glaven[, 1]))

tau_change_refl <- 100*c(mean(agg_data[[2]]$tau_general_l_glaven[, 1]/agg_data[[1]]$tau_general_l_glaven[, 1]),
			mean(agg_data[[2]]$tau_general_m_glaven[, 1]/agg_data[[1]]$tau_general_m_glaven[, 1]),
			mean(agg_data[[2]]$tau_general_s_glaven[, 1]/agg_data[[1]]$tau_general_s_glaven[, 1]))

# let's look at how RMC, RSD, and tau_general change when dark pixels are excluded
rmc_change_dark <- 100*c(mean(agg_data[[4]]$rmc_l_glaven[, 1]/agg_data[[1]]$rmc_l_glaven[, 1]),
			mean(agg_data[[4]]$rmc_m_glaven[, 1]/agg_data[[1]]$rmc_m_glaven[, 1]),
			mean(agg_data[[4]]$rmc_s_glaven[, 1]/agg_data[[1]]$rmc_s_glaven[, 1]))

rsd_change_dark <- 100*c(mean(agg_data[[4]]$rsd_l_glaven[, 1]/agg_data[[1]]$rsd_l_glaven[, 1]),
			mean(agg_data[[4]]$rsd_m_glaven[, 1]/agg_data[[1]]$rsd_m_glaven[, 1]),
			mean(agg_data[[4]]$rsd_s_glaven[, 1]/agg_data[[1]]$rsd_s_glaven[, 1]))

tau_change_dark <- 100*c(mean(agg_data[[4]]$tau_general_l_glaven[, 1]/agg_data[[1]]$tau_general_l_glaven[, 1]),
			mean(agg_data[[4]]$tau_general_m_glaven[, 1]/agg_data[[1]]$tau_general_m_glaven[, 1]),
			mean(agg_data[[4]]$tau_general_s_glaven[, 1]/agg_data[[1]]$tau_general_s_glaven[, 1]))

print("rmc_change_high -> ")
print(round(abs(100 - mean(rmc_change_high)), digits=3))
print(round(sd(rmc_change_high), digits=3))
print("rsd_change_high -> ")
print(round(abs(100 - mean(rsd_change_high)), digits=3))
print(round(sd(rsd_change_high), digits=3))
print("tau_change_high -> ")
print(round(abs(100 - mean(tau_change_high)), digits=3))
print(round(sd(tau_change_high), digits=3))

print("rmc_change_refl -> ")
print(round(abs(100 - mean(rmc_change_refl)), digits=3))
print(round(sd(rmc_change_refl), digits=3))
print("rsd_change_refl -> ")
print(round(abs(100 - mean(rsd_change_refl)), digits=3))
print(round(sd(rsd_change_refl), digits=3))
print("tau_change_refl -> ")
print(round(abs(100 - mean(tau_change_refl)), digits=3))
print(round(sd(tau_change_refl), digits=3))

print("rmc_change_dark -> ")
print(round(abs(100 - mean(rmc_change_dark)), digits=3))
print(round(sd(rmc_change_dark), digits=3))
print("rsd_change_dark -> ")
print(round(abs(100 - mean(rsd_change_dark)), digits=3))
print(round(sd(rsd_change_dark), digits=3))
print("tau_change_dark -> ")
print(round(abs(100 - mean(tau_change_dark)), digits=3))
print(round(sd(tau_change_dark), digits=3))

# [1] "rmc_change_high -> "
# [1] 0.048
# [1] 0.013
# [1] "rsd_change_high -> "
# [1] 0.157
# [1] 0.067
# [1] "tau_change_high -> "
# [1] 0.105
# [1] 0.012
# [1] "rmc_change_refl -> "
# [1] 0.393
# [1] 0.02
# [1] "rsd_change_refl -> "
# [1] 0.137
# [1] 0.026
# [1] "tau_change_refl -> "
# [1] 0.568
# [1] 0.024