source("../aux/plotFuncs.R")

# let's plot data from experiment with the chosen five and their brighter counterparts

data.from.exp <- read.csv('../data/exp5_dark_pxs_excl_obs_data.csv')

# exclude me, to clean up messiness of data

data.from.exp <- data.from.exp[data.from.exp["sub_id"] != "re", ]

obs.mean.data <- aggregate(. ~ sub_id + img + img_name + illum + rg + by + lighter_darker + bkgd, data.from.exp, mean)

# look at observers that did experiment with EXTREMMEEE distributions

sub_ids <- obs.mean.data["sub_id"]
obs.mean.dataEXTREMMEEE <- obs.mean.data[sub_ids == "05" | sub_ids == "re2" | sub_ids == "08a", ]
obs.mean.dataEXTREMMEEE <- subset(obs.mean.dataEXTREMMEEE, select = -c(sub_id))
obs.mean.dataEXTREMMEEE <- subset(obs.mean.dataEXTREMMEEE, select = -c(img_name))
obs.mean.dataEXTREMMEEE <- subset(obs.mean.dataEXTREMMEEE, select = -c(illum))
obs.mean.dataEXTREMMEEE <- subset(obs.mean.dataEXTREMMEEE, select = -c(rg))
obs.mean.dataEXTREMMEEE <- subset(obs.mean.dataEXTREMMEEE, select = -c(by))
obs.mean.dataEXTREMMEEE <- subset(obs.mean.dataEXTREMMEEE, select = -c(lighter_darker))
obs.mean.dataEXTREMMEEE <- subset(obs.mean.dataEXTREMMEEE, select = -c(bkgd))
pop.mean.dataEXTREMMEEE <- aggregate(. ~ img, obs.mean.dataEXTREMMEEE, function(x) { c(MN=mean(x), SD=sd(x)) })

pdf("../figures/rmc_rsd_comp_EXTREMMEEE_dark_pxs_excl.pdf", height = 3, width = 9)
mypar(mfrow=c(1,3))
plotRMCandRSD(pop.mean.dataEXTREMMEEE, obs.mean.dataEXTREMMEEE, 3, 3, "topleft")
dev.off()


# now, let's zoom in and paint the interesting image in red

pdf("../figures/rmc_rsd_comp_ZOOM_dark_pxs_excl.pdf", height = 3, width = 9)
maxx <- 1.2
maxy <- 1.2
mypar(mfrow=c(1,3))
for (i in c("l", "m", "s")) {
	pop_rmc_x <- pop.mean.dataEXTREMMEEE[[paste0("rmc_", i)]][,1]
	pop_rmc_y <- pop.mean.dataEXTREMMEEE[[paste0("mean_ratio_", i)]][,1]
	pop_rsd_x <- pop.mean.dataEXTREMMEEE[[paste0("rsd_", i)]][,1]
	pop_rsd_y <- pop.mean.dataEXTREMMEEE[[paste0("std_ratio_", i)]][,1]
	pop_rmc_sd <- pop.mean.dataEXTREMMEEE[[paste0("mean_ratio_", i)]][,2]
	pop_rsd_sd <- pop.mean.dataEXTREMMEEE[[paste0("std_ratio_", i)]][,2]

	obs_rmc_x <- obs.mean.dataEXTREMMEEE[[paste0("rmc_", i)]]
	obs_rmc_y <- obs.mean.dataEXTREMMEEE[[paste0("mean_ratio_", i)]]
	obs_rsd_x <- obs.mean.dataEXTREMMEEE[[paste0("rsd_", i)]]
	obs_rsd_y <- obs.mean.dataEXTREMMEEE[[paste0("std_ratio_", i)]]

	red_id <- pop.mean.dataEXTREMMEEE$img_name == "mitsuba_comp_rmc_rsd_illum_4_rg_3_by_1_ld_2_bkgd_voronoi_diff_dist.png"

	myplot(pop_rmc_x,
		pop_rmc_y,
		xlim=c(0,maxx),
		ylim=c(0,maxy),
		col = 1,
		pch = 2)
	# mypoints(pop_rmc_x[red_id],
	# 	pop_rmc_y[red_id],
	# 	col = 4,
	# 	pch = 2)
	mypoints(pop_rsd_x,
		pop_rsd_y,
		col = 6,
		pch = 6)
	my_V_ebar(pop_rmc_x, pop_rmc_y, pop_rmc_sd, col=1)
	# my_V_ebar(pop_rmc_x[red_id], pop_rmc_y[red_id], pop_rmc_sd[red_id], col=4)
	my_V_ebar(pop_rsd_x, pop_rsd_y, pop_rsd_sd, col=6)
	abline(0, 1)
	abline(r_rmc <- lm(obs_rmc_y ~ obs_rmc_x), lty=2, col=1)
	abline(r_rsd <- lm(obs_rsd_y ~ obs_rsd_x), lty=2, col=6)
	myaxis(1, at = seq(0, maxx, maxx/5))
	myaxis(2, at = seq(0, maxy, maxy/5))

	rmc_rsq <- round(summary(r_rmc)$adj.r.squared, digits=2)
	rsd_rsq <- round(summary(r_rsd)$adj.r.squared, digits=2)
	mymtext("Adjusted R Squared =      ", line = 0.6)
	mytitle(substitute(list(x, phantom(y)), list(x = rmc_rsq, y = rsd_rsq)),
		xlab = paste0("Ratio ", toupper(i), " - Image"),
		ylab = paste0("Ratio ", toupper(i), " - Filter Settings"),
		col.main = 1)
	mytitle(substitute(list(phantom(x), y), list(x = rmc_rsq, y = rsd_rsq)), col.main = 6)

	if (i == "l") {
		mylegend("topleft",
			c("Mean Cone Exc.", "Sta. Dev. Cone Exc."),
			pch=c(2, 6),
			col=c(1, 6))
	}
}
dev.off()

# without the the one poorly rated image... they still both don't explain it...

non_outliers <- which(obs.mean.dataEXTREMMEEE$img_name != "mitsuba_comp_rmc_rsd_illum_4_rg_3_by_1_ld_2_bkgd_voronoi_diff_dist.png")
obs.mean.dataEXTREMMEEE.no.weird <- obs.mean.dataEXTREMMEEE[non_outliers, ]
pop.mean.dataEXTREMMEEE.no.weird <- aggregate(. ~ img + img_name + illum + rg + by + lighter_darker + bkgd, obs.mean.dataEXTREMMEEE.no.weird, function(x) { c(MN=mean(x), SD=sd(x)) })

pdf("../figures/rmc_rsd_comp_wout_poorly_rated_dist_EXTREMMMEEE_dark_pxs_excl.pdf", height = 3, width = 9)
mypar(mfrow=c(1,3))
plotRMCandRSD(pop.mean.dataEXTREMMEEE.no.weird, obs.mean.dataEXTREMMEEE.no.weird, 1.2, 1.2, "topleft", "Adjusted R Squared =      ", 0.6)
dev.off()
