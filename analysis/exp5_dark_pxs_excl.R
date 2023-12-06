# we are going to check if excluding dark pixels makes RMC or RSD do better for the darker stimuli,
# so we need to get the data from both with and without those dark pixels to compare
# data.from.exp.dark.excl <- read.csv('../data/exp5_dark_pxs_excl_obs_data.csv')
data <- read.csv("../data/data_table_dark_exc.csv")
data.from.exp.dark.excl <- data[data$exp_name == "rmc_rsd_comp" & data$mask_name == "obj_mask", ]

# data.from.exp <- read.csv('../data/exp5_obs_data.csv')
data <- read.csv("../data/data_table.csv")
data.from.exp <- data[data$exp_name == "rmc_rsd_comp" & data$mask_name == "obj_mask", ]

# exclude me, to clean up messiness of data
data.from.exp.dark.excl <- data.from.exp.dark.excl[data.from.exp.dark.excl["obs_name"] != "re", ]
data.from.exp <- data.from.exp[data.from.exp["obs_name"] != "re", ]

# obs.mean.data.dark.excl <- aggregate(. ~ sub_id + img + img_name + illum + rg + by + lighter_darker + bkgd, data.from.exp.dark.excl, mean)
obs.mean.data.dark.excl <- aggregate(data.from.exp.dark.excl, list(img_name = data.from.exp.dark.excl$img_name, rg = data.from.exp.dark.excl$rg_glaven, by = data.from.exp.dark.excl$by_glaven, bkgd = data.from.exp.dark.excl$bkgd_glaven, illum = data.from.exp.dark.excl$illum_glaven, sub_id = data.from.exp.dark.excl$obs_name, lighter_darker = data.from.exp.dark.excl$hilo_glaven), mean)

# obs.mean.data <- aggregate(. ~ sub_id + img + img_name + illum + rg + by + lighter_darker + bkgd, data.from.exp, mean)
obs.mean.data <- aggregate(data.from.exp, list(img_name = data.from.exp$img_name, rg = data.from.exp$rg_glaven, by = data.from.exp$by_glaven, bkgd = data.from.exp$bkgd_glaven, illum = data.from.exp$illum_glaven, sub_id = data.from.exp$obs_name, lighter_darker = data.from.exp$hilo_glaven), mean)

# look at observers that did experiment with EXTREMMEEE distributions
# first, compute the means for the dark excluded data
sub_ids <- obs.mean.data.dark.excl["sub_id"]
obs.mean.dataEXTREMMEEE.dark.excl <- obs.mean.data.dark.excl[sub_ids == "05" | sub_ids == "re2" | sub_ids == "08a", ]
obs.mean.dataEXTREMMEEE.dark.excl <- subset(obs.mean.dataEXTREMMEEE.dark.excl, select = -c(obs_name, rg_glaven, by_glaven, bkgd_glaven, illum_glaven, hilo_glaven, exp_name, mask_name))
obs.mean.dataEXTREMMEEE.dark.excl <- subset(obs.mean.dataEXTREMMEEE.dark.excl, select = -c(rg, by, bkgd, illum, lighter_darker))

# pop.mean.dataEXTREMMEEE.dark.excl <- aggregate(. ~ img, obs.mean.dataEXTREMMEEE.dark.excl, function(x) { c(MN=mean(x), SD=sd(x)) })
pop.mean.dataEXTREMMEEE.dark.excl <- aggregate(obs.mean.dataEXTREMMEEE.dark.excl, list(img_name = obs.mean.dataEXTREMMEEE.dark.excl$img_name), function(x) {
  c(MN = mean(x), SD = sd(x))
})

# and then, compute them for the full data
sub_ids <- obs.mean.data["sub_id"]
obs.mean.dataEXTREMMEEE <- obs.mean.data[sub_ids == "05" | sub_ids == "re2" | sub_ids == "08a", ]
obs.mean.dataEXTREMMEEE <- subset(obs.mean.dataEXTREMMEEE, select = -c(obs_name, rg_glaven, by_glaven, bkgd_glaven, illum_glaven, hilo_glaven, exp_name, mask_name))
obs.mean.dataEXTREMMEEE <- subset(obs.mean.dataEXTREMMEEE, select = -c(rg, by, bkgd, illum, lighter_darker))

# pop.mean.dataEXTREMMEEE <- aggregate(. ~ img, obs.mean.dataEXTREMMEEE, function(x) { c(MN=mean(x), SD=sd(x)) })
pop.mean.dataEXTREMMEEE <- aggregate(obs.mean.dataEXTREMMEEE, list(img_name = obs.mean.dataEXTREMMEEE$img_name), function(x) {
  c(MN = mean(x), SD = sd(x))
})

# we really just want to see if the RMC/RSD for these two images improves when we exclude the dark pixels.
# all of the other images/objects are bright enough to not warrant rejecting pixels.

# so, grab their data for the dark excluded images first
idxs <- pop.mean.dataEXTREMMEEE.dark.excl$img_name == "mitsuba_comp_rmc_rsd_illum_4_rg_3_by_1_ld_2_bkgd_voronoi_diff_dist.png" | pop.mean.dataEXTREMMEEE.dark.excl$img_name == "mitsuba_comp_rmc_rsd_illum_4_rg_3_by_4_ld_1_bkgd_voronoi_yellow.png"
pop.mean.thetwoimages.dark.excl <- pop.mean.dataEXTREMMEEE.dark.excl[idxs, ]

# and, now for the full data
idxs <- pop.mean.dataEXTREMMEEE$img_name == "mitsuba_comp_rmc_rsd_illum_4_rg_3_by_1_ld_2_bkgd_voronoi_diff_dist.png" | pop.mean.dataEXTREMMEEE$img_name == "mitsuba_comp_rmc_rsd_illum_4_rg_3_by_4_ld_1_bkgd_voronoi_yellow.png"
pop.mean.thetwoimages <- pop.mean.dataEXTREMMEEE[idxs, ]

maxx <- 10.0
maxy <- 0.4
legendloc <- "topright"

mypalette()
pdf("../figures/rmc_rsd_comp_EXTREMMEEE_dark_pxs_excl.pdf", height = 3, width = 9)
mypar(mfrow=c(1,3))
for (i in c("l", "m", "s")) {
	pop_rmc_x.wodark <- pop.mean.thetwoimages.dark.excl[[paste0("rmc_", i, "_glaven")]][, 1]
	pop_rmc_y.wodark <- pop.mean.thetwoimages.dark.excl[[paste0("rmc_", i, "_match")]][, 1]
	pop_rmc_sd.wodark <- pop.mean.thetwoimages.dark.excl[[paste0("rmc_", i, "_match")]][, 2]

	pop_rmc_x.second <- pop.mean.thetwoimages[[paste0("rmc_", i, "_glaven")]][, 1]
	pop_rmc_y.second <- pop.mean.thetwoimages[[paste0("rmc_", i, "_match")]][, 1]
	pop_rmc_sd.second <- pop.mean.thetwoimages[[paste0("rmc_", i, "_match")]][, 2]
	
	pop_rsd_x.wodark <- pop.mean.thetwoimages.dark.excl[[paste0("rsd_", i, "_glaven")]][, 1]
	pop_rsd_y.wodark <- pop.mean.thetwoimages.dark.excl[[paste0("rsd_", i, "_filter")]][, 1]
	pop_rsd_sd.wodark <- pop.mean.thetwoimages.dark.excl[[paste0("rsd_", i, "_filter")]][, 2]
	
	pop_rsd_x.second <- pop.mean.thetwoimages[[paste0("rsd_", i, "_glaven")]][, 1]
	pop_rsd_y.second <- pop.mean.thetwoimages[[paste0("rsd_", i, "_filter")]][, 1]
	pop_rsd_sd.second <- pop.mean.thetwoimages[[paste0("rsd_", i, "_filter")]][, 2]

	myplot(pop_rmc_x.wodark,
	       pop_rmc_y.wodark,
		xlim=c(0,maxx),
		ylim=c(0,maxy),
		col = 1,
		pch = 1,)
	mypoints(pop_rmc_x.second,
		pop_rmc_y.second,
		xlim=c(0,maxx),
		ylim=c(0,maxy),
		col = 1,
		pch = 1,)
	my_V_ebar(pop_rmc_x.wodark, pop_rmc_y.wodark, pop_rmc_sd.wodark, col=1)
	my_V_ebar(pop_rmc_x.second, pop_rmc_y.second, pop_rmc_sd.second, col=1)
	
	mypoints(pop_rsd_x.wodark,
	         pop_rsd_y.wodark,
	         xlim=c(0,maxx),
	         ylim=c(0,maxy),
	         col = 6,
	         pch = 1,)
	mypoints(pop_rsd_x.second,
	         pop_rsd_y.second,
	         xlim=c(0,maxx),
	         ylim=c(0,maxy),
	         col = 6,
	         pch = 1,)
	my_V_ebar(pop_rsd_x.wodark, pop_rsd_y.wodark, pop_rsd_sd.wodark, col=6)
	my_V_ebar(pop_rsd_x.second, pop_rsd_y.second, pop_rsd_sd.second, col=6)
	
	abline(0, 1)
	myaxis(1, at = seq(0, maxx, maxx/5))
	myaxis(2, at = seq(0, maxy, maxy/5))
	
	arrows(pop_rmc_x.second, pop_rmc_y.second, pop_rmc_x.wodark, pop_rmc_y.wodark, length = 0.12, lwd = 1, col = 4)
	arrows(pop_rsd_x.second, pop_rsd_y.second, pop_rsd_x.wodark, pop_rsd_y.wodark, length = 0.12, lwd = 1, col = 4)

	# mymtext(adjtext, line = lineval)
	mytitle(xlab = paste0("Ratio ", toupper(i), " - Image"),
		ylab = paste0("Ratio ", toupper(i), " - Filter Settings"),
		col.main = 1)

	if (i == "l") {
		mylegend(legendloc,
			c("Mean Cone Exc.", "Std. Dev. Cone Exc."),
			pch=c(1, 1),
			col=c(1, 6))
	}
}
dev.off()
