source("../aux/plotFuncs.R")

data <- read.csv("../data/data_table.csv")
data.from.exp <- data[data$exp_name == "robust_rsd" & data$mask_name == "obj_mask", ]
# data.from.exp <- data.from.exp[data.from.exp$obs_name != "JF" & data.from.exp$obs_name != "AM", ]

obs.mean.data <- aggregate(data.from.exp,
	list(img_name = data.from.exp$img_name,
		rg = data.from.exp$rg_glaven,
		by = data.from.exp$by_glaven,
		bkgd = data.from.exp$bkgd_glaven,
		illum = data.from.exp$illum_glaven,
		sub_id = data.from.exp$obs_name,
		lighter_darker = data.from.exp$hilo_glaven),
	mean)

obs.mean.data <- subset(obs.mean.data,
	select = -c(sub_id,
		obs_name,
		rg_glaven,
		by_glaven,
		bkgd_glaven,
		illum_glaven,
		hilo_glaven,
		exp_name,
		mask_name))

pop.mean.data <- aggregate(obs.mean.data,
	list(img_name = obs.mean.data$img_name,
		rg = obs.mean.data$rg,
		by = obs.mean.data$by,
		lighter_darker = obs.mean.data$lighter_darker,
		illum = obs.mean.data$illum,
		bkgd = obs.mean.data$bkgd),
	function(x) {
		c(MN = mean(x), SD = sd(x))
	})

pop.data <- pop.mean.data
obs.data <- obs.mean.data

mypalette()

maxx <- 1
maxy <- 1
adjtext <- "Adjusted R Squared =     "
lineval <- 0.5
legendloc <- "topleft"
pdf("../figures/robust_rsd_comp.pdf", height = 3, width = 9)
mypar(mfrow = c(1, 3))
for (i in c("l", "m", "s")) {
  pop_rmc_x <- pop.data[[paste0("rmc_", i, "_glaven")]][, 1]
  pop_rmc_y <- pop.data[[paste0("rmc_", i, "_match")]][, 1]
  pop_rsd_x <- pop.data[[paste0("tau_general_", i, "_glaven")]][, 1]
  pop_rsd_y <- pop.data[[paste0("tau_general_", i, "_filter")]][, 1]
  pop_rmc_sd <- pop.data[[paste0("rmc_", i, "_match")]][, 2]
  pop_rsd_sd <- pop.data[[paste0("tau_general_", i, "_filter")]][, 2]
  
  obs_rmc_x <- obs.data[[paste0("rmc_", i, "_glaven")]]
  obs_rmc_y <- obs.data[[paste0("rmc_", i, "_match")]]
  obs_rsd_x <- obs.data[[paste0("tau_general_", i, "_glaven")]]
  obs_rsd_y <- obs.data[[paste0("tau_general_", i, "_filter")]]
  
  myplot(pop_rmc_x,
         pop_rmc_y,
         xlim = c(0, maxx),
         ylim = c(0, maxy),
         col = 1,
         pch = 2,
  )
  mypoints(pop_rsd_x,
           pop_rsd_y,
           col = 6,
           pch = 6
  )
  my_V_ebar(pop_rmc_x, pop_rmc_y, pop_rmc_sd, col = 1)
  my_V_ebar(pop_rsd_x, pop_rsd_y, pop_rsd_sd, col = 6)
  abline(0, 1)
  abline(r_rmc <- lm(obs_rmc_y ~ obs_rmc_x), lty = 2, col = 1)
  abline(r_rsd <- lm(obs_rsd_y ~ obs_rsd_x), lty = 2, col = 6)
  myaxis(1, at = seq(0, maxx, maxx / 5))
  myaxis(2, at = seq(0, maxy, maxy / 5))
  
  rmc_rsq <- round(summary(r_rmc)$adj.r.squared, digits = 2)
  rsd_rsq <- round(summary(r_rsd)$adj.r.squared, digits = 2)
  mymtext(adjtext, line = lineval)
  mytitle(substitute(list(x, phantom(y)), list(x = rmc_rsq, y = rsd_rsq)),
          xlab = paste0("Ratio ", toupper(i), " - Image"),
          ylab = paste0("Ratio ", toupper(i), " - Filter Settings"),
          col.main = 1
  )
  mytitle(substitute(list(phantom(x), y), list(x = rmc_rsq, y = rsd_rsq)), col.main = 6)
  
  if (i == "l") {
    mylegend(legendloc,
             c("Mean Cone Exc.", "Robust Std. Dev. Cone Exc."),
             pch = c(2, 6),
             col = c(1, 6)
    )
  }
}
dev.off()

# let's also plot out the data once for each observer
obs.mean.data <- aggregate(data.from.exp,
	list(img_name = data.from.exp$img_name,
		rg = data.from.exp$rg_glaven,
		by = data.from.exp$by_glaven,
		bkgd = data.from.exp$bkgd_glaven,
		illum = data.from.exp$illum_glaven,
		sub_id = data.from.exp$obs_name,
		lighter_darker = data.from.exp$hilo_glaven),
	function(x) {
		c(MN = mean(x), SD = sd(x))
	})

obs.mean.data <- subset(obs.mean.data,
	select = -c(obs_name,
		rg_glaven,
		by_glaven,
		bkgd_glaven,
		illum_glaven,
		hilo_glaven,
		exp_name,
		mask_name))

for (obs in unique(obs.mean.data$sub_id)) {
	pdf(paste0("../figures/robust_rsd_comp_", obs, ".pdf"), height = 3, width = 9)
	mypar(mfrow = c(1, 3))

	subdata <- obs.mean.data[obs.mean.data$sub_id == obs, ]
	print(paste0(obs, " ", subdata$qualityc))
	for (i in c("l", "m", "s")) {
		obs_rmc_x <- subdata[[paste0("rmc_", i, "_glaven")]][, 1]
		obs_rmc_y <- subdata[[paste0("rmc_", i, "_match")]][, 1]
		obs_rmc_sd <- subdata[[paste0("rmc_", i, "_match")]][, 2]
		obs_rsd_x <- subdata[[paste0("tau_general_", i, "_glaven")]][, 1]
		obs_rsd_y <- subdata[[paste0("tau_general_", i, "_filter")]][, 1]
		obs_rsd_sd <- subdata[[paste0("tau_general_", i, "_filter")]][, 2]
  
		myplot(obs_rmc_x,
			obs_rmc_y,
			xlim = c(0, maxx),
			ylim = c(0, maxy),
			col = 1,
			pch = 2,
		)
		mypoints(obs_rsd_x,
			obs_rsd_y,
			col = 6,
			pch = 6
		)
		my_V_ebar(obs_rmc_x, obs_rmc_y, obs_rmc_sd, col = 1)
		my_V_ebar(obs_rsd_x, obs_rsd_y, obs_rsd_sd, col = 6)
		abline(0, 1)
		abline(r_rmc <- lm(obs_rmc_y ~ obs_rmc_x), lty = 2, col = 1)
		abline(r_rsd <- lm(obs_rsd_y ~ obs_rsd_x), lty = 2, col = 6)
		myaxis(1, at = seq(0, maxx, maxx / 5))
		myaxis(2, at = seq(0, maxy, maxy / 5))
  
		rmc_rsq <- round(summary(r_rmc)$adj.r.squared, digits = 2)
		rsd_rsq <- round(summary(r_rsd)$adj.r.squared, digits = 2)
		mymtext(adjtext, line = lineval)
		mytitle(substitute(list(x, phantom(y)), list(x = rmc_rsq, y = rsd_rsq)),
			xlab = paste0("Ratio ", toupper(i), " - Image"),
			ylab = paste0("Ratio ", toupper(i), " - Filter Settings"),
			col.main = 1
		)
		mytitle(substitute(list(phantom(x), y), list(x = rmc_rsq, y = rsd_rsq)), col.main = 6)
  
		if (i == "l") {
			mylegend(legendloc,
				c("Mean Cone Exc.", "Robust Std. Dev. Cone Exc."),
				pch = c(2, 6),
				col = c(1, 6)
			)
		}
	}
	dev.off()
}

# let's save out some of the matches that the observers made and see what they look like
#for (obs in unique(obs.mean.data$sub_id)) {
#	subdata <- obs.mean.data[obs.mean.data$sub_id == obs, ]
#
#	for (ic in 1:nrow(subdata)) {
#		ld <- subdata$ld_or_r[ic]
#		rg <- subdata$rg_or_g[ic]
#		yv <- subdata$by_or_b[ic]
#
#		system(paste0("voronoi_filters -saveImg=true -printStats=false",
#						" -red=../base_stimuli/spectra/munsell_red_EXTREEEMMMEE.spd",
#						" -green=../base_stimuli/spectra/munsell_green_EXTREEEMMMEE.spd",
#						" -blue=../base_stimuli/spectra/munsell_blue_EXTREEEMMMEE.spd",
#						" -yellow=../base_stimuli/spectra/munsell_yellow_EXTREEEMMMEE.spd",
#						" -ld=", as.character(ld),
#						" -rg=", as.character(rg),
#						" -by=", as.character(yv),
#						" -spectraFile=../calibration/exp5.spectra.csv",
#						" -chromaFile=../calibration/exp5.chroma.csv",
#						" -gammaFile=../calibration/exp5.gamma.csv"))
#
#		l_or_d <- ""
#		if (subdata$lighter_darker[ic] == "darker") {
#			l_or_d <- "1"
#		} else {
#			l_or_d <- "2"
#		}
#
#		g_or_r <- ""
#		if (subdata$illum[ic] == "red") {
#			g_or_r <- "3"
#		} else {
#			g_or_r <- "4"
#		}
#
#		system(paste0("mv",
#				" voronoi_filter_ld_*.png",
#				" ",
#				"../images/exp6/voronoi_matches/",
#				"voronoi_filter_subID_", obs,
#				"_illum_", subdata$illum[ic],
#				"_rg_", subdata$rg[ic],
#				"_by_", subdata$by[ic],
#				"_ld_", l_or_d,
#				"_bkgd_", subdata$bkgd[ic],
#				".png"))
#	}
#}