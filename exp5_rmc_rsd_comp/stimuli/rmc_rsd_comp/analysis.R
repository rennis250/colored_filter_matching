source("../../../aux/plotFuncs.R")

data <- read.csv('scenes/eizo_chrom_stats.csv')

data$color <- NA

for(ic in 1:length(data$img_name)) {
	fparts <- strsplit(as.character(data$img_name[ic]), "_")
	if (fparts[[1]][6] == "1") {
		data$color[ic] <- 6
	} else if (fparts[[1]][6] == "2") {
		data$color[ic] <- 7
	} else if (fparts[[1]][6] == "3") {
		data$color[ic] <- 4
	} else if (fparts[[1]][6] == "4") {
		data$color[ic] <- 5
	}
}

data$pch <- NA

for(ic in 1:length(data$img_name)) {
	fparts <- strsplit(as.character(data$img_name[ic]), "_")
	bkgd <- strsplit(fparts[[1]][length(fparts[[1]])], "p")
	l_or_d <- fparts[[1]][12]
	if (l_or_d == "1") {
		if (bkgd[[1]][1] == "blue.") {
			data$pch[ic] <- 0
		} else if (bkgd[[1]][1] == "yellow.") {
			data$pch[ic] <- 1
		} else if (bkgd[[1]][1] == "red.") {
			data$pch[ic] <- 2
		} else if (bkgd[[1]][1] == "green.") {
			data$pch[ic] <- 5
		} else if (bkgd[[1]][1] == "dist.") {
			data$pch[ic] <- 3
		}
	} else if (l_or_d == "2") {
		if (bkgd[[1]][1] == "blue") {
			data$pch[ic] <- 15
		} else if (bkgd[[1]][1] == "yellow.") {
			data$pch[ic] <- 16
		} else if (bkgd[[1]][1] == "red.") {
			data$pch[ic] <- 17
		} else if (bkgd[[1]][1] == "green.") {
			data$pch[ic] <- 18
		} else if (bkgd[[1]][1] == "dist.") {
			data$pch[ic] <- 4
		}
	}
}

mypalette()

pdf("../../../figures/rmc_rsd_from_simulation.pdf", height = 3, width = 9)
mypar(mfrow=c(1,3))
for(d in c("l", "m", "s")) {
	myplot(data[[paste0("rmc_", d)]],
		data[[paste0("rsd_", d)]],
		xlim=c(0,1),
		ylim=c(0,3.5),
		col=data$color,
		pch=data$pch)
	abline(0, 1)
	myaxis(1, at = seq(0, 1, 1/5))
	myaxis(2, at = seq(0, 3.5, 3.5/5))
	mytitle(toupper(d), xlab = "Ratio Mean Cone Exc.", ylab = "Ratio Sta. Dev. Cone Exc.")
}
dev.off()

data.from.exp <- read.csv('../../data/all_obs_data_together.csv')
obs.mean.data <- aggregate(. ~ sub_id + img + illum + rg + by + lighter_darker + bkgd, data.from.exp, mean)
pop.mean.data <- aggregate(. ~ img + illum + rg + by + lighter_darker + bkgd, obs.mean.data, function(x) { c(MN=mean(x), SD=sd(x)) })

pdf("../../../figures/rmc_rsd_comp.pdf", height = 3, width = 9)
mypar(mfrow=c(1,3))
plotRMCandRSD(pop.mean.data, obs.mean.data, 3, 3, "topleft")
dev.off()

# without the two weird darker images... they still both don't explain it...
obs.mean.data.no.weird <- obs.mean.data[13:nrow(obs.mean.data), ]
pop.mean.data.no.weird <- pop.mean.data[3:nrow(pop.mean.data), ]

pdf("../../../figures/rmc_rsd_comp_wout_darker_diff_dist.pdf", height = 3, width = 9)
mypar(mfrow=c(1,3))
plotRMCandRSD(pop.mean.data.no.weird, obs.mean.data.no.weird, 3, 3, "topleft")
dev.off()

