source("../aux/plotFuncs.R")

data <- list(read.csv('../data/exp1_obs_data.csv'),
			read.csv('../data/exp2_obs_data.csv'),
			read.csv('../data/exp3_obs_data.csv'),
			read.csv('../data/exp11a_obs_data.csv'),
			read.csv('../data/exp11b_obs_data.csv'))

pop.mean.data <- list()
obs.mean.data <- list()

for (cc in 1:length(data)) {
	if (cc != 2) {
		obs.data <- aggregate(. ~ sub_id + img + illum + body + highlow, data[[cc]], mean)
		tmp.data <- aggregate(. ~ img + illum + body + highlow, obs.data, function(x) { c(MN=mean(x), SD=sd(x)) })
	} else {
		obs.data <- aggregate(. ~ sub_id + img + illum + body + highlow, data[[cc]], mean)
		tmp.data <- aggregate(. ~ img + illum + body + highlow, obs.data, function(x) { c(MN=mean(x), SD=sd(x)) })
	}

	tmp.data$color <- NA

	tmp.data[tmp.data$body == "red" & tmp.data$highlow == 0, ]$color <- 8
	tmp.data[tmp.data$body == "green" & tmp.data$highlow == 0, ]$color <- 9
	tmp.data[tmp.data$body == "blue" & tmp.data$highlow == 0, ]$color <- 10
	tmp.data[tmp.data$body == "yellow" & tmp.data$highlow == 0, ]$color <- 11

	tmp.data[tmp.data$body == "red" & tmp.data$highlow == 1, ]$color <- 4
	tmp.data[tmp.data$body == "green" & tmp.data$highlow == 1, ]$color <- 5
	tmp.data[tmp.data$body == "blue" & tmp.data$highlow == 1, ]$color <- 6
	tmp.data[tmp.data$body == "yellow" & tmp.data$highlow == 1, ]$color <- 14

	tmp.data$pch <- NA

	if (cc == 4 || cc == 5) {
		tmp.data$pch <- 19
	} else {
		tmp.data[tmp.data$illum == "blue", ]$pch <- 19
		tmp.data[tmp.data$illum == "yellow", ]$pch <- 6
	}

	pop.mean.data[[cc]] <- as.data.frame(tmp.data)
	obs.mean.data[[cc]] <- as.data.frame(obs.data)
}

rg <- list()

mypalette()

# plotting results for whole obj in terms of mean LAB
pdf("../figures/all_exps_mean_lab.pdf", height = 3.5, width = 11)
mypar(mfrow=c(1,3))

cols <- c(1, 6, 4, 5)
for (d in c("l", "a", "b")) {
	for (i in 1:4) { # doesn't actually access data from set 4. see hardcoded below
		if (i < 3) {
			pop_x <- pop.mean.data[[i]][[paste0("gw", d)]][,1]
			pop_y <- pop.mean.data[[i]][[d]][,1]
			pop_sd <- pop.mean.data[[i]][[d]][,2]

			obs_x <- obs.mean.data[[i]][[paste0("gw", d)]]
			obs_y <- obs.mean.data[[i]][[d]]
			plotch <- pop.mean.data[[i]]$pch
		} else if (i == 3) {
			pop_x <- pop.mean.data[[i]][[paste0("gw", d)]][,1]
			pop_y <- pop.mean.data[[i]][[paste0("mean_", d)]][,1]
			pop_sd <- pop.mean.data[[i]][[paste0("mean_", d)]][,2]

			obs_x <- obs.mean.data[[i]][[paste0("gw", d)]]
			obs_y <- obs.mean.data[[i]][[paste0("mean_", d)]]
			plotch <- pop.mean.data[[i]]$pch
		} else {
			pop_x <- pop.mean.data[[3]][[paste0("gw", d)]][,1]
			pop_y <- pop.mean.data[[3]][[paste0("max_", d)]][,1]
			pop_sd <- pop.mean.data[[3]][[paste0("max_", d)]][,2]

			obs_x <- obs.mean.data[[3]][[paste0("gw", d)]]
			obs_y <- obs.mean.data[[3]][[paste0("max_", d)]]
			plotch <- pop.mean.data[[3]]$pch
		}

		if (d == "l") {
			minlim <- 0
			maxlim <- 100
		} else {
			minlim <- -100
			maxlim <- 100
		}

		if(i == 1) {
			myplot(pop_x,
				pop_y,
				xlim=c(minlim,maxlim),
				ylim=c(minlim,maxlim),
				col = cols[i],
				pch = plotch)
		} else {
			mypoints(pop_x,
				pop_y,
				col = cols[i],
				pch = plotch)
		}
		abline(rg[[i]] <- lm(obs_y ~ obs_x), lty=2, col=cols[i])
	}
	abline(0, 1)
	if (d == "l") {
		mytitle(xlab = paste0("Mean Image ", toupper(d), "*"), ylab = paste0("Observer Match ", toupper(d), "*"))
	} else {
		mytitle(xlab = paste0("Mean Image ", d, "*"), ylab = paste0("Observer Match ", d, "*"))
	}

	if (d == "l") {
		myaxis(1, at = seq(0, 100, 20))
		myaxis(2, at = seq(0, 100, 20))
		mylegend("bottomright",
			c("Blue Illuminant", "Yellow Illuminant", "Mean Proximal Match", "Mean Dye Match", "Mean Filter Match", "White Point Filter Match"),
			pch=c(19, 6, NA, NA, NA, NA),
			lty=c(NA, NA, 2, 2, 2, 2),
			col = c(1, 1, 1, 6, 4, 5))
	} else if (d == "a") {
		myaxis(1, at = seq(-50, 50, 50))
		myaxis(2, at = seq(-50, 50, 50))
		myaxis(1, at = -100, labels = "-100", col.axis = 5)
		myaxis(1, at = 100, labels = "100", col.axis = 4)
		myaxis(2, at = -100, labels = "-100", col.axis = 5)
		myaxis(2, at = 100, labels = "100", col.axis = 4)
	} else if (d == "b") {
		myaxis(1, at = seq(-50, 50, 50))
		myaxis(2, at = seq(-50, 50, 50))
		myaxis(1, at = -100, labels = "-100", col.axis = 7)
		myaxis(1, at = 100, labels = "100", col.axis = 6)
		myaxis(2, at = -100, labels = "-100", col.axis = 7)
		myaxis(2, at = 100, labels = "100", col.axis = 6)
	}
}

### bar charts showing R^2

# mybarplot(c(summary(prox_r[[1]])$adj.r.squared,
# 	summary(dye_r[[1]])$adj.r.squared,
# 	summary(mean_filt_r[[1]])$adj.r.squared,
# 	summary(max_filt_r[[1]])$adj.r.squared),
# 	col=c(1, 6, 4, 5),
# 	ylim=c(0,1))
# mytitle(ylab="Adjusted R Squared")
#
# mybarplot(c(summary(prox_r[[2]])$adj.r.squared,
# 	summary(dye_r[[2]])$adj.r.squared,
# 	summary(mean_filt_r[[2]])$adj.r.squared,
# 	summary(max_filt_r[[2]])$adj.r.squared),
# 	col=c(1, 6, 4, 5),
# 	names.arg=c("Prox", "Dye", "Mean Filt", "Max Filt"), ylim=c(0,1))
# mytitle(xlab="Experimental Condition")
#
# mybarplot(c(summary(prox_r[[3]])$adj.r.squared,
# 	summary(dye_r[[3]])$adj.r.squared,
# 	summary(mean_filt_r[[3]])$adj.r.squared,
# 	summary(max_filt_r[[3]])$adj.r.squared),
# 	col=c(1, 6, 4, 5),
# 	ylim=c(0,1))
dev.off()

# plotting results for white bkg and illum in terms of mean LAB for uniform patch
pdf("../figures/white_wall_results_exp11a.pdf", height = 3, width = 9)
pop.mean.data[[4]]$color <- 1
pop.mean.data[[4]]$pch <- 16
obs.mean.data[[4]]$color <- 1
obs.mean.data[[4]]$pch <- 16
mypar(mfrow=c(1,3))
plotLABMean(pop.mean.data[[4]], obs.mean.data[[4]], "", TRUE)
dev.off()

# plotting results for white bkg and illum in terms of mean LAB for filter
pdf("../figures/white_wall_results_exp11b.pdf", height = 3, width = 9)
pop.mean.data[[5]]$color <- 1
pop.mean.data[[5]]$pch <- 16
obs.mean.data[[5]]$color <- 1
obs.mean.data[[5]]$pch <- 16
mypar(mfrow=c(1,3))
plotLABMean(pop.mean.data[[5]], obs.mean.data[[5]], "mean_", TRUE)
dev.off()