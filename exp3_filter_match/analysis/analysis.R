source("../../aux/plotFuncs.R")

data <- read.csv('../data/all_obs_data_together.csv')

obs.mean.data <- aggregate(. ~ sub_id + img + illum + body + highlow, data, mean)
pop.mean.data <- aggregate(. ~ img + illum + body + highlow, obs.mean.data, function(x) { c(MN=mean(x), SD=sd(x)) })

pop.mean.data$color <- NA

pop.mean.data[pop.mean.data$body == "red" & pop.mean.data$highlow == 0, ]$color <- 8
pop.mean.data[pop.mean.data$body == "green" & pop.mean.data$highlow == 0, ]$color <- 9
pop.mean.data[pop.mean.data$body == "blue" & pop.mean.data$highlow == 0, ]$color <- 10
pop.mean.data[pop.mean.data$body == "yellow" & pop.mean.data$highlow == 0, ]$color <- 11

pop.mean.data[pop.mean.data$body == "red" & pop.mean.data$highlow == 1, ]$color <- 4
pop.mean.data[pop.mean.data$body == "green" & pop.mean.data$highlow == 1, ]$color <- 5
pop.mean.data[pop.mean.data$body == "blue" & pop.mean.data$highlow == 1, ]$color <- 6
pop.mean.data[pop.mean.data$body == "yellow" & pop.mean.data$highlow == 1, ]$color <- 14

pop.mean.data$pch <- NA

pop.mean.data[pop.mean.data$illum == "blue", ]$pch <- 19
pop.mean.data[pop.mean.data$illum == "yellow", ]$pch <- 6

mypalette()

# plotting results for whole obj in terms of mean LAB
pdf("../../figures/exp3_mean_lab.pdf", height = 3, width = 9)
mypar(mfrow=c(1,3))
plotLABMean(pop.mean.data, obs.mean.data, "mean_", FALSE)
dev.off()

# plotting results for whole obj in terms of max LAB
pdf("../../figures/exp3_max_lab.pdf", height = 3, width = 9)
mypar(mfrow=c(1,3))
for (i in c("L", "A", "B")) {
	pop_x <- pop.mean.data[[paste0("wp", i)]][,1]
	pop_y <- pop.mean.data[[paste0("max_", tolower(i))]][,1]
	pop_sd <- pop.mean.data[[paste0("max_", tolower(i))]][,2]

	obs_x <- obs.mean.data[[paste0("wp", i)]]
	obs_y <- obs.mean.data[[paste0("max_", tolower(i))]]

	if (i == "L") {
		minlim <- 0
		maxlim <- 100
	} else {
		minlim <- -100
		maxlim <- 100
	}

	myplot(pop_x,
		pop_y,
		xlim=c(minlim,maxlim),
		ylim=c(minlim,maxlim),
		col = pop.mean.data$color,
		pch = pop.mean.data$pch)
	my_V_ebar(pop_x, pop_y, pop_sd, col=pop.mean.data$color)
	abline(0, 1)
	abline(rg <- lm(obs_y ~ obs_x), lty=2, col="darkslategray")

	rsq <- round(summary(rg)$adj.r.squared, digits=2)
	if (i == "L") {
		mytitle(paste0("Adjusted R Squared = ", rsq),
			xlab = paste0("White Point Image ", i, "*"),
			ylab = paste0("White Point Observer Match ", i, "*"))
	} else {
		mytitle(paste0("Adjusted R Squared = ", rsq),
			xlab = paste0("White Point Image ", tolower(i), "*"),
			ylab = paste0("White Point Observer Match ", tolower(i), "*"))
	}

	if (i == "L") {
		myaxis(1, at = seq(0, 100, 20))
		myaxis(2, at = seq(0, 100, 20))
		mylegend("topleft", c("Blue Illuminant", "Yellow Illuminant"), pch=c(19, 6))
	} else if (i == "A") {
		myaxis(1, at = seq(-50, 50, 50))
		myaxis(2, at = seq(-50, 50, 50))
		myaxis(1, at = -100, labels = "-100", col.axis = 5)
		myaxis(1, at = 100, labels = "100", col.axis = 4)
		myaxis(2, at = -100, labels = "-100", col.axis = 5)
		myaxis(2, at = 100, labels = "100", col.axis = 4)
	} else if (i == "B") {
		myaxis(1, at = seq(-50, 50, 50))
		myaxis(2, at = seq(-50, 50, 50))
		myaxis(1, at = -100, labels = "-100", col.axis = 6)
		myaxis(1, at = 100, labels = "100", col.axis = 7)
		myaxis(2, at = -100, labels = "-100", col.axis = 6)
		myaxis(2, at = 100, labels = "100", col.axis = 7)
	}
}
dev.off()

# plotting results for ratios of mean and sd. of cone excitations
pdf("../../figures/cer_and_sd_results.pdf", height = 3, width = 9)
mypar(mfrow=c(1,3))
plotRMCandRSD(pop.mean.data, obs.mean.data, 2, 2, "topleft")
dev.off()
