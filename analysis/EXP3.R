source("../aux/plotFuncs.R")

# data <- read.csv('../data/exp3_obs_data.csv')
data <- read.csv("../data/data_table.csv")
data <- data[data$exp_name == "filter" & data$mask_name == "obj_mask", ]

obs.mean.data <- aggregate(data, list(body = data$body_glaven, illum = data$illum_glaven, obs = data$obs_name, hilo = data$hilo_glaven), mean)
obs.mean.data <- subset(obs.mean.data, select = -c(obs, obs_name, body_glaven, illum_glaven, hilo_glaven, exp_name, mask_name))
pop.mean.data <- aggregate(obs.mean.data, list(body_glaven = obs.mean.data$body, illum_glaven = obs.mean.data$illum, hilo_glaven = obs.mean.data$hilo), function(x) {
  c(MN = mean(x), SD = sd(x))
})
pop.mean.data <- subset(pop.mean.data, select = -c(body, illum, hilo))

# pop.mean.data$color <- NA

# pop.mean.data[pop.mean.data$body == "red" & pop.mean.data$highlow == 0, ]$color <- 8
# pop.mean.data[pop.mean.data$body == "green" & pop.mean.data$highlow == 0, ]$color <- 9
# pop.mean.data[pop.mean.data$body == "blue" & pop.mean.data$highlow == 0, ]$color <- 10
# pop.mean.data[pop.mean.data$body == "yellow" & pop.mean.data$highlow == 0, ]$color <- 11

# pop.mean.data[pop.mean.data$body == "red" & pop.mean.data$highlow == 1, ]$color <- 4
# pop.mean.data[pop.mean.data$body == "green" & pop.mean.data$highlow == 1, ]$color <- 5
# pop.mean.data[pop.mean.data$body == "blue" & pop.mean.data$highlow == 1, ]$color <- 6
# pop.mean.data[pop.mean.data$body == "yellow" & pop.mean.data$highlow == 1, ]$color <- 14

# pop.mean.data$pch <- NA

# pop.mean.data[pop.mean.data$illum == "blue", ]$pch <- 19
# pop.mean.data[pop.mean.data$illum == "yellow", ]$pch <- 6

mypalette()

# # plotting results for whole obj in terms of mean LAB
# pdf("../figures/exp3_mean_lab.pdf", height = 3, width = 9)
# mypar(mfrow=c(1,3))
# for (i in c("l", "a", "b")) {
# 	pop_x <- pop.mean.data[[paste0("gw", i)]][,1]
# 	pop_y <- pop.mean.data[[paste0("mean_", i)]][,1]
# 	pop_sd <- pop.mean.data[[paste0("mean_", i)]][,2]
# 
# 	obs_x <- obs.mean.data[[paste0("gw", i)]]
# 	obs_y <- obs.mean.data[[paste0("mean_", i)]]
# 
# 	if (i == "l") {
# 		minlim <- 10
# 		maxlim <- 70
# 		i <- toupper(i)
# 	} else if (i == "a") {
# 		minlim <- -35
# 		maxlim <- 65
# 	} else {
# 		minlim <- -25
# 		maxlim <- 45
# 	}
# 
# 	myplot(pop_x,
# 		pop_y,
# 		xlim=c(minlim,maxlim),
# 		ylim=c(minlim,maxlim),
# 		col = pop.mean.data$color,
# 		pch = pop.mean.data$pch)
# 	my_V_ebar(pop_x, pop_y, pop_sd, col=pop.mean.data$color)
# 	abline(0, 1)
# 	abline(rg <- lm(obs_y ~ obs_x), lty=2, col="darkslategray")
# 
# 	rsq <- round(summary(rg)$adj.r.squared, digits=2)
# 	mytitle(paste0("Adjusted R Squared = ", rsq),
# 		xlab = paste0("Mean Image ", i, "*"),
# 		ylab = paste0("Mean Observer Match ", i, "*"))
# 
# 	if (i == "L") {
# 		myaxis(1, at = seq(minlim, maxlim, (maxlim - minlim)/4))
# 		myaxis(2, at = seq(minlim, maxlim, (maxlim - minlim)/4))
# 		mylegend("bottomright", c("Blue Illuminant", "White Illuminant"), pch=c(19, 6))
# 	} else if (i == "a") {
# 		myaxis(1, at = seq(minlim, maxlim, (maxlim - minlim)/4))
# 		myaxis(2, at = seq(minlim, maxlim, (maxlim - minlim)/4))
# 		myaxis(1, at = -100, labels = "-100", col.axis = 5)
# 		myaxis(1, at = 100, labels = "100", col.axis = 4)
# 		myaxis(2, at = -100, labels = "-100", col.axis = 5)
# 		myaxis(2, at = 100, labels = "100", col.axis = 4)
# 	} else if (i == "b") {
# 		myaxis(1, at = seq(minlim, maxlim, (maxlim - minlim)/5))
# 		myaxis(2, at = seq(minlim, maxlim, (maxlim - minlim)/5))
# 		myaxis(1, at = -100, labels = "-100", col.axis = 6)
# 		myaxis(1, at = 100, labels = "100", col.axis = 7)
# 		myaxis(2, at = -100, labels = "-100", col.axis = 6)
# 		myaxis(2, at = 100, labels = "100", col.axis = 7)
# 	}
# }
# dev.off()
# 
# # plotting results for whole obj in terms of max LAB
# pdf("../figures/exp3_max_lab.pdf", height = 3, width = 9)
# mypar(mfrow=c(1,3))
# for (i in c("L", "A", "B")) {
# 	pop_x <- pop.mean.data[[paste0("wp", i)]][,1]
# 	pop_y <- pop.mean.data[[paste0("max_", tolower(i))]][,1]
# 	pop_sd <- pop.mean.data[[paste0("max_", tolower(i))]][,2]
# 
# 	obs_x <- obs.mean.data[[paste0("wp", i)]]
# 	obs_y <- obs.mean.data[[paste0("max_", tolower(i))]]
# 
# 	if (i == "L") {
# 		minlim <- 0
# 		maxlim <- 100
# 	} else {
# 		minlim <- -100
# 		maxlim <- 100
# 	}
# 
# 	myplot(pop_x,
# 		pop_y,
# 		xlim=c(minlim,maxlim),
# 		ylim=c(minlim,maxlim),
# 		col = pop.mean.data$color,
# 		pch = pop.mean.data$pch)
# 	my_V_ebar(pop_x, pop_y, pop_sd, col=pop.mean.data$color)
# 	abline(0, 1)
# 	abline(rg <- lm(obs_y ~ obs_x), lty=2, col="darkslategray")
# 
# 	rsq <- round(summary(rg)$adj.r.squared, digits=2)
# 	if (i == "L") {
# 		mytitle(paste0("Adjusted R Squared = ", rsq),
# 			xlab = paste0("White Point Image ", i, "*"),
# 			ylab = paste0("White Point Observer Match ", i, "*"))
# 	} else {
# 		mytitle(paste0("Adjusted R Squared = ", rsq),
# 			xlab = paste0("White Point Image ", tolower(i), "*"),
# 			ylab = paste0("White Point Observer Match ", tolower(i), "*"))
# 	}
# 
# 	if (i == "L") {
# 		myaxis(1, at = seq(0, 100, 20))
# 		myaxis(2, at = seq(0, 100, 20))
# 		mylegend("topleft", c("Blue Illuminant", "White Illuminant"), pch=c(19, 6))
# 	} else if (i == "A") {
# 		myaxis(1, at = seq(-50, 50, 50))
# 		myaxis(2, at = seq(-50, 50, 50))
# 		myaxis(1, at = -100, labels = "-100", col.axis = 5)
# 		myaxis(1, at = 100, labels = "100", col.axis = 4)
# 		myaxis(2, at = -100, labels = "-100", col.axis = 5)
# 		myaxis(2, at = 100, labels = "100", col.axis = 4)
# 	} else if (i == "B") {
# 		myaxis(1, at = seq(-50, 50, 50))
# 		myaxis(2, at = seq(-50, 50, 50))
# 		myaxis(1, at = -100, labels = "-100", col.axis = 6)
# 		myaxis(1, at = 100, labels = "100", col.axis = 7)
# 		myaxis(2, at = -100, labels = "-100", col.axis = 6)
# 		myaxis(2, at = 100, labels = "100", col.axis = 7)
# 	}
# }
# dev.off()

# plotting results for ratios of mean and sd. of cone excitations
# including the results for the white filter matching element
# data_white_walls <- read.csv("../data/exp11b_obs_data.csv")
data <- read.csv("../data/data_table.csv")
data_white_walls <- data[data$exp_name == "white_walls_filter" & data$mask_name == "obj_mask", ]

# apparently obs 3 and 9 are bad observers for the robust ratio model
data_white_walls <- data_white_walls[data_white_walls$obs_name != "03" & data_white_walls$obs_name != "09", ]

obs_mean_data_white_walls <- aggregate(data_white_walls, list(body = data_white_walls$body_glaven, illum = data_white_walls$illum_glaven, obs = data_white_walls$obs_name, hilo = data_white_walls$hilo_glaven), mean)
obs_mean_data_white_walls <- subset(obs_mean_data_white_walls, select = -c(obs, obs_name, body_glaven, illum_glaven, hilo_glaven, exp_name, mask_name))
pop_mean_data_white_walls <- aggregate(obs_mean_data_white_walls, list(body_glaven = obs_mean_data_white_walls$body, illum_glaven = obs_mean_data_white_walls$illum, hilo_glaven = obs_mean_data_white_walls$hilo), function(x) {
  c(MN = mean(x), SD = sd(x))
})
pop_mean_data_white_walls <- subset(pop_mean_data_white_walls, select = -c(body, illum, hilo))
pop_mean_data_white_walls$color <- 1
pop_mean_data_white_walls$pch <- 16

pdf("../figures/cer_and_sd_results.pdf", height = 6, width = 9)
mypar(mfrow=c(2,3))
plotRMCandRSD(pop.mean.data, obs.mean.data, 2, 2, "topleft")
plotRMCandRSD(pop_mean_data_white_walls, obs_mean_data_white_walls, 2, 2, "topleft")
dev.off()

# plot robust ratio model results

maxx <- 3
maxy <- 3
adjtext <- "Adjusted R Squared =     "
lineval <- 0.5
legendloc <- "bottomright"
pdf("../figures/cer_and_robust_sd_results.pdf", height = 6, width = 9)
mypar(mfrow = c(2, 3))

# first for exp 3
pop.data <- pop.mean.data
obs.data <- obs.mean.data
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
             c("Mean Cone Exc.", "Robust Ratio"),
             pch = c(2, 6),
             col = c(1, 6)
    )
  }
}

# then for exp 11b
pop.data <- pop_mean_data_white_walls
obs.data <- obs_mean_data_white_walls
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
             c("Mean Cone Exc.", "Robust Ratio"),
             pch = c(2, 6),
             col = c(1, 6)
    )
  }
}
dev.off()