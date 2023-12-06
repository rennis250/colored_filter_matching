source("../aux/plotFuncs.R")

pop.mean.data <- read.csv('../data/popstats_for_R.csv')
obs.mean.data <- read.csv('../data/obsstats_for_R.csv')

# wout_high mask is necessary here!
pop.mean.data <- pop.mean.data[pop.mean.data$exp_name == "patch_dye" & pop.mean.data$mask_name == "obj_wout_high_mask", ]
obs.mean.data <- obs.mean.data[obs.mean.data$exp_name == "patch_dye" & obs.mean.data$mask_name == "obj_wout_high_mask", ]

pop.mean.data$color <- NA

pop.mean.data[pop.mean.data$body_glaven  == "red" & pop.mean.data$hilo_glaven  == 0, ]$color <- 8
pop.mean.data[pop.mean.data$body_glaven  == "green" & pop.mean.data$hilo_glaven  == 0, ]$color <- 9
pop.mean.data[pop.mean.data$body_glaven  == "blue" & pop.mean.data$hilo_glaven  == 0, ]$color <- 10
pop.mean.data[pop.mean.data$body_glaven  == "yellow" & pop.mean.data$hilo_glaven  == 0, ]$color <- 11

pop.mean.data[pop.mean.data$body_glaven  == "red" & pop.mean.data$hilo_glaven  == 1, ]$color <- 4
pop.mean.data[pop.mean.data$body_glaven  == "green" & pop.mean.data$hilo_glaven  == 1, ]$color <- 5
pop.mean.data[pop.mean.data$body_glaven  == "blue" & pop.mean.data$hilo_glaven  == 1, ]$color <- 6
pop.mean.data[pop.mean.data$body_glaven  == "yellow" & pop.mean.data$hilo_glaven  == 1, ]$color <- 14

pop.mean.data$pch <- NA

pop.mean.data[pop.mean.data$illum_glaven  == "blue", ]$pch <- 19
pop.mean.data[pop.mean.data$illum_glaven  == "yellow", ]$pch <- 6

obs.mean.data$color <- NA

obs.mean.data[obs.mean.data$body_glaven  == "red" & obs.mean.data$hilo_glaven  == 0, ]$color <- 8
obs.mean.data[obs.mean.data$body_glaven  == "green" & obs.mean.data$hilo_glaven  == 0, ]$color <- 9
obs.mean.data[obs.mean.data$body_glaven  == "blue" & obs.mean.data$hilo_glaven  == 0, ]$color <- 10
obs.mean.data[obs.mean.data$body_glaven  == "yellow" & obs.mean.data$hilo_glaven  == 0, ]$color <- 11

obs.mean.data[obs.mean.data$body_glaven  == "red" & obs.mean.data$hilo_glaven  == 1, ]$color <- 4
obs.mean.data[obs.mean.data$body_glaven  == "green" & obs.mean.data$hilo_glaven  == 1, ]$color <- 5
obs.mean.data[obs.mean.data$body_glaven  == "blue" & obs.mean.data$hilo_glaven  == 1, ]$color <- 6
obs.mean.data[obs.mean.data$body_glaven  == "yellow" & obs.mean.data$hilo_glaven  == 1, ]$color <- 14

obs.mean.data$pch <- NA

obs.mean.data[obs.mean.data$illum_glaven  == "blue", ]$pch <- 19
obs.mean.data[obs.mean.data$illum_glaven  == "yellow", ]$pch <- 6

mypalette()

# plotting results for whole obj in terms of most saturated and most frequent color
pdf("../figures/exp3_most_sat_most_freq.pdf", height = 9, width = 9)
mypar(mfrow=c(3,3))
for (i in c("L", "a", "b")) {
  pop_x <- pop.mean.data[[paste0("mean_mean_msat_", i, "_glaven")]]
  pop_x <- pop_x[!is.nan(pop_x)]
  pop_y <- pop.mean.data[[paste0("mean_mean_gw", i, "_match")]]
  pop_y <- pop_y[!is.nan(pop_y)]
  pop_sd <- pop.mean.data[[paste0("std_mean_gw", i, "_match")]]
  
  obs_x <- obs.mean.data[[paste0("mean_msat_", i, "_glaven")]]
  obs_x <- obs_x[!is.nan(obs_x)]
  obs_y <- obs.mean.data[[paste0("mean_gw", i, "_match")]]
  obs_y <- obs_y[!is.nan(obs_y)]
  
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
            xlab = paste0("Most Sat. Image ", i, "*"),
            ylab = paste0("Patch Match ", i, "*"))
  } else {
    mytitle(paste0("Adjusted R Squared = ", rsq),
            xlab = paste0("Most Sat. Image ", tolower(i), "*"),
            ylab = paste0("Patch Match ", tolower(i), "*"))
  }
  
  if (i == "L") {
    myaxis(1, at = seq(0, 100, 20))
    myaxis(2, at = seq(0, 100, 20))
    mylegend("bottomright", c("Blue Illuminant", "White Illuminant"), pch=c(19, 6))
  } else if (i == "a") {
    myaxis(1, at = seq(-50, 50, 50))
    myaxis(2, at = seq(-50, 50, 50))
    myaxis(1, at = -100, labels = "-100", col.axis = 5)
    myaxis(1, at = 100, labels = "100", col.axis = 4)
    myaxis(2, at = -100, labels = "-100", col.axis = 5)
    myaxis(2, at = 100, labels = "100", col.axis = 4)
  } else if (i == "b") {
    myaxis(1, at = seq(-50, 50, 50))
    myaxis(2, at = seq(-50, 50, 50))
    myaxis(1, at = -100, labels = "-100", col.axis = 6)
    myaxis(1, at = 100, labels = "100", col.axis = 7)
    myaxis(2, at = -100, labels = "-100", col.axis = 6)
    myaxis(2, at = 100, labels = "100", col.axis = 7)
  }
}

for (i in c("L", "a", "b")) {
  pop_x <- pop.mean.data[[paste0("mean_mean_mfreq_", i, "_glaven")]]
  pop_x <- pop_x[!is.nan(pop_x)]
  pop_y <- pop.mean.data[[paste0("mean_mean_gw", i, "_match")]]
  pop_y <- pop_y[!is.nan(pop_y)]
  pop_sd <- pop.mean.data[[paste0("std_mean_gw", i, "_match")]]
  
  obs_x <- obs.mean.data[[paste0("mean_mfreq_", i, "_glaven")]]
  obs_x <- obs_x[!is.nan(obs_x)]
  obs_y <- obs.mean.data[[paste0("mean_gw", i, "_match")]]
  obs_y <- obs_y[!is.nan(obs_y)]
  
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
            xlab = paste0("Most Freq. Image ", i, "*"),
            ylab = paste0("Patch Match ", i, "*"))
  } else {
    mytitle(paste0("Adjusted R Squared = ", rsq),
            xlab = paste0("Most Freq. Image ", tolower(i), "*"),
            ylab = paste0("Patch Match ", tolower(i), "*"))
  }
  
  if (i == "L") {
    myaxis(1, at = seq(0, 100, 20))
    myaxis(2, at = seq(0, 100, 20))
    mylegend("bottomright", c("Blue Illuminant", "White Illuminant"), pch=c(19, 6))
  } else if (i == "a") {
    myaxis(1, at = seq(-50, 50, 50))
    myaxis(2, at = seq(-50, 50, 50))
    myaxis(1, at = -100, labels = "-100", col.axis = 5)
    myaxis(1, at = 100, labels = "100", col.axis = 4)
    myaxis(2, at = -100, labels = "-100", col.axis = 5)
    myaxis(2, at = 100, labels = "100", col.axis = 4)
  } else if (i == "b") {
    myaxis(1, at = seq(-50, 50, 50))
    myaxis(2, at = seq(-50, 50, 50))
    myaxis(1, at = -100, labels = "-100", col.axis = 6)
    myaxis(1, at = 100, labels = "100", col.axis = 7)
    myaxis(2, at = -100, labels = "-100", col.axis = 6)
    myaxis(2, at = 100, labels = "100", col.axis = 7)
  }
}

for (i in c("L", "a", "b")) {
  pop_x <- pop.mean.data[[paste0("mean_mean_wp_top5", i, "_glaven")]]
  pop_x <- pop_x[!is.nan(pop_x)]
  pop_y <- pop.mean.data[[paste0("mean_mean_gw", i, "_match")]]
  pop_y <- pop_y[!is.nan(pop_y)]
  pop_sd <- pop.mean.data[[paste0("std_mean_gw", i, "_match")]]
  
  obs_x <- obs.mean.data[[paste0("mean_wp_top5", i, "_glaven")]]
  obs_x <- obs_x[!is.nan(obs_x)]
  obs_y <- obs.mean.data[[paste0("mean_gw", i, "_match")]]
  obs_y <- obs_y[!is.nan(obs_y)]
  
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
            ylab = paste0("Patch Match ", i, "*"))
  } else {
    mytitle(paste0("Adjusted R Squared = ", rsq),
            xlab = paste0("White Point Image ", tolower(i), "*"),
            ylab = paste0("Patch Match ", tolower(i), "*"))
  }
  
  if (i == "L") {
    myaxis(1, at = seq(0, 100, 20))
    myaxis(2, at = seq(0, 100, 20))
    mylegend("bottomright", c("Blue Illuminant", "White Illuminant"), pch=c(19, 6))
  } else if (i == "a") {
    myaxis(1, at = seq(-50, 50, 50))
    myaxis(2, at = seq(-50, 50, 50))
    myaxis(1, at = -100, labels = "-100", col.axis = 5)
    myaxis(1, at = 100, labels = "100", col.axis = 4)
    myaxis(2, at = -100, labels = "-100", col.axis = 5)
    myaxis(2, at = 100, labels = "100", col.axis = 4)
  } else if (i == "b") {
    myaxis(1, at = seq(-50, 50, 50))
    myaxis(2, at = seq(-50, 50, 50))
    myaxis(1, at = -100, labels = "-100", col.axis = 6)
    myaxis(1, at = 100, labels = "100", col.axis = 7)
    myaxis(2, at = -100, labels = "-100", col.axis = 6)
    myaxis(2, at = 100, labels = "100", col.axis = 7)
  }
}
dev.off()

# and plot out the rmc for uniform patch

# plotting results for whole obj in terms of most saturated and most frequent color
pdf("../figures/exp3_rmc_patch.pdf", height = 3, width = 9)
mypar(mfrow=c(1,3))
for (i in c("l", "m", "s")) {
  pop_x <- pop.mean.data[[paste0("mean_mean_rmc_", i, "_glaven")]]
  pop_x <- pop_x[!is.nan(pop_x)]
  pop_y <- pop.mean.data[[paste0("mean_mean_rmc_", i, "_match")]]
  pop_y <- pop_y[!is.nan(pop_y)]
  pop_sd <- pop.mean.data[[paste0("std_mean_rmc_", i, "_match")]]
  
  obs_x <- obs.mean.data[[paste0("mean_rmc_", i, "_glaven")]]
  obs_x <- obs_x[!is.nan(obs_x)]
  obs_y <- obs.mean.data[[paste0("mean_rmc_", i, "_match")]]
  obs_y <- obs_y[!is.nan(obs_y)]
  
  minlim <- 0
  maxlim <- 2
  
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
            xlab = paste0("RMC Image ", toupper(i), "*"),
            ylab = paste0("RMC Patch Match ", toupper(i), "*"))
  } else {
    mytitle(paste0("Adjusted R Squared = ", rsq),
            xlab = paste0("RMC Image ", toupper(i), "*"),
            ylab = paste0("RMC Patch Match ", toupper(i), "*"))
  }
  
  myaxis(1, at = seq(0, 2, 2/5))
  myaxis(2, at = seq(0, 2, 2/5))
  if (i == "l") {
    mylegend("bottomright", c("Blue Illuminant", "White Illuminant"), pch=c(19, 6))
  }
}
dev.off()


# make a plot showing the color matches in a chromaticity diagram
# pdf("../figures/chrom_diagram_exps_1_2_3.pdf", height = 6.5, width = 6)
# exps <- c("Patch proximal match", "Patch dye match", "Flat filter - mean color", "Flat filter - White Point")
# mypar(mfrow = c(2, 2), pty="s")
# for (i in 1:4) {
#   if (i < 3) {
#     img_a <- pop_mean_data[[i]][["gwa"]][, 1]
#     img_b <- pop_mean_data[[i]][["gwb"]][, 1]
#     
#     pop_a <- pop_mean_data[[i]][["a"]][, 1]
#     pop_b <- pop_mean_data[[i]][["b"]][, 1]
#     
#     pop_a_sd <- pop_mean_data[[i]][["a"]][, 2]
#     pop_b_sd <- pop_mean_data[[i]][["b"]][, 2]
#   } else if (i == 3) {
#     img_a <- pop_mean_data[[i]][["gwa"]][, 1]
#     img_b <- pop_mean_data[[i]][["gwb"]][, 1]
#     
#     pop_a <- pop_mean_data[[i]][["mean_a"]][, 1]
#     pop_b <- pop_mean_data[[i]][["mean_b"]][, 1]
#     
#     pop_a_sd <- pop_mean_data[[i]][["mean_a"]][, 2]
#     pop_b_sd <- pop_mean_data[[i]][["mean_b"]][, 2]
#   } else {
#     img_a <- pop_mean_data[[3]][["gwa"]][, 1]
#     img_b <- pop_mean_data[[3]][["gwb"]][, 1]
#     
#     pop_a <- pop_mean_data[[3]][["max_a"]][, 1]
#     pop_b <- pop_mean_data[[3]][["max_b"]][, 1]
#     
#     pop_a_sd <- pop_mean_data[[3]][["max_a"]][, 2]
#     pop_b_sd <- pop_mean_data[[3]][["max_b"]][, 2]
#   }
#   
#   myplot(0, 0, xlim = c(-100, 100), ylim = c(-100, 100), col = "white", cex = 1)
#   mypoints(img_a,
#            img_b,
#            xlim = c(-100, 100),
#            ylim = c(-100, 100),
#            col = hex(as(LAB(50, img_a * 0.7, img_b * 0.7), "RGB")),
#            pch = 1,
#            cex = 1
#   )
#   mypoints(pop_a,
#            pop_b,
#            col = "black",
#            bg = hex(as(LAB(50, pop_a * 0.7, pop_b * 0.7), "RGB")),
#            pch = 21,
#            cex = 1
#   )
#   draw.ellipse(pop_a, pop_b, a = pop_a_sd, b = pop_b_sd)
#   arrows(img_a, img_b, pop_a, pop_b, length = 0.05, lwd = 1)
#   
#   mytitle(xlab = "a*", ylab = "b*", main = exps[[i]])
#   
#   myaxis(1, at = seq(-50, 50, 50))
#   myaxis(1, at = 100, labels = "100", col.axis = 4)
#   myaxis(1, at = -100, labels = "-100", col.axis = 5)
#   
#   myaxis(2, at = seq(-50, 50, 50))
#   myaxis(2, at = 100, labels = "100", col.axis = 7)
#   myaxis(2, at = -100, labels = "-100", col.axis = 6)
# }
# dev.off()

# let's see what the robust RSD looks like for rmc_rsd_comp vs. normal rsd
pop.mean.data <- read.csv('../data/popstats_for_R.csv')
obs.mean.data <- read.csv('../data/obsstats_for_R.csv')

pop.mean.data <- pop.mean.data[pop.mean.data$exp_name == "rmc_rsd_comp" & pop.mean.data$mask_name == "obj_mask", ]
obs.mean.data <- obs.mean.data[obs.mean.data$exp_name == "rmc_rsd_comp" & obs.mean.data$mask_name == "obj_mask", ]

mypalette()

maxx <- 1.5
maxy <- 1.5
pdf("../figures/exp3_robust_rsd.pdf", height = 3, width = 9)
mypar(mfrow=c(1,3))
for (i in c("l", "m", "s")) {
  pop_taug_x <- pop.mean.data[[paste0("mean_mean_tau_general_", i, "_glaven")]]
  pop_taug_y <- pop.mean.data[[paste0("mean_mean_tau_general_", i, "_filter")]]
  pop_taus_x <- pop.mean.data[[paste0("mean_mean_tau_simpler_", i, "_glaven")]]
  pop_taus_y <- pop.mean.data[[paste0("mean_mean_tau_simpler_", i, "_filter")]]
  pop_taug_sd <- pop.mean.data[[paste0("std_mean_tau_general_", i, "_filter")]]
  pop_taus_sd <- pop.mean.data[[paste0("std_mean_tau_simpler_", i, "_filter")]]
  
  obs_taug_x <- obs.mean.data[[paste0("mean_tau_general_", i, "_glaven")]]
  obs_taug_y <- obs.mean.data[[paste0("mean_tau_general_", i, "_filter")]]
  obs_taus_x <- obs.mean.data[[paste0("mean_tau_simpler_", i, "_glaven")]]
  obs_taus_y <- obs.mean.data[[paste0("mean_tau_simpler_", i, "_filter")]]
  
  pop_rmc_x <- pop.mean.data[[paste0("mean_mean_rmc_", i, "_glaven")]]
  pop_rmc_y <- pop.mean.data[[paste0("mean_mean_rmc_", i, "_match")]]
  pop_rsd_x <- pop.mean.data[[paste0("mean_mean_rsd_", i, "_glaven")]]
  pop_rsd_y <- pop.mean.data[[paste0("mean_mean_rsd_", i, "_filter")]]
  pop_rmc_sd <- pop.mean.data[[paste0("std_mean_rmc_", i, "_match")]]
  pop_rsd_sd <- pop.mean.data[[paste0("std_mean_rsd_", i, "_filter")]]
  
  obs_rmc_x <- obs.mean.data[[paste0("mean_rmc_", i, "_glaven")]]
  obs_rmc_y <- obs.mean.data[[paste0("mean_rmc_", i, "_match")]]
  obs_rsd_x <- obs.mean.data[[paste0("mean_rsd_", i, "_glaven")]]
  obs_rsd_y <- obs.mean.data[[paste0("mean_rsd_", i, "_filter")]]
  
  myplot(pop_taug_x,
         pop_taug_y,
         xlim = c(0, maxx),
         ylim = c(0, maxy),
         col = 1,
         pch = 2,
  )
  mypoints(pop_taus_x,
           pop_taus_y,
           col = 6,
           pch = 6
  )
  mypoints(pop_rsd_x,
           pop_rsd_y,
           col = 7,
           pch = 8
  )
  mypoints(pop_rmc_x,
           pop_rmc_y,
           col = 4,
           pch = 8
  )
  my_V_ebar(pop_taug_x, pop_taug_y, pop_taug_sd, col = 1)
  my_V_ebar(pop_taus_x, pop_taus_y, pop_taus_sd, col = 6)
  my_V_ebar(pop_rsd_x, pop_rsd_y, pop_rsd_sd, col = 7)
  my_V_ebar(pop_rmc_x, pop_rmc_y, pop_rmc_sd, col = 4)
  abline(0, 1)
  abline(r_taug <- lm(obs_taug_y ~ obs_taug_x), lty = 2, col = 1)
  abline(r_taus <- lm(obs_taus_y ~ obs_taus_x), lty = 2, col = 6)
  abline(r_rsd <- lm(obs_rsd_y ~ obs_rsd_x), lty = 2, col = 7)
  abline(r_rmc <- lm(obs_rmc_y ~ obs_rmc_x), lty = 2, col = 4)
  myaxis(1, at = seq(0, maxx, maxx / 5))
  myaxis(2, at = seq(0, maxy, maxy / 5))
  
  taug_rsq <- round(summary(r_taug)$adj.r.squared, digits = 2)
  taus_rsq <- round(summary(r_taus)$adj.r.squared, digits = 2)
  rsd_rsq <- round(summary(r_rsd)$adj.r.squared, digits = 2)
  rmc_rsq <- round(summary(r_rmc)$adj.r.squared, digits = 2)
  mymtext("Adjusted R Squared = ", line = 0.6)
  mytitle(substitute(list(x, phantom(y)), list(x = taug_rsq, y = taus_rsq)),
          xlab = paste0("Robust Ratio ", toupper(i), " - Image"),
          ylab = paste0("Robust Ratio ", toupper(i), " - Filter Settings"),
          col.main = 1
  )
  mytitle(substitute(list(phantom(x), y), list(x = taug_rsq, y = taus_rsq)), col.main = 6)
  
  if (i == "l") {
    mylegend("topleft",
             c("Tau General", "Tau Simpler", "Std. Dev. Cone Exc.", "Mean Cone Exc."),
             pch = c(2, 6, 8, 8),
             col = c(1, 6, 7, 4)
    )
  }
}
dev.off()