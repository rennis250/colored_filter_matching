data <- read.csv('test.csv')

obj.data <- data[data$mask == "obj_mask.png", ]
bkg.data <- data[data$mask == "bkg_mask.png", ]

obj.data$rmc_l <- obj.data$mean_cone_l/bkg.data$mean_cone_l
obj.data$rmc_m <- obj.data$mean_cone_m/bkg.data$mean_cone_m
obj.data$rmc_s <- obj.data$mean_cone_s/bkg.data$mean_cone_s

obj.data$rsd_l <- obj.data$sd_cone_l/bkg.data$sd_cone_l
obj.data$rsd_m <- obj.data$sd_cone_m/bkg.data$sd_cone_m
obj.data$rsd_s <- obj.data$sd_cone_s/bkg.data$sd_cone_s

data.from.exp <- read.csv('../../../exp3_filter_match/data/all_obs_data_together.csv')

obs.mean.data <- aggregate(. ~ sub_id + img + illum + body + highlow + mask, data.from.exp, mean)
pop.mean.data <- aggregate(. ~ img + illum + body + highlow + mask, obs.mean.data, mean)

obj.data.from.exp <- as.data.frame(pop.mean.data[pop.mean.data$mask == "obj_mask.png", ])
bkg.data.from.exp <- as.data.frame(pop.mean.data[pop.mean.data$mask == "bkg_mask.png", ])

obj.data.from.exp$color <- NA

obj.data.from.exp[obj.data.from.exp$body == "red" & obj.data.from.exp$highlow == 0, ]$color <- rgb(255, 0, 0, maxColorValue = 255)
obj.data.from.exp[obj.data.from.exp$body == "green" & obj.data.from.exp$highlow == 0, ]$color <- rgb(0, 255, 0, maxColorValue = 255)
obj.data.from.exp[obj.data.from.exp$body == "blue" & obj.data.from.exp$highlow == 0, ]$color <- rgb(0, 0, 255, maxColorValue = 255)
obj.data.from.exp[obj.data.from.exp$body == "yellow" & obj.data.from.exp$highlow == 0, ]$color <- rgb(255, 126, 0, maxColorValue = 255)

obj.data.from.exp[obj.data.from.exp$body == "red" & obj.data.from.exp$highlow == 1, ]$color <- rgb(126, 0, 0, maxColorValue = 255)
obj.data.from.exp[obj.data.from.exp$body == "green" & obj.data.from.exp$highlow == 1, ]$color <- rgb(0, 126, 0, maxColorValue = 255)
obj.data.from.exp[obj.data.from.exp$body == "blue" & obj.data.from.exp$highlow == 1, ]$color <- rgb(0, 0, 126, maxColorValue = 255)
obj.data.from.exp[obj.data.from.exp$body == "yellow" & obj.data.from.exp$highlow == 1, ]$color <- rgb(126, 63, 0, maxColorValue = 255)

obj.data.from.exp$pch <- NA

obj.data.from.exp[obj.data.from.exp$illum == "blue", ]$pch <- 19
obj.data.from.exp[obj.data.from.exp$illum == "yellow", ]$pch <- 6

obj.data.from.exp$rmc_l <- obj.data.from.exp$mean_cone_l/bkg.data.from.exp$mean_cone_l
obj.data.from.exp$rmc_m <- obj.data.from.exp$mean_cone_m/bkg.data.from.exp$mean_cone_m
obj.data.from.exp$rmc_s <- obj.data.from.exp$mean_cone_s/bkg.data.from.exp$mean_cone_s

obj.data.from.exp$rsd_l <- obj.data.from.exp$sd_cone_l/bkg.data.from.exp$sd_cone_l
obj.data.from.exp$rsd_m <- obj.data.from.exp$sd_cone_m/bkg.data.from.exp$sd_cone_m
obj.data.from.exp$rsd_s <- obj.data.from.exp$sd_cone_s/bkg.data.from.exp$sd_cone_s

pdf("../../../figures/rmc_rsd_comp.pdf", height = 5, width = 15)

par(mfrow=c(1,3), pty="s", cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
plot(obj.data$rmc_l, obj.data$rsd_l, xlim=c(0,1.3), ylim=c(0,1.3), ann = FALSE, cex = 1.5)
points(obj.data.from.exp$rmc_l, obj.data.from.exp$rsd_l, col=obj.data.from.exp$color, pch=obj.data.from.exp$pch, cex = 1.5)
abline(0,1, lty=2,)
title(ylab = "Ratio SD CE", main = "L cone")
rmcs <- c(obj.data$rmc_l, obj.data.from.exp$rmc_l)
rsds <- c(obj.data$rsd_l, obj.data.from.exp$rsd_l)
df <- data.frame(rmcs, rsds)
df <- df[is.finite(rowSums(df)),]
cr <- cor(df$rmcs, df$rsds)
text(1.2, 0.0, paste0(round(cr, digit=3)), cex=2.0)

plot(obj.data$rmc_m, obj.data$rsd_m, xlim=c(0,1.3), ylim=c(0,1.3), ann = FALSE, cex = 1.5)
points(obj.data.from.exp$rmc_m, obj.data.from.exp$rsd_m, col=obj.data.from.exp$color, pch=obj.data.from.exp$pch, cex = 1.5)
abline(0,1, lty=2,)
title(xlab = "Ratio Mean CE", main = "M cone")
rmcs <- c(obj.data$rmc_m, obj.data.from.exp$rmc_m)
rsds <- c(obj.data$rsd_m, obj.data.from.exp$rsd_m)
df <- data.frame(rmcs, rsds)
df <- df[is.finite(rowSums(df)),]
cr <- cor(df$rmcs, df$rsds)
text(1.2, 0.0, paste0(round(cr, digit=3)), cex=2.0)

plot(obj.data$rmc_s, obj.data$rsd_s, xlim=c(0,1.3), ylim=c(0,1.3), ann = FALSE, cex = 1.5)
points(obj.data.from.exp$rmc_s, obj.data.from.exp$rsd_s, col=obj.data.from.exp$color, pch=obj.data.from.exp$pch, cex = 1.5)
title(main = "S cone")
abline(0,1, lty=2,)
rmcs <- c(obj.data$rmc_s, obj.data.from.exp$rmc_s)
rsds <- c(obj.data$rsd_s, obj.data.from.exp$rsd_s)
df <- data.frame(rmcs, rsds)
df <- df[is.finite(rowSums(df)),]
cr <- cor(df$rmcs, df$rsds)
text(1.2, 0.0, paste0(round(cr, digit=3)), cex=2.0)

dev.off()

