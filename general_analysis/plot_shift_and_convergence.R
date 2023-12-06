library(shape)
library(EBImage)

data <- read.csv('../exp3_filter_match/data/all_obs_data_together.csv')
obs.mean.data <- aggregate(. ~ sub_id + img + illum + body + highlow + mask, data, mean)
tmp.data <- aggregate(. ~ img + illum + body + highlow + mask, obs.mean.data, mean)
obj.data <- as.data.frame(tmp.data[tmp.data$mask == "obj_mask.png", ])

lab_vals <- read.csv('../images/patterned/mitsuba_caustics_rgb_trans_lightest_illum_blue_refl_munsell_blue_lightest_lab_coords.csv')
img <- readImage('../images/empty_scene_blue_illum.png')
img_w_filter <- readImage('../images/patterned/mitsuba_caustics_rgb_trans_lightest_illum_blue_refl_munsell_blue_lightest.png')

lab_filter <- lab_vals[ , 1:3]
lab_wout_filter <- lab_vals[ , 4:6]
lab_dark_filter <- lab_vals[ , 7:9]

random_pts <- sample(1:nrow(lab_filter), 1000)

pdf("../figures/dzmura_3d.pdf", height = 7, width = 25)
par(mar=c(2,5,1,1), mfrow=c(1,4), pty="s")

display(img, method="raster", margin=12)
display(img_w_filter, method="raster", margin=12)

plot(NA, NA, xlim = c(-120,120), ylim = c(-120,120), ylab = "b*", xlab = "a*",
	cex.lab=2, cex.axis=2)
Arrows(lab_wout_filter[random_pts,2], lab_wout_filter[random_pts,3], lab_filter[random_pts,2], lab_filter[random_pts,3], arr.col = "black", lcol = "grey")
Arrows(lab_wout_filter[random_pts,2], lab_wout_filter[random_pts,3], lab_dark_filter[random_pts,2], lab_dark_filter[random_pts,3], arr.col = "blue", lcol = "red")

points(obj.data$mean_a[[1]], obj.data$mean_b[[1]], cex=3, col="black", pch=19)
points(obj.data$mean_a[[9]], obj.data$mean_b[[9]], cex=3, col="black", pch=18)

plot(NA, NA, xlim = c(-120,120), ylim = c(0,100), ylab = "L*", xlab = "a*",
	cex.lab=2, cex.axis=2)
Arrows(lab_wout_filter[random_pts,2], lab_wout_filter[random_pts,1], lab_filter[random_pts,2], lab_filter[random_pts,1], arr.col = "black", lcol = "grey")
Arrows(lab_wout_filter[random_pts,2], lab_wout_filter[random_pts,1], lab_dark_filter[random_pts,2], lab_dark_filter[random_pts,1], arr.col = "blue", lcol = "red")

points(obj.data$mean_a[[1]], obj.data$mean_l[[1]], cex=3, col="black", pch=19)
points(obj.data$mean_a[[9]], obj.data$mean_l[[9]], cex=3, col="black", pch=18)

axis(1, at = -100, labels = "-100", col.axis = rgb(0, 255, 0, maxColorValue = 255), cex.axis=2)
axis(1, at = 100, labels = "100", col.axis = rgb(255, 0, 0, maxColorValue = 255), cex.axis=2)
axis(2, at = -100, labels = "-100", col.axis = rgb(0, 0, 255, maxColorValue = 255), cex.axis=2)
axis(2, at = 100, labels = "100", col.axis = rgb(255, 126, 0, maxColorValue = 255), cex.axis=2)

dev.off()
