library(colorspace)

source("../aux/plotFuncs.R")

# data <- read.csv("../data/exp11b_obs_data.csv")
data <- read.csv("../data/data_table.csv")
data <- data[data$exp_name == "white_walls_filter" & data$mask_name == "obj_mask", ]

# apparently obs 3 and 9 are bad observers for the robust ratio model
data <- data[data$obs_name != "03" & data$obs_name != "09", ]

obs_mean_data <- aggregate(data, list(body = data$body_glaven, illum = data$illum_glaven, obs = data$obs_name, hilo = data$hilo_glaven), mean)
obs_mean_data <- subset(obs_mean_data, select = -c(obs, obs_name, body_glaven, illum_glaven, hilo_glaven, exp_name, mask_name))
pop_mean_data <- aggregate(obs_mean_data, list(body_glaven = obs_mean_data$body, illum_glaven = obs_mean_data$illum, hilo_glaven = obs_mean_data$hilo), function(x) {
    c(MN = mean(x), SD = sd(x))
})
pop_mean_data <- subset(pop_mean_data, select = -c(body, illum, hilo))

pop_mean_data$color <- 1

pop_mean_data$pch <- 16

mypalette()

# loop over different background types for the matching stimulus
# for (bg in c("black", "white", "gray")) {
# myplotting results for whole obj in terms of mean LAB
# pdf(paste0("../figures/exp11b_mean_lab_", bg,".pdf"), height = 3, width = 9)
# pdf(paste0("../figures/exp11b_mean_lab_gray.pdf"), height = 3, width = 9)
# mypar(mfrow = c(1, 3))
# plotLABMean(pop_mean_data, obs_mean_data, "mean_", TRUE)
# dev.off()
# # }
# 
# img_a <- pop_mean_data[["gwa"]][, 1]
# img_b <- pop_mean_data[["gwb"]][, 1]
# 
# pop_a <- pop_mean_data[["mean_a"]][, 1]
# pop_b <- pop_mean_data[["mean_b"]][, 1]
# 
# myplot(0, 0, xlim = c(-100, 100), ylim = c(-100, 100), col = "white")
# mypoints(img_a,
#     img_b,
#     xlim = c(-100, 100),
#     ylim = c(-100, 100),
#     col = hex(as(LAB(50, img_a * 0.65, img_b * 0.65), "RGB")),
#     pch = 1,
#     cex = 1
# )
# mypoints(pop_a,
#     pop_b,
#     col = "black",
#     bg = hex(as(LAB(50, pop_a * 0.65, pop_b * 0.65), "RGB")),
#     pch = 21,
#     cex = 1
# )
# arrows(img_a, img_b, pop_a, pop_b, length = 0.05, lwd = 1)
# 
# myaxis(1, at = seq(-50, 50, 50))
# myaxis(1, at = 100, labels = "100", col.axis = 4)
# myaxis(1, at = -100, labels = "-100", col.axis = 5)
# 
# myaxis(2, at = seq(-50, 50, 50))
# myaxis(2, at = 100, labels = "100", col.axis = 7)
# myaxis(2, at = -100, labels = "-100", col.axis = 6)

# plotting results for ratios of mean and sd. of cone excitations
pdf("../figures/cer_and_sd_results_white_walls.pdf", height = 3, width = 9)
mypar(mfrow=c(1,3))
plotRMCandRSD(pop_mean_data, obs_mean_data, 2, 2, "topleft")
dev.off()


# let's look at what happens with the different mask
data <- read.csv("../data/data_table_white_walls.csv")
data <- data[data$exp_name == "white_walls_filter" & data$mask_name == "white_walls_mask", ]

# apparently obs 3 and 9 are bad observers for the robust ratio model
data <- data[data$obs_name != "03" & data$obs_name != "09", ]

obs_mean_data <- aggregate(data, list(body = data$body_glaven, illum = data$illum_glaven, obs = data$obs_name, hilo = data$hilo_glaven), mean)
obs_mean_data <- subset(obs_mean_data, select = -c(obs, obs_name, body_glaven, illum_glaven, hilo_glaven, exp_name, mask_name))
pop_mean_data <- aggregate(obs_mean_data, list(body_glaven = obs_mean_data$body, illum_glaven = obs_mean_data$illum, hilo_glaven = obs_mean_data$hilo), function(x) {
    c(MN = mean(x), SD = sd(x))
})
pop_mean_data <- subset(pop_mean_data, select = -c(body, illum, hilo))

pop_mean_data$color <- 1

pop_mean_data$pch <- 16

# plotting results for ratios of mean and sd. of cone excitations
pdf("../figures/cer_and_sd_results_white_walls_Special_Mask.pdf", height = 3, width = 9)
mypar(mfrow=c(1,3))
plotRMCandRSD(pop_mean_data, obs_mean_data, 2, 2, "topleft")
dev.off()