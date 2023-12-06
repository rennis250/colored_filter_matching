source("../aux/plotFuncs.R")

# data <- read.csv("../data/exp1_obs_data.csv")
data <- read.csv("../data/data_table.csv")
data <- data[data$exp_name == "patch_dye" & data$mask_name == "obj_mask", ]

# obs.mean.data <- aggregate(. ~ sub_id + img + illum + body + highlow, data, mean)
# obs.mean.data <- subset(obs.mean.data, select = -c(sub_id))
# pop.mean.data <- aggregate(. ~ img + illum + body + highlow, obs.mean.data, function(x) {
#     c(MN = mean(x), SD = sd(x))
# })

obs.mean.data <- aggregate(data, list(body = data$body_glaven,
		illum = data$illum_glaven,
		sub_id = data$obs_name,
		lighter_darker = data$hilo_glaven),
	mean)
obs.mean.data <- subset(obs.mean.data, select = -c(img_name,
		obs_name,
		sub_id,
		body_glaven,
		illum_glaven,
		hilo_glaven))

pop.mean.data <- aggregate(obs.mean.data, list(illum = obs.mean.data$illum,
		body = obs.mean.data$body,
		highlow = obs.mean.data$lighter_darker),
		function(x) {
			c(MN = mean(x), SD = sd(x))
		})

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

# plotting results for whole obj in terms of mean LAB

mypalette()

# pdf("../figures/exp1_mean_lab.pdf", height = 3, width = 9)
# mypar(mfrow = c(1, 3))
# plotLABMean(pop.mean.data, obs.mean.data, FALSE)
# dev.off()

# pdf("../figures/exp1_mean_lms.pdf", height = 3, width = 9)
# mypar(mfrow = c(1, 3))
# plotLMSMean(pop.mean.data, obs.mean.data, "cone_", FALSE)
# dev.off()

# pdf("../figures/exp1_mean_dkl.pdf", height = 3, width = 9)
# mypar(mfrow = c(1, 3))
# plotDKLMean(pop.mean.data, obs.mean.data, "", FALSE)
# dev.off()

# let's also make DE2000 maps for each grand mean setting

for (ic in 1:nrow(pop.mean.data)) {
    l_or_d <- ""
    l_or_d_trans <- "_"
    if (pop.mean.data$highlow[ic] == 1) {
        l_or_d <- "lighter"
    } else {
        l_or_d <- "lightest"
        l_or_d_trans <- "_lightest_"
    }

    img_name <- paste0(
        "../images/exp1/mitsuba_caustics_rgb_trans",
        l_or_d_trans,
        "illum_",
        pop.mean.data$illum[ic],
        "_refl_munsell_",
        pop.mean.data$body[ic],
        "_",
        l_or_d,
        ".png"
    )

    system(paste0(
        "de2000_map -monName=eizo",
        " -chromaFile=../calibration/exp1.chroma.csv",
        " -gammaFile=../calibration/exp1.gamma.csv",
        " -maskDir=../images/masks",
        " -anchorColorL=", as.character(pop.mean.data$gwL_match[[ic, 1]]),
        " -anchorColorA=", as.character(pop.mean.data$gwa_match[[ic, 1]]),
        " -anchorColorB=", as.character(pop.mean.data$gwb_match[[ic, 1]]),
        " ", img_name
    ))
}

system("mv de2000_map*.png ../images/de2000_maps/exp1/")