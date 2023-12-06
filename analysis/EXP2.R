source("../alx/plotFuncs.R")

data <- read.csv("../data/exp2_obs_data.csv")

obs.mean.data <- aggregate(. ~ sub_id + img + illum + body + highlow, data, mean)
obs.mean.data <- subset(obs.mean.data, select = -c(sub_id))
pop.mean.data <- aggregate(. ~ img + illum + body + highlow, obs.mean.data, function(x) {
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

mypalette()

# loop over different background types for the matching stimulus
# for (bg in c("black", "white", "gray")) {
# plotting results for whole obj in terms of mean LAB
# pdf(paste0("../figures/exp2_mean_lab_", bg,".pdf"), height = 3, width = 9)
pdf(paste0("../figures/exp2_mean_lab_gray.pdf"), height = 3, width = 9)
mypar(mfrow = c(1, 3))
plotLABMean(pop.mean.data, obs.mean.data, "", FALSE)
dev.off()
# }