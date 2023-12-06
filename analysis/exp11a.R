source("../aux/plotFuncs.R")

data <- read.csv('../data/exp11a_obs_data.csv')

obs.mean.data <- aggregate(. ~ sub_id + img + illum + body + highlow, data, mean)
pop.mean.data <- aggregate(. ~ img + illum + body + highlow, obs.mean.data, function(x) { c(MN=mean(x), SD=sd(x)) })

pop.mean.data$color <- 1

pop.mean.data$pch <- 16

mypalette()

# loop over different background types for the matching stimulus
# for (bg in c("black", "white", "gray")) {
    # myplotting results for whole obj in terms of mean LAB
    # pdf(paste0("../figures/exp11a_mean_lab_", bg,".pdf"), height = 3, width = 9)
    pdf(paste0("../figures/exp11a_mean_lab_gray.pdf"), height = 3, width = 9)
    mypar(mfrow=c(1,3))
    plotLABMean(pop.mean.data, obs.mean.data, "", TRUE)
    dev.off()
# }
