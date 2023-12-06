library("colorscience")

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

# loop over different background types for the matching stimulus
# for (bg in c("black", "white", "gray")) {
    # plotting results for whole obj in terms of mean LAB
    # pdf(paste0("../../figures/exp2_mean_lab_", bg,".pdf"), height = 3, width = 9)
    pdf(paste0("../../figures/exp2_mean_lab_gray.pdf"), height = 3, width = 9)
    mypar(mfrow=c(1,3))
    plotLABMean(pop.mean.data, obs.mean.data, "", FALSE)
    dev.off()
# }

# let us also plot the average responses in LAB space and compare the responses across illumination conditions
mypalette()
pdf("../../figures/exp2_ab_cc.pdf", height = 6, width = 6)
mypar()
myplot(pop.mean.data[["a"]][,1],
       pop.mean.data[["b"]][,1],
       xlim=c(-80, 80),
       ylim=c(-80, 80),
       col = pop.mean.data$color,
       pch = pop.mean.data$pch)

pop_x_blue <- pop.mean.data[pop.mean.data$illum == "blue", ][["a"]][,1]
pop_y_blue <- pop.mean.data[pop.mean.data$illum == "blue", ][["b"]][,1]

pop_x_yellow <- pop.mean.data[pop.mean.data$illum == "yellow", ][["a"]][,1]
pop_y_yellow <- pop.mean.data[pop.mean.data$illum == "yellow", ][["b"]][,1]

line_colors <- pop.mean.data[pop.mean.data$illum == "blue", ]$color
for (cc in 1:length(pop_x_blue)) {
    lines(c(pop_x_blue[cc], pop_x_yellow[cc]), c(pop_y_blue[cc], pop_y_yellow[cc]), col=line_colors[cc])
}

pop_x_blue_sd <- pop.mean.data[pop.mean.data$illum == "blue", ][["a"]][,2]
pop_y_blue_sd <- pop.mean.data[pop.mean.data$illum == "blue", ][["b"]][,2]

pop_x_yellow_sd <- pop.mean.data[pop.mean.data$illum == "yellow", ][["a"]][,2]
pop_y_yellow_sd <- pop.mean.data[pop.mean.data$illum == "yellow", ][["b"]][,2]
theta <- seq(0, 2 * pi, length = 200)
for (cc in 1:length(pop_x_blue)) {
    radius_x <- pop_x_blue_sd[[cc]]
    radius_y <- pop_y_blue_sd[[cc]]
    lines(x = pop_x_blue[[cc]] + radius_x * cos(theta), y = pop_y_blue[[cc]] + radius_y * sin(theta), col=line_colors[cc])
    
    radius_x <- pop_x_yellow_sd[[cc]]
    radius_y <- pop_y_yellow_sd[[cc]]
    lines(x = pop_x_yellow[[cc]] + radius_x * cos(theta), y = pop_y_yellow[[cc]] + radius_y * sin(theta), col=line_colors[cc])
}

rsq <- 1
mytitle(paste0("Adjusted R Squared = ", rsq),
        xlab = "Mean Observer Match a*",
        ylab = "Mean Observer Match b*")
myaxis(1, at = seq(-80, 80, 20))
myaxis(2, at = seq(-80, 80, 20))
# a* axis colors
myaxis(1, at = -80, labels = "-80", col.axis = 5)
myaxis(1, at = 80, labels = "80", col.axis = 4)
# b* axis colors
myaxis(2, at = -80, labels = "-80", col.axis = 6)
myaxis(2, at = 80, labels = "80", col.axis = 7)
dev.off()

# manova to see if l, a, b locations are affected by illuminant, body color, and transparency level
res.man <- manova(cbind(l, a, b) ~ illum * body * highlow, data = obs.mean.data)
summary(res.man)

# let's try computing DE2000 between obs matches under blue and yellow illums and see what we get
bodies <- c("blue", "yellow", "red", "green")
body <- c()
highlow <- c()
sub_id <- c()
trial <- c()
val <- c()
dfc <- 1
for (bc in 1:length(bodies)) {
    obs_blue <- data[data$illum == "blue" & data$body == bodies[[bc]], ]
    obs_yellow <- data[data$illum == "yellow" & data$body == bodies[[bc]], ]
    
    for (cc in 1:nrow(obs_blue)) {
        obs_l_blue <- obs_blue$l[[cc]]
        obs_a_blue <- obs_blue$a[[cc]]
        obs_b_blue <- obs_blue$b[[cc]]
        
        obs_l_yellow <- obs_yellow$l[[cc]]
        obs_a_yellow <- obs_yellow$a[[cc]]
        obs_b_yellow <- obs_yellow$b[[cc]]
        
        lab_blue <- c(obs_l_blue, obs_a_blue, obs_b_blue)
        lab_yellow <- c(obs_l_yellow, obs_a_yellow, obs_b_yellow)
        body[[dfc]] <- bodies[[bc]]
        highlow[[dfc]] <- obs_blue$highlow[[cc]]
        sub_id[[dfc]] <- obs_blue$sub_id[[cc]]
        trial[[dfc]] <- obs_blue$trial[[cc]]
        val[[dfc]] <- deltaE2000(lab_blue, lab_yellow)
        
        dfc <- dfc + 1
    }
}
de2000_vals <- data.frame(body, highlow, sub_id, trial, val)
de2000_vals.avg <- aggregate(. ~ body + highlow, de2000_vals, function(x) { c(MN=mean(x), SD=sd(x)) })