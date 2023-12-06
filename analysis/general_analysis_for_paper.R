source("../aux/plotFuncs.R")

library("colorspace")
library("plotrix")

# data <- list(
#     read.csv("../data/exp1_obs_data.csv"),
#     read.csv("../data/exp2_obs_data.csv"),
#     read.csv("../data/exp3_obs_data.csv"),
#     read.csv("../data/exp11a_obs_data.csv"),
#     read.csv("../data/exp11b_obs_data.csv")
# )
data <- read.csv("../data/data_table.csv")

pop_mean_data <- list()
obs_mean_data <- list()
obs_covms_mean <- list()
obs_pas_mean <- list()
obs_covms_max <- list()
obs_pas_max <- list()

exp_names <- c("patch_app", "patch_dye", "filter", "white_walls_patch", "white_walls_filter")
mask_names <- unique(data$mask_name)

for (cc in 1:5) {
    expn <- exp_names[[cc]]
    # d <- data[[cc]]
    d <- data[data$mask_name == "obj_mask" & data$exp_name == expn, ]
    if (expn == "white_walls_filter") {
    	# apparently, these obs had really strange values for robust ratio model, so i
    	# don't trust any of their data anymore
    	d <- d[d$obs_name != "03" & d$obs_name != "09", ]
    }

    obs_data <- data.frame()
    for (obs_name in unique(d$obs_name)) {
        for (illum_glaven in unique(d$illum_glaven)) {
            for (body_glaven in unique(d$body_glaven)) {
                for (hilo_glaven in unique(d$hilo_glaven)) {
                    sub_data <- d[d$obs_name == obs_name & d$illum_glaven == illum_glaven & d$body_glaven == body_glaven & d$hilo_glaven == hilo_glaven, ]
                    if (is.null(sub_data)) {
                        next
                    }

                    # from image
                    t <- data.frame(
                        gwL_glaven.MN = mean(sub_data$gwL_glaven), gwL_glaven.SD = sd(sub_data$gwL_glaven),
                        gwa_glaven.MN = mean(sub_data$gwa_glaven), gwa_glaven.SD = sd(sub_data$gwa_glaven),
                        gwb_glaven.MN = mean(sub_data$gwb_glaven), gwb_glaven.SD = sd(sub_data$gwb_glaven)
                    )

                    t <- cbind(t, data.frame(illum_glaven = illum_glaven, body_glaven = body_glaven, hilo_glaven = hilo_glaven))

                    # from obs
                    t <- cbind(t, data.frame(
                        gwL_match.MN = mean(sub_data$gwL_match), gwL_match.SD = sd(sub_data$gwL_match),
                        gwa_match.MN = mean(sub_data$gwa_match), gwa_match.SD = sd(sub_data$gwa_match),
                        gwb_match.MN = mean(sub_data$gwb_match), gwb_match.SD = sd(sub_data$gwb_match)
                    ))

                    if (expn == "filter" | expn == "white_walls_filter") { # cc == 3 | cc == 5) {
                        t <- cbind(t, data.frame(
                            wp_top5L_filter.MN = mean(sub_data$wp_top5L_filter), wp_top5L_filter.SD = sd(sub_data$wp_top5L_filter),
                            wp_top5a_filter.MN = mean(sub_data$wp_top5a_filter), wp_top5a_filter.SD = sd(sub_data$wp_top5a_filter),
                            wp_top5b_filter.MN = mean(sub_data$wp_top5b_filter), wp_top5b_filter.SD = sd(sub_data$wp_top5b_filter)
                        ))
                    }

                    obs_data <- rbind(obs_data, t)
                }
            }
        }
    }

    ct <- 1
    covms_mean <- list()
    pas_mean <- list()
    covms_max <- list()
    pas_max <- list()
    tmp_data <- data.frame()
    for (illum_glaven in unique(obs_data$illum_glaven)) {
        for (body_glaven in unique(obs_data$body_glaven)) {
            for (hilo_glaven in unique(obs_data$hilo_glaven)) {
                sub_data <- obs_data[obs_data$illum_glaven == illum_glaven & obs_data$body_glaven == body_glaven & obs_data$hilo_glaven == hilo_glaven, ]
                if (is.null(sub_data)) {
                    next
                }

                # from image
                t <- data.frame(
                    gwL_glaven.MN = mean(sub_data$gwL_glaven.MN), gwL_glaven.SD = sd(sub_data$gwL_glaven.MN),
                    gwa_glaven.MN = mean(sub_data$gwa_glaven.MN), gwa_glaven.SD = sd(sub_data$gwa_glaven.MN),
                    gwb_glaven.MN = mean(sub_data$gwb_glaven.MN), gwb_glaven.SD = sd(sub_data$gwb_glaven.MN)
                )

                t <- cbind(t, data.frame(illum_glaven = illum_glaven, body_glaven = body_glaven, hilo_glaven = hilo_glaven))

                # from obs
                t <- cbind(t, data.frame(
                    gwL_match.MN = mean(sub_data$gwL_match.MN), gwL_match.SD = sd(sub_data$gwL_match.MN),
                    gwa_match.MN = mean(sub_data$gwa_match.MN), gwa_match.SD = sd(sub_data$gwa_match.MN),
                    gwb_match.MN = mean(sub_data$gwb_match.MN), gwb_match.SD = sd(sub_data$gwb_match.MN)
                ))

                covms_mean[[ct]] <- cov(data.frame(sub_data["gwa_match.MN"], sub_data["gwb_match.MN"]))
                pas_mean[[ct]] <- prcomp(cbind(sub_data["gwa_match.MN"], sub_data["gwb_match.MN"]), center = T, scale = T)
                if (expn == "filter" | expn == "white_walls_filter") { # (cc == 3 | cc == 5) {
                    t <- cbind(t, data.frame(
                        wp_top5L_filter.MN = mean(sub_data$wp_top5L_filter.MN), wp_top5L_filter.SD = sd(sub_data$wp_top5L_filter.MN),
                        wp_top5a_filter.MN = mean(sub_data$wp_top5a_filter.MN), wp_top5a_filter.SD = sd(sub_data$wp_top5a_filter.MN),
                        wp_top5b_filter.MN = mean(sub_data$wp_top5b_filter.MN), wp_top5b_filter.SD = sd(sub_data$wp_top5b_filter.MN)
                    ))

                    covms_max[[ct]] <- cov(data.frame(sub_data["wp_top5a_filter.MN"], sub_data["wp_top5b_filter.MN"]))
                    pas_max[[ct]] <- prcomp(cbind(sub_data["wp_top5a_filter.MN"], sub_data["wp_top5b_filter.MN"]), center = T, scale = T)
                }

                tmp_data <- rbind(tmp_data, t)

                ct <- ct + 1
            }
        }
    }

    # if (cc != 2) {
    #     obs_data <- aggregate(. ~ obs_name + img + illum_glaven + body_glaven + hilo_glaven, data[[cc]], mean)
    #     obs_data <- subset(obs_data, select = -c(obs_name))
    #     tmp_data <- aggregate(. ~ img + illum_glaven + body_glaven + hilo_glaven, obs_data, function(x) {
    #         c(MN = mean(x), SD = sd(x))
    #     })
    # } else {
    #     obs_data <- aggregate(. ~ obs_name + img + illum_glaven + body_glaven + hilo_glaven, data[[cc]], mean)
    #     obs_data <- subset(obs_data, select = -c(obs_name))
    #     tmp_data <- aggregate(. ~ img + illum_glaven + body_glaven + hilo_glaven, obs_data, function(x) {
    #         c(MN = mean(x), SD = sd(x))
    #     })
    # }

    # covms <- list()
    # ct <- 1
    # for(img in unique(obs_data$img)) {
    #     d <- obs_data[obs_data$img == img, ]
    #     if(is.null(d)) {
    #         next
    #     }
    #     if (cc == 3 | cc == 5) {
    #         covm <- cov(data.frame(d["gwa_match"], d["gwb_glaven"]))
    #     } else {
    #         covm <- cov(data.frame(d["a"], d["b"]))
    #     }
    #     covms[[ct]] <- list(img = img, m = covm)
    #     ct <- ct + 1
    # }

    tmp_data$color <- NA

    tmp_data[tmp_data$body_glaven == "red" & tmp_data$hilo_glaven == 0, ]$color <- 8
    tmp_data[tmp_data$body_glaven == "green" & tmp_data$hilo_glaven == 0, ]$color <- 9
    tmp_data[tmp_data$body_glaven == "blue" & tmp_data$hilo_glaven == 0, ]$color <- 10
    tmp_data[tmp_data$body_glaven == "yellow" & tmp_data$hilo_glaven == 0, ]$color <- 11

    if (expn != "white_walls_patch") {
        tmp_data[tmp_data$body_glaven == "red" & tmp_data$hilo_glaven == 1, ]$color <- 4
        tmp_data[tmp_data$body_glaven == "green" & tmp_data$hilo_glaven == 1, ]$color <- 5
        tmp_data[tmp_data$body_glaven == "blue" & tmp_data$hilo_glaven == 1, ]$color <- 6
        tmp_data[tmp_data$body_glaven == "yellow" & tmp_data$hilo_glaven == 1, ]$color <- 14
    }

    tmp_data$pch <- NA

    if (cc == 4 || cc == 5) {
        tmp_data$pch <- 19
    } else {
        tmp_data[tmp_data$illum_glaven == "blue", ]$pch <- 19
        tmp_data[tmp_data$illum_glaven == "yellow", ]$pch <- 6
    }

    pop_mean_data[[cc]] <- as.data.frame(tmp_data)
    obs_mean_data[[cc]] <- as.data.frame(obs_data)
    obs_covms_mean[[cc]] <- covms_mean
    obs_pas_mean[[cc]] <- pas_mean
    obs_covms_max[[cc]] <- covms_max
    obs_pas_max[[cc]] <- pas_max
}

rg <- list()

mypalette()

# plotting results for whole obj in terms of mean LAB
pdf("../figures/all_exps_mean_lab.pdf", height = 7, width = 11)
mypar(mfrow = c(2, 3))

# for (i in 1:4) { # doesn't actually access data from set 4. see hardcoded
# below
cols <- c(1, 6, 4, 5)
for (d in c("L", "a", "b")) {
    for (i in 1:2) {
        pop_x <- pop_mean_data[[i]][[paste0("gw", d, "_glaven.MN")]]
        pop_y <- pop_mean_data[[i]][[paste0("gw", d, "_match.MN")]]
        pop_sd <- pop_mean_data[[i]][[paste0("gw", d, "_match.SD")]]

        obs_x <- obs_mean_data[[i]][[paste0("gw", d, "_glaven.MN")]]
        obs_y <- obs_mean_data[[i]][[paste0("gw", d, "_match.MN")]]
        plotch <- pop_mean_data[[i]]$pch

        if (d == "L") {
            minlim <- 0
            maxlim <- 100
        } else {
            minlim <- -100
            maxlim <- 100
        }

        if (i == 1) {
            myplot(pop_x,
                pop_y,
                xlim = c(minlim, maxlim),
                ylim = c(minlim, maxlim),
                col = cols[i],
                pch = plotch
            )
        } else {
            mypoints(pop_x,
                pop_y,
                col = cols[i],
                pch = plotch
            )
        }
        my_V_ebar(pop_x, pop_y, pop_sd, col = cols[i])
        abline(rg[[i]] <- lm(obs_y ~ obs_x), lty = 2, col = cols[i])
    }
    abline(0, 1)
    if (d == "L") {
        mytitle(xlab = paste0("Mean Image ", toupper(d), "*"), ylab = paste0("Observer Match ", toupper(d), "*"))
    } else {
        mytitle(xlab = paste0("Mean Image ", d, "*"), ylab = paste0("Observer Match ", d, "*"))
    }

    if (d == "L") {
        myaxis(1, at = seq(0, 100, 20))
        myaxis(2, at = seq(0, 100, 20))
        mylegend("bottomright",
            c("Blue Illuminant", "White Illuminant", "Mean Proximal Match", "Mean Dye Match"),
            pch = c(19, 6, NA, NA, NA, NA),
            lty = c(NA, NA, 2, 2, 2, 2),
            col = c(1, 1, 1, 6, 4, 5)
        )
    } else if (d == "a") {
        myaxis(1, at = seq(-50, 50, 50))
        myaxis(2, at = seq(-50, 50, 50))
        myaxis(1, at = -100, labels = "-100", col.axis = 5)
        myaxis(1, at = 100, labels = "100", col.axis = 4)
        myaxis(2, at = -100, labels = "-100", col.axis = 5)
        myaxis(2, at = 100, labels = "100", col.axis = 4)
    } else if (d == "b") {
        myaxis(1, at = seq(-50, 50, 50))
        myaxis(2, at = seq(-50, 50, 50))
        myaxis(1, at = -100, labels = "-100", col.axis = 7)
        myaxis(1, at = 100, labels = "100", col.axis = 6)
        myaxis(2, at = -100, labels = "-100", col.axis = 7)
        myaxis(2, at = 100, labels = "100", col.axis = 6)
    }
}

for (d in c("L", "a", "b")) {
    for (i in 3:4) {
        if (i == 3) {
            pop_x <- pop_mean_data[[i]][[paste0("gw", d, "_glaven.MN")]]
            pop_y <- pop_mean_data[[i]][[paste0("gw", d, "_match.MN")]]
            pop_sd <- pop_mean_data[[i]][[paste0("gw", d, "_match.SD")]]

            obs_x <- obs_mean_data[[i]][[paste0("gw", d, "_glaven.MN")]]
            obs_y <- obs_mean_data[[i]][[paste0("gw", d, "_match.MN")]]
            plotch <- pop_mean_data[[i]]$pch
        } else {
            pop_x <- pop_mean_data[[3]][[paste0("gw", d, "_glaven.MN")]]
            pop_y <- pop_mean_data[[3]][[paste0("wp_top5", d, "_filter.MN")]]
            pop_sd <- pop_mean_data[[3]][[paste0("wp_top5", d, "_filter.SD")]]

            obs_x <- obs_mean_data[[3]][[paste0("gw", d, "_glaven.MN")]]
            obs_y <- obs_mean_data[[3]][[paste0("wp_top5", d, "_filter.MN")]]
            plotch <- pop_mean_data[[3]]$pch
        }

        if (d == "L") {
            minlim <- 0
            maxlim <- 100
        } else {
            minlim <- -100
            maxlim <- 100
        }

        if (i == 3) {
            myplot(pop_x,
                pop_y,
                xlim = c(minlim, maxlim),
                ylim = c(minlim, maxlim),
                col = cols[i],
                pch = plotch
            )
        } else {
            mypoints(pop_x,
                pop_y,
                col = cols[i],
                pch = plotch
            )
        }
        my_V_ebar(pop_x, pop_y, pop_sd, col = cols[i])
        abline(rg[[i]] <- lm(obs_y ~ obs_x), lty = 2, col = cols[i])
    }
    abline(0, 1)
    if (d == "L") {
        mytitle(xlab = paste0("Mean Image ", toupper(d), "*"), ylab = paste0("Observer Match ", toupper(d), "*"))
    } else {
        mytitle(xlab = paste0("Mean Image ", d, "*"), ylab = paste0("Observer Match ", d, "*"))
    }

    if (d == "L") {
        myaxis(1, at = seq(0, 100, 20))
        myaxis(2, at = seq(0, 100, 20))
        mylegend("bottomright",
            c("Blue Illuminant", "White Illuminant", "Mean Filter Match", "White Point Filter Match"),
            pch = c(19, 6, NA, NA, NA, NA),
            lty = c(NA, NA, 2, 2, 2, 2),
            col = c(1, 1, 4, 5)
        )
    } else if (d == "a") {
        myaxis(1, at = seq(-50, 50, 50))
        myaxis(2, at = seq(-50, 50, 50))
        myaxis(1, at = -100, labels = "-100", col.axis = 5)
        myaxis(1, at = 100, labels = "100", col.axis = 4)
        myaxis(2, at = -100, labels = "-100", col.axis = 5)
        myaxis(2, at = 100, labels = "100", col.axis = 4)
    } else if (d == "b") {
        myaxis(1, at = seq(-50, 50, 50))
        myaxis(2, at = seq(-50, 50, 50))
        myaxis(1, at = -100, labels = "-100", col.axis = 7)
        myaxis(1, at = 100, labels = "100", col.axis = 6)
        myaxis(2, at = -100, labels = "-100", col.axis = 7)
        myaxis(2, at = 100, labels = "100", col.axis = 6)
    }
}

### bar charts showing R^2

# mybarplot(c(summary(prox_r[[1]])$adj.r.squared,
#     summary(dye_r[[1]])$adj.r.squared,
#     summary(mean_filt_r[[1]])$adj.r.squared,
#     summary(max_filt_r[[1]])$adj.r.squared),
#     col=c(1, 6, 4, 5),
#     ylim=c(0,1))
# mytitle(ylab="Adjusted R Squared")
#
# mybarplot(c(summary(prox_r[[2]])$adj.r.squared,
#     summary(dye_r[[2]])$adj.r.squared,
#     summary(mean_filt_r[[2]])$adj.r.squared,
#     summary(max_filt_r[[2]])$adj.r.squared),
#     col=c(1, 6, 4, 5),
#     names.arg=c("Prox", "Dye", "Mean Filt", "Max Filt"), ylim=c(0,1))
# mytitle(xlab="Experimental Condition")
#
# mybarplot(c(summary(prox_r[[3]])$adj.r.squared,
#     summary(dye_r[[3]])$adj.r.squared,
#     summary(mean_filt_r[[3]])$adj.r.squared,
#     summary(max_filt_r[[3]])$adj.r.squared),
#     col=c(1, 6, 4, 5),
#     ylim=c(0,1))
dev.off()

# make a plot showing the color matches in a chromaticity diagram
pdf("../figures/chrom_diagram_exps_1_2_3.pdf", height = 6.5, width = 6)
exps <- c("Patch proximal match", "Patch dye match", "Flat filter - mean color", "Flat filter - White Point")
mypar(mfrow = c(2, 2), pty = "s")
for (i in 1:4) {
    if (i < 3) {
        img_a <- pop_mean_data[[i]][["gwa_glaven.MN"]]
        img_b <- pop_mean_data[[i]][["gwb_glaven.MN"]]

        pop_a <- pop_mean_data[[i]][["gwa_match.MN"]]
        pop_b <- pop_mean_data[[i]][["gwb_match.MN"]]

        obs_a <- obs_mean_data[[i]][["gwa_match.MN"]]
        obs_b <- obs_mean_data[[i]][["gwb_match.MN"]]

        pop_a_sd <- pop_mean_data[[i]][["gwa_match.SD"]]
        pop_b_sd <- pop_mean_data[[i]][["gwb_match.SD"]]

        the_covm <- obs_covms_mean[[i]]
        the_pas <- obs_pas_mean[[i]]
    } else if (i == 3) {
        img_a <- pop_mean_data[[i]][["gwa_glaven.MN"]]
        img_b <- pop_mean_data[[i]][["gwb_glaven.MN"]]

        pop_a <- pop_mean_data[[i]][["gwa_match.MN"]]
        pop_b <- pop_mean_data[[i]][["gwb_match.MN"]]

        obs_a <- obs_mean_data[[i]][["gwa_match.MN"]]
        obs_b <- obs_mean_data[[i]][["gwb_match.MN"]]

        pop_a_sd <- pop_mean_data[[i]][["gwa_match.SD"]]
        pop_b_sd <- pop_mean_data[[i]][["gwb_match.SD"]]

        the_covm <- obs_covms_mean[[i]]
        the_pas <- obs_pas_mean[[i]]
    } else {
        img_a <- pop_mean_data[[3]][["gwa_glaven.MN"]]
        img_b <- pop_mean_data[[3]][["gwb_glaven.MN"]]

        pop_a <- pop_mean_data[[3]][["wp_top5a_filter.MN"]]
        pop_b <- pop_mean_data[[3]][["wp_top5b_filter.MN"]]

        obs_a <- obs_mean_data[[3]][["wp_top5a_filter.MN"]]
        obs_b <- obs_mean_data[[3]][["wp_top5b_filter.MN"]]

        pop_a_sd <- pop_mean_data[[3]][["wp_top5a_filter.SD"]]
        pop_b_sd <- pop_mean_data[[3]][["wp_top5b_filter.SD"]]

        the_covm <- obs_covms_mean[[3]]
        the_pas <- obs_pas_max[[3]]
    }

    myplot(0, 0, xlim = c(-100, 100), ylim = c(-100, 100), col = "white", cex = 1)
    # plotrix::draw.ellipse(pop_a, pop_b, a = pop_a_sd, b = pop_b_sd)
    # dataEllipse(obs_a, obs_b, add = TRUE, xlim = c(-100, 100), ylim = c(-100, 100), levels = c(0.68), center.pch = FALSE)
    for (cc in 1:length(the_pas)) {
        lines(ellipse::ellipse(the_covm[[cc]], centre = c(pop_a[[cc]], pop_b[[cc]]), level = 0.68), col='black')
        # plotrix::draw.ellipse(pop_a[[cc]], pop_b[[cc]], a = the_pas[[cc]]$scale[[1]], b = the_pas[[cc]]$scale[[2]], angle = atan2(the_pas[[cc]]$rotation[2, 1], the_pas[[cc]]$rotation[1, 1]) * 180 / pi)
    }
    cols <- as(LAB(50, img_a * 0.7, img_b * 0.7), "RGB")
    cols <- coords(cols)
    cols[cols < 0] <- 0
    cols[cols > 1] <- 1
    cols <- RGB(cols[, 1], cols[, 2], cols[, 3])
    mypoints(img_a,
        img_b,
        xlim = c(-100, 100),
        ylim = c(-100, 100),
        col = hex(cols),
        pch = 1,
        cex = 1
    )
    cols <- as(LAB(50, pop_a * 0.7, pop_b * 0.7), "RGB")
    cols <- coords(cols)
    cols[cols < 0] <- 0
    cols[cols > 1] <- 1
    cols <- RGB(cols[, 1], cols[, 2], cols[, 3])
    mypoints(pop_a,
        pop_b,
        col = "black",
        bg = hex(cols),
        pch = 21,
        cex = 1
    )

    arrows(img_a, img_b, pop_a, pop_b, length = 0.05, lwd = 1)

    mytitle(xlab = "a*", ylab = "b*", main = exps[[i]])

    myaxis(1, at = seq(-50, 50, 50))
    myaxis(1, at = 100, labels = "100", col.axis = 4)
    myaxis(1, at = -100, labels = "-100", col.axis = 5)

    myaxis(2, at = seq(-50, 50, 50))
    myaxis(2, at = 100, labels = "100", col.axis = 7)
    myaxis(2, at = -100, labels = "-100", col.axis = 6)
}
dev.off()

# plotting results for white bkg and illum_glaven in terms of mean LAB for uniform
# patch
pdf("../figures/white_wall_results_exp11a.pdf", height = 3, width = 9)
pop_mean_data[[4]]$color <- 1
pop_mean_data[[4]]$pch <- 16
obs_mean_data[[4]]$color <- 1
obs_mean_data[[4]]$pch <- 16
mypar(mfrow = c(1, 3))
plotLABMean(pop_mean_data[[4]], obs_mean_data[[4]], TRUE)
dev.off()

# plotting results for white bkg and illum_glaven in terms of mean LAB for filter
pdf("../figures/white_wall_results_exp11b.pdf", height = 3, width = 9)
pop_mean_data[[5]]$color <- 1
pop_mean_data[[5]]$pch <- 16
obs_mean_data[[5]]$color <- 1
obs_mean_data[[5]]$pch <- 16
mypar(mfrow = c(1, 3))
plotLABMean(pop_mean_data[[5]], obs_mean_data[[5]], TRUE)
dev.off()

# try removing the one annoying point to see if a* improves to be like before
# (no, doesn't really improve anything)
pdf("../figures/white_wall_results_exp11b_wo_bad_a.pdf", height = 3, width = 9)
bad_a <- pop_mean_data[[5]]$gwa_glaven.MN[[2]]
pop_mean_wwalls_wo_bad_a <- pop_mean_data[[5]]
pop_mean_wwalls_wo_bad_a <- pop_mean_wwalls_wo_bad_a[pop_mean_wwalls_wo_bad_a$gwa_glaven.MN != bad_a, ]
obs_mean_wwalls_wo_bad_a <- obs_mean_data[[5]]
obs_mean_wwalls_wo_bad_a <- obs_mean_wwalls_wo_bad_a[obs_mean_wwalls_wo_bad_a$gwa_glaven.MN != bad_a, ]
pop_mean_wwalls_wo_bad_a$color <- 1
pop_mean_wwalls_wo_bad_a$pch <- 16
obs_mean_wwalls_wo_bad_a$color <- 1
obs_mean_wwalls_wo_bad_a$pch <- 16
mypar(mfrow = c(1, 3))
plotLABMean(pop_mean_wwalls_wo_bad_a, obs_mean_wwalls_wo_bad_a, TRUE)
dev.off()


# white walls - flat filter
# make a plot showing the color matches in a chromaticity diagram
pdf("../figures/chrom_diagram_exp11b.pdf", height = 6.5, width = 6)
mypar(pty = "s")
img_a <- pop_mean_data[[5]][["gwa_glaven.MN"]]
img_b <- pop_mean_data[[5]][["gwb_glaven.MN"]]

pop_a <- pop_mean_data[[5]][["gwa_match.MN"]]
pop_b <- pop_mean_data[[5]][["gwb_match.MN"]]

obs_a <- obs_mean_data[[5]][["gwa_match.MN"]]
obs_b <- obs_mean_data[[5]][["gwb_match.MN"]]

pop_a_sd <- pop_mean_data[[5]][["gwa_match.SD"]]
pop_b_sd <- pop_mean_data[[5]][["gwb_match.SD"]]

#the_covm <- obs_covms_mean[[5]]
the_pas <- obs_pas_mean[[5]]

myplot(0, 0, xlim = c(-100, 100), ylim = c(-100, 100), col = "white", cex = 1)
# plotrix::draw.ellipse(pop_a, pop_b, a = pop_a_sd, b = pop_b_sd)
# dataEllipse(obs_a, obs_b, add = TRUE, xlim = c(-100, 100), ylim = c(-100, 100), levels = c(0.68), center.pch = FALSE)
for (cc in 1:length(the_pas)) {
	lines(ellipse::ellipse(the_covm[[cc]], centre = c(pop_a[[cc]], pop_b[[cc]]), level = 0.68), col='black')
	# plotrix::draw.ellipse(pop_a[[cc]], pop_b[[cc]], a = the_pas[[cc]]$scale[[1]], b = the_pas[[cc]]$scale[[2]], angle = atan2(the_pas[[cc]]$rotation[2, 1], the_pas[[cc]]$rotation[1, 1]) * 180 / pi)
}
cols <- as(LAB(50, img_a * 0.7, img_b * 0.7), "RGB")
cols <- coords(cols)
cols[cols < 0] <- 0
cols[cols > 1] <- 1
cols <- RGB(cols[, 1], cols[, 2], cols[, 3])
mypoints(img_a,
	img_b,
	xlim = c(-100, 100),
	ylim = c(-100, 100),
	col = hex(cols),
	pch = 1,
	cex = 1
)
cols <- as(LAB(50, pop_a * 0.7, pop_b * 0.7), "RGB")
cols <- coords(cols)
cols[cols < 0] <- 0
cols[cols > 1] <- 1
cols <- RGB(cols[, 1], cols[, 2], cols[, 3])
mypoints(pop_a,
	pop_b,
	col = "black",
	bg = hex(cols),
	pch = 21,
	cex = 1
)

arrows(img_a, img_b, pop_a, pop_b, length = 0.05, lwd = 1)

mytitle(xlab = "a*", ylab = "b*", main = "Flat filter - mean color - White Walls")

myaxis(1, at = seq(-50, 50, 50))
myaxis(1, at = 100, labels = "100", col.axis = 4)
myaxis(1, at = -100, labels = "-100", col.axis = 5)

myaxis(2, at = seq(-50, 50, 50))
myaxis(2, at = 100, labels = "100", col.axis = 7)
myaxis(2, at = -100, labels = "-100", col.axis = 6)
dev.off()