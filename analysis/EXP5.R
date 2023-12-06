source("../aux/plotFuncs.R")

data <- read.csv("../images/exp5_rmc_rsd_comp/exp5_rmc_rsd_comp_chrom_stats.csv")

data$color <- NA

for (ic in 1:length(data$img_name)) {
    fparts <- strsplit(as.character(data$img_name[ic]), "_")
    if (fparts[[1]][6] == "1") {
        data$color[ic] <- 6
    } else if (fparts[[1]][6] == "2") {
        data$color[ic] <- 7
    } else if (fparts[[1]][6] == "3") {
        data$color[ic] <- 4
    } else if (fparts[[1]][6] == "4") {
        data$color[ic] <- 5
    }
}

data$pch <- NA

for (ic in 1:length(data$img_name)) {
    fparts <- strsplit(as.character(data$img_name[ic]), "_")
    bkgd <- strsplit(fparts[[1]][length(fparts[[1]])], "p")
    l_or_d <- fparts[[1]][12]
    if (l_or_d == "1") {
        if (bkgd[[1]][1] == "blue.") {
            data$pch[ic] <- 0
        } else if (bkgd[[1]][1] == "yellow.") {
            data$pch[ic] <- 1
        } else if (bkgd[[1]][1] == "red.") {
            data$pch[ic] <- 2
        } else if (bkgd[[1]][1] == "green.") {
            data$pch[ic] <- 5
        } else if (bkgd[[1]][1] == "dist.") {
            data$pch[ic] <- 3
        }
    } else if (l_or_d == "2") {
        if (bkgd[[1]][1] == "blue") {
            data$pch[ic] <- 15
        } else if (bkgd[[1]][1] == "yellow.") {
            data$pch[ic] <- 16
        } else if (bkgd[[1]][1] == "red.") {
            data$pch[ic] <- 17
        } else if (bkgd[[1]][1] == "green.") {
            data$pch[ic] <- 18
        } else if (bkgd[[1]][1] == "dist.") {
            data$pch[ic] <- 4
        }
    }
}

mypalette()

pdf("../figures/rmc_rsd_from_simulation.pdf", height = 3, width = 9)
#png("../figures/rmc_rsd_from_simulation.png", height = 300, width = 900)
mypar(mfrow = c(1, 3))
for (d in c("l", "m", "s")) {
    myplot(data[[paste0("rmc_", d)]],
        data[[paste0("rsd_", d)]],
        xlim = c(0, 1),
        ylim = c(0, 3.5),
        col = data$color,
        pch = data$pch
    )
    abline(0, 1)
    myaxis(1, at = seq(0, 1, 1 / 5))
    myaxis(2, at = seq(0, 3.5, 3.5 / 5))
    mytitle(toupper(d), xlab = "Ratio Mean Cone Exc.", ylab = "Ratio Sta. Dev. Cone Exc.")
}
dev.off()

# plot just the five that we chose

imgs <- Sys.glob("../exp5_rmc_rsd_comp/stimuli/images/*.png")
idxs <- c()
ic <- 1
for (i in imgs) {
    if (substr(i, nchar(i) - 11, nchar(i)) != "brighter.png") {
        the.name <- strsplit(i, "/")
        idxs[ic] <- which(data["img_name"] == the.name[[1]][5])
        ic <- ic + 1
    }
}

data.the.five <- data[idxs, ]

pdf("../figures/rmc_rsd_from_simulation_the_chosen_five.pdf", height = 3, width = 9)
mypar(mfrow = c(1, 3))
for (d in c("l", "m", "s")) {
    myplot(data.the.five[[paste0("rmc_", d)]],
        data.the.five[[paste0("rsd_", d)]],
        xlim = c(0, 1),
        ylim = c(0, 3.5),
        col = data.the.five$color,
        pch = data.the.five$pch
    )
    abline(0, 1)
    myaxis(1, at = seq(0, 1, 1 / 5))
    myaxis(2, at = seq(0, 3.5, 3.5 / 5))
    mytitle(toupper(d), xlab = "Ratio Mean Cone Exc.", ylab = "Ratio Sta. Dev. Cone Exc.")
}
dev.off()

# let's plot data from experiment with the chosen five and their brighter counterparts

# data.from.exp <- read.csv("../data/exp5_obs_data.csv")
data <- read.csv("../data/data_table.csv")
data.from.exp <- data[data$exp_name == "rmc_rsd_comp" & data$mask_name == "obj_mask", ]

# exclude me, to clean up messiness of data

data.from.exp <- data.from.exp[data.from.exp["obs_name"] != "re", ]

obs.mean.data <- aggregate(data.from.exp, list(img_name = data.from.exp$img_name, rg = data.from.exp$rg_glaven, by = data.from.exp$by_glaven, bkgd = data.from.exp$bkgd_glaven, illum = data.from.exp$illum_glaven, sub_id = data.from.exp$obs_name, lighter_darker = data.from.exp$hilo_glaven), mean)
obs.mean.data <- subset(obs.mean.data, select = -c(obs_name, rg_glaven, by_glaven, bkgd_glaven, illum_glaven, hilo_glaven, exp_name, mask_name))
# obs.mean.data <- aggregate(. ~ sub_id + img + img_name + illum + rg + by + lighter_darker + bkgd, data.from.exp, mean)

# look at observers that did experiment with EXTREMMEEE distributions

sub_ids <- obs.mean.data["sub_id"]
obs.mean.dataEXTREMMEEE <- obs.mean.data[sub_ids == "05" | sub_ids == "re2" | sub_ids == "08a", ]
obs.mean.dataEXTREMMEEE <- subset(obs.mean.dataEXTREMMEEE, select = -c(sub_id))
# pop.mean.dataEXTREMMEEE <- aggregate(. ~ img + img_name + illum + rg + by + lighter_darker + bkgd, obs.mean.dataEXTREMMEEE, function(x) {
#    c(MN = mean(x), SD = sd(x))
# })
pop.mean.dataEXTREMMEEE <- aggregate(obs.mean.dataEXTREMMEEE, list(img_name = obs.mean.dataEXTREMMEEE$img_name, rg = obs.mean.dataEXTREMMEEE$rg, by = obs.mean.dataEXTREMMEEE$by, lighter_darker = obs.mean.dataEXTREMMEEE$lighter_darker, illum = obs.mean.dataEXTREMMEEE$illum, bkgd = obs.mean.dataEXTREMMEEE$bkgd), function(x) {
    c(MN = mean(x), SD = sd(x))
})

pdf("../figures/rmc_rsd_comp_EXTREMMEEE.pdf", height = 3, width = 9)
mypar(mfrow = c(1, 3))
plotRMCandRSD(pop.mean.dataEXTREMMEEE, obs.mean.dataEXTREMMEEE, 3, 3, "topleft", "Adjusted R Squared =      ", 0.6)
dev.off()

# now, let's zoom in and paint the interesting image in red

pdf("../figures/rmc_rsd_comp_ZOOM.pdf", height = 3, width = 9)
maxx <- 1.2
maxy <- 1.2
mypar(mfrow = c(1, 3))
for (i in c("l", "m", "s")) {
    pop_rmc_x <- pop.mean.dataEXTREMMEEE[[paste0("rmc_", i, "_glaven")]][, 1]
    pop_rmc_y <- pop.mean.dataEXTREMMEEE[[paste0("rmc_", i, "_match")]][, 1]
    pop_rsd_x <- pop.mean.dataEXTREMMEEE[[paste0("rsd_", i, "_glaven")]][, 1]
    pop_rsd_y <- pop.mean.dataEXTREMMEEE[[paste0("rsd_", i, "_filter")]][, 1]
    pop_rmc_sd <- pop.mean.dataEXTREMMEEE[[paste0("rmc_", i, "_match")]][, 2]
    pop_rsd_sd <- pop.mean.dataEXTREMMEEE[[paste0("rsd_", i, "_filter")]][, 2]

    obs_rmc_x <- obs.mean.dataEXTREMMEEE[[paste0("rmc_", i, "_glaven")]]
    obs_rmc_y <- obs.mean.dataEXTREMMEEE[[paste0("rmc_", i, "_match")]]
    obs_rsd_x <- obs.mean.dataEXTREMMEEE[[paste0("rsd_", i, "_glaven")]]
    obs_rsd_y <- obs.mean.dataEXTREMMEEE[[paste0("rsd_", i, "_filter")]]

    red_id <- pop.mean.dataEXTREMMEEE$img_name == "mitsuba_comp_rmc_rsd_illum_4_rg_3_by_1_ld_2_bkgd_voronoi_diff_dist.png"

    myplot(pop_rmc_x,
        pop_rmc_y,
        xlim = c(0, maxx),
        ylim = c(0, maxy),
        col = 1,
        pch = 2
    )
    # mypoints(pop_rmc_x[red_id],
    # 	pop_rmc_y[red_id],
    # 	col = 4,
    # 	pch = 2)
    mypoints(pop_rsd_x,
        pop_rsd_y,
        col = 6,
        pch = 6
    )
    my_V_ebar(pop_rmc_x, pop_rmc_y, pop_rmc_sd, col = 1)
    # my_V_ebar(pop_rmc_x[red_id], pop_rmc_y[red_id], pop_rmc_sd[red_id], col=4)
    my_V_ebar(pop_rsd_x, pop_rsd_y, pop_rsd_sd, col = 6)
    abline(0, 1)
    abline(r_rmc <- lm(obs_rmc_y ~ obs_rmc_x), lty = 2, col = 1)
    abline(r_rsd <- lm(obs_rsd_y ~ obs_rsd_x), lty = 2, col = 6)
    myaxis(1, at = seq(0, maxx, maxx / 5))
    myaxis(2, at = seq(0, maxy, maxy / 5))

    rmc_rsq <- round(summary(r_rmc)$adj.r.squared, digits = 2)
    rsd_rsq <- round(summary(r_rsd)$adj.r.squared, digits = 2)
    mymtext("Adjusted R Squared =         ", line = 0.6)
    mytitle(substitute(list(x, phantom(y)), list(x = rmc_rsq, y = rsd_rsq)),
        xlab = paste0("Ratio ", toupper(i), " - Image"),
        ylab = paste0("Ratio ", toupper(i), " - Filter Settings"),
        col.main = 1
    )
    mytitle(substitute(list(phantom(x), y), list(x = rmc_rsq, y = rsd_rsq)), col.main = 6)

    if (i == "l") {
        mylegend("topleft",
            c("Mean Cone Exc.", "Sta. Dev. Cone Exc."),
            pch = c(2, 6),
            col = c(1, 6)
        )
    }
}
dev.off()

# without the the one poorly rated image... they still both don't explain it...

# non_outliers <- which(obs.mean.dataEXTREMMEEE$img_name != "mitsuba_comp_rmc_rsd_illum_4_rg_3_by_1_ld_2_bkgd_voronoi_diff_dist.png")
non_outliers <- which(!(obs.mean.dataEXTREMMEEE$illum == "green" & obs.mean.dataEXTREMMEEE$rg == 3 & obs.mean.dataEXTREMMEEE$by == 1 & obs.mean.dataEXTREMMEEE$lighter_darker == 1 & obs.mean.dataEXTREMMEEE$bkgd == "diff_dist.png"))
obs.mean.dataEXTREMMEEE.no.weird <- obs.mean.dataEXTREMMEEE[non_outliers, ]
# pop.mean.dataEXTREMMEEE.no.weird <- aggregate(. ~ img + img_name + illum + rg + by + lighter_darker + bkgd, obs.mean.dataEXTREMMEEE.no.weird, function(x) {
#    c(MN = mean(x), SD = sd(x))
# })
pop.mean.dataEXTREMMEEE.no.weird <- aggregate(obs.mean.dataEXTREMMEEE.no.weird, list(img_name = obs.mean.dataEXTREMMEEE.no.weird$img_name, rg = obs.mean.dataEXTREMMEEE.no.weird$rg, by = obs.mean.dataEXTREMMEEE.no.weird$by, lighter_darker = obs.mean.dataEXTREMMEEE.no.weird$lighter_darker, illum = obs.mean.dataEXTREMMEEE.no.weird$illum, bkgd = obs.mean.dataEXTREMMEEE.no.weird$bkgd), function(x) {
    c(MN = mean(x), SD = sd(x))
})

pdf("../figures/rmc_rsd_comp_wout_poorly_rated_dist_EXTREMMMEEE.pdf", height = 3, width = 9)
mypar(mfrow = c(1, 3))
plotRMCandRSD(pop.mean.dataEXTREMMEEE.no.weird, obs.mean.dataEXTREMMEEE.no.weird, 1.2, 1.2, "topleft", "Adjusted R Squared =      ", 0.6)
dev.off()

# plot out data for each individual observer

for (sub_id in unique(sub_ids[[1]])) {
    obs.mean.data.sub_id <- obs.mean.data[sub_ids == sub_id, ]
    obs.mean.data.sub_id <- subset(obs.mean.data.sub_id, select = -c(sub_id))
    # pop.mean.data.sub_id <- aggregate(. ~ img + img_name + illum + rg + by + lighter_darker + bkgd, obs.mean.data.sub_id, function(x) {
    #     c(MN = mean(x), SD = sd(x))
    # })
    pop.mean.data.sub_id <- aggregate(obs.mean.data.sub_id, list(img_name = obs.mean.data.sub_id$img_name, rg = obs.mean.data.sub_id$rg, by = obs.mean.data.sub_id$by, lighter_darker = obs.mean.data.sub_id$lighter_darker, illum = obs.mean.data.sub_id$illum, bkgd = obs.mean.data.sub_id$bkgd), function(x) {
        c(MN = mean(x), SD = sd(x))
    })

    pdf(paste0("../figures/rmc_rsd_comp_", sub_id, ".pdf"), height = 3, width = 9)
    mypar(mfrow = c(1, 3))
    plotRMCandRSD(pop.mean.data.sub_id, obs.mean.data.sub_id, 3, 3, "topleft")
    dev.off()
}

# make images of the average voronoi matches that observers made with EXTREEEMMMEE filters

for (ic in 1:nrow(pop.mean.dataEXTREMMEEE)) {
    ld <- pop.mean.dataEXTREMMEEE$ld[ic, 1]
    rg <- pop.mean.dataEXTREMMEEE$rg_mix[ic, 1]
    yv <- pop.mean.dataEXTREMMEEE$by_mix[ic, 1]
# 
#     system(paste0(
#         "voronoi_filters -saveImg=true -printStats=false",
#         " -red=../base_stimuli/spectra/munsell_red_EXTREEEMMMEE.spd",
#         " -green=../base_stimuli/spectra/munsell_green_EXTREEEMMMEE.spd",
#         " -blue=../base_stimuli/spectra/munsell_blue_EXTREEEMMMEE.spd",
#         " -yellow=../base_stimuli/spectra/munsell_yellow_EXTREEEMMMEE.spd",
#         " -ld=", as.character(ld),
#         " -rg=", as.character(rg),
#         " -by=", as.character(yv),
#         " -spectraFile=../calibration/exp5.spectra.csv",
#         " -chromaFile=../calibration/exp5.chroma.csv",
#         " -gammaFile=../calibration/exp5.gamma.csv"
#     ))

    l_or_d <- ""
    if (pop.mean.dataEXTREMMEEE$lighter_darker[ic] == "darker") {
        l_or_d <- "1"
    } else {
        l_or_d <- "2"
    }

    g_or_r <- ""
    if (pop.mean.dataEXTREMMEEE$illum[ic] == "red") {
        g_or_r <- "3"
    } else {
        g_or_r <- "4"
    }
# 
#     system(paste0(
#         "mv",
#         " voronoi_filter_ld_*.png",
#         " ",
#         "../images/exp5/voronoi_matches/",
#         "voronoi_filter_illum_", g_or_r,
#         "_rg_", pop.mean.dataEXTREMMEEE$rg[ic],
#         "_by_", pop.mean.dataEXTREMMEEE$by[ic],
#         "_ld_", l_or_d,
#         "_bkgd_", pop.mean.dataEXTREMMEEE$bkgd[ic],
#         ".png"
#     ))
}

# make images of the individual voronoi matches that observers made with EXTREEEMMMEE filters

# sub_ids_EXTREMMEEE <- obs.mean.dataEXTREMMEEE["sub_id"]
# for (sub_id in unique(sub_ids_EXTREMMEEE[[1]])) {
#     obs.mean.data.sub_idEXTREMMEEE <- obs.mean.dataEXTREMMEEE[sub_ids_EXTREMMEEE == sub_id, ]
#     pop.mean.data.sub_idEXTREMMEEE <- aggregate(. ~ img + img_name + illum + rg + by + lighter_darker + bkgd, obs.mean.data.sub_idEXTREMMEEE, function(x) {
#         c(MN = mean(x), SD = sd(x))
#     })
# 
#     for (ic in 1:nrow(pop.mean.data.sub_idEXTREMMEEE)) {
#         ld <- pop.mean.data.sub_idEXTREMMEEE$ld[ic, 1]
#         rg <- pop.mean.data.sub_idEXTREMMEEE$rg_mix[ic, 1]
#         yv <- pop.mean.data.sub_idEXTREMMEEE$by_mix[ic, 1]
#
#         system(paste0(
#             "voronoi_filters -saveImg=true -printStats=false",
#             " -red=../base_stimuli/spectra/munsell_red_EXTREEEMMMEE.spd",
#             " -green=../base_stimuli/spectra/munsell_green_EXTREEEMMMEE.spd",
#             " -blue=../base_stimuli/spectra/munsell_blue_EXTREEEMMMEE.spd",
#             " -yellow=../base_stimuli/spectra/munsell_yellow_EXTREEEMMMEE.spd",
#             " -ld=", as.character(ld),
#             " -rg=", as.character(rg),
#             " -by=", as.character(yv),
#             " -spectraFile=../calibration/exp5.spectra.csv",
#             " -chromaFile=../calibration/exp5.chroma.csv",
#             " -gammaFile=../calibration/exp5.gamma.csv"
#         ))
# 
#         l_or_d <- ""
#         if (pop.mean.data.sub_idEXTREMMEEE$lighter_darker[ic] == "darker") {
#             l_or_d <- "1"
#         } else {
#             l_or_d <- "2"
#         }
# 
#         g_or_r <- ""
#         if (pop.mean.data.sub_idEXTREMMEEE$illum[ic] == "red") {
#             g_or_r <- "3"
#         } else {
#             g_or_r <- "4"
#         }
# 
#         system(paste0(
#             "mv",
#             " voronoi_filter_ld_*.png",
#             " ",
#             "../images/exp5/voronoi_matches/",
#             "voronoi_filter_subID_", sub_id,
#             "_illum_", g_or_r,
#             "_rg_", pop.mean.data.sub_idEXTREMMEEE$rg[ic],
#             "_by_", pop.mean.data.sub_idEXTREMMEEE$by[ic],
#             "_ld_", l_or_d,
#             "_bkgd_", pop.mean.data.sub_idEXTREMMEEE$bkgd[ic],
#             ".png"
#         ))
#     }
# }

# what do obs do with LAB stats?

pop.mean.dataEXTREMMEEE$color <- NA
pop.mean.dataEXTREMMEEE$pch <- NA

pop.mean.dataEXTREMMEEE$color <- 1
pop.mean.dataEXTREMMEEE$pch <- 16

mypalette()

# plotting results for whole obj in terms of mean LAB
# pdf("../figures/exp5_mean_lab.pdf", height = 3, width = 9)
# mypar(mfrow = c(1, 3))
# plotLABMean(pop.mean.dataEXTREMMEEE, obs.mean.dataEXTREMMEEE, FALSE)
# dev.off()

# plot the data of exp3 and exp5 together... maybe that becomes conclusive?

# data.exp3 <- read.csv("../data/exp3_obs_data.csv")
data.exp3 <- data[data$exp_name == "filter" & data$mask_name == "obj_mask", ]

cols_to_stack <- c("rmc_l_glaven", "rmc_m_glaven", "rmc_s_glaven", "rsd_l_glaven", "rsd_m_glaven", "rsd_s_glaven", "rmc_l_match", "rmc_m_match", "rmc_s_match", "rsd_l_filter", "rsd_m_filter", "rsd_s_filter")

data_exp3_exp5 <- rbind(data.exp3[, cols_to_stack], data.from.exp[, cols_to_stack])

pdf("../figures/exp3_and_exp5_results_together.pdf", height = 3, width = 9)
maxx <- 4
maxy <- 4
mypar(mfrow = c(1, 3))
for (i in c("l", "m", "s")) {
    pop_rmc_x <- data_exp3_exp5[[paste0("rmc_", i, "_glaven")]]
    pop_rmc_y <- data_exp3_exp5[[paste0("rmc_", i, "_match")]]
    pop_rsd_x <- data_exp3_exp5[[paste0("rsd_", i, "_glaven")]]
    pop_rsd_y <- data_exp3_exp5[[paste0("rsd_", i, "_filter")]]

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
    abline(0, 1)
    abline(r_rmc <- lm(pop_rmc_y ~ pop_rmc_x), lty = 2, col = 1)
    abline(r_rsd <- lm(pop_rsd_y ~ pop_rsd_x), lty = 2, col = 6)
    myaxis(1, at = seq(0, maxx, maxx / 5))
    myaxis(2, at = seq(0, maxy, maxy / 5))

    rmc_rsq <- round(summary(r_rmc)$adj.r.squared, digits = 2)
    rsd_rsq <- round(summary(r_rsd)$adj.r.squared, digits = 2)
    mymtext("Adjusted R Squared = ", line = 0.5)
    mytitle(substitute(list(x, phantom(y)), list(x = rmc_rsq, y = rsd_rsq)),
        xlab = paste0("Ratio ", toupper(i), " - Image"),
        ylab = paste0("Ratio ", toupper(i), " - Filter Settings"),
        col.main = 1
    )
    mytitle(substitute(list(phantom(x), y), list(x = rmc_rsq, y = rsd_rsq)), col.main = 6)

    if (i == "l") {
        mylegend("topleft",
            c("Mean Cone Exc.", "Sta. Dev. Cone Exc."),
            pch = c(2, 6),
            col = c(1, 6)
        )
    }
}
dev.off()

# print out RMC and RSD for each tested image, to put in plot

for (x in 1:nrow(pop.mean.dataEXTREMMEEE)) {
    img_name <- pop.mean.dataEXTREMMEEE$img_name[x]

    rmc_l <- pop.mean.dataEXTREMMEEE$rmc_l_glaven[x]
    rmc_m <- pop.mean.dataEXTREMMEEE$rmc_m_glaven[x]
    rmc_s <- pop.mean.dataEXTREMMEEE$rmc_s_glaven[x]

    rsd_l <- pop.mean.dataEXTREMMEEE$rsd_l_glaven[x]
    rsd_m <- pop.mean.dataEXTREMMEEE$rsd_m_glaven[x]
    rsd_s <- pop.mean.dataEXTREMMEEE$rsd_s_glaven[x]

    print(paste0(
        img_name,
        " ---> RMC = (", rmc_l, ", ", rmc_m, ", ", rmc_s,
        "), RSD = (", rsd_l, ", ", rsd_m, ", ", rsd_s, ")"
    ))
}

for (x in 1:nrow(pop.mean.dataEXTREMMEEE)) {
    qualityc <- pop.mean.dataEXTREMMEEE$qualityc[x, ]

    img_name <- pop.mean.dataEXTREMMEEE$img_name[x]
    print(paste0(img_name, " ---> QualityC = ", round(qualityc[1], digit = 3), " +/- ", round(qualityc[2], digit = 3)))
}