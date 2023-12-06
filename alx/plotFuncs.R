plotLABMean <- function(pop.data, obs.data, neutral) {
    for (i in c("L", "a", "b")) {
        pop_x <- pop.data[[paste0("gw", i, "_glaven.MN")]]
        pop_y <- pop.data[[paste0("gw", i, "_match.MN")]]
        pop_sd <- pop.data[[paste0("gw", i, "_match.SD")]]

        obs_x <- obs.data[[paste0("gw", i, "_glaven.MN")]]
        obs_y <- obs.data[[paste0("gw", i, "_match.MN")]]

        if (i == "L") {
            minlim <- 0
            maxlim <- 100
            i <- toupper(i)
        } else {
            minlim <- -100
            maxlim <- 100
        }

        myplot(pop_x,
            pop_y,
            xlim = c(minlim, maxlim),
            ylim = c(minlim, maxlim),
            col = pop.data$color,
            pch = pop.data$pch
        )
        my_V_ebar(pop_x, pop_y, pop_sd, col = pop.data$color)
        abline(0, 1)
        abline(rg <- lm(obs_y ~ obs_x), lty = 2, col = "darkslategray")

        rsq <- round(summary(rg)$adj.r.squared, digits = 2)
        mytitle(paste0("Adjusted R Squared = ", rsq),
            xlab = paste0("Mean Image ", i, "*"),
            ylab = paste0("Mean Observer Match ", i, "*")
        )

        if (i == "L") {
            myaxis(1, at = seq(0, 100, 20))
            myaxis(2, at = seq(0, 100, 20))
            if (neutral) {
                mylegend("bottomright", c("Neutral Illuminant"), pch = c(19, 6))
            } else {
                mylegend("bottomright", c("Blue Illuminant", "Yellow Illuminant"), pch = c(19, 6))
            }
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
}

plotRMCandRSD <- function(pop.data, obs.data, maxx, maxy, legendloc, adjtext = "Adjusted R Squared = ", lineval = 0.5) {
    for (i in c("l", "m", "s")) {
        pop_rmc_x <- pop.data[[paste0("rmc_", i, "_glaven")]][, 1]
        pop_rmc_y <- pop.data[[paste0("rmc_", i, "_match")]][, 1]
        pop_rsd_x <- pop.data[[paste0("rsd_", i, "_glaven")]][, 1]
        pop_rsd_y <- pop.data[[paste0("rsd_", i, "_filter")]][, 1]
        pop_rmc_sd <- pop.data[[paste0("rmc_", i, "_match")]][, 2]
        pop_rsd_sd <- pop.data[[paste0("rsd_", i, "_filter")]][, 2]

        obs_rmc_x <- obs.data[[paste0("rmc_", i, "_glaven")]]
        obs_rmc_y <- obs.data[[paste0("rmc_", i, "_match")]]
        obs_rsd_x <- obs.data[[paste0("rsd_", i, "_glaven")]]
        obs_rsd_y <- obs.data[[paste0("rsd_", i, "_filter")]]

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
        my_V_ebar(pop_rmc_x, pop_rmc_y, pop_rmc_sd, col = 1)
        my_V_ebar(pop_rsd_x, pop_rsd_y, pop_rsd_sd, col = 6)
        abline(0, 1)
        abline(r_rmc <- lm(obs_rmc_y ~ obs_rmc_x), lty = 2, col = 1)
        abline(r_rsd <- lm(obs_rsd_y ~ obs_rsd_x), lty = 2, col = 6)
        myaxis(1, at = seq(0, maxx, maxx / 5))
        myaxis(2, at = seq(0, maxy, maxy / 5))

        rmc_rsq <- round(summary(r_rmc)$adj.r.squared, digits = 2)
        rsd_rsq <- round(summary(r_rsd)$adj.r.squared, digits = 2)
        mymtext(adjtext, line = lineval)
        mytitle(substitute(list(x, phantom(y)), list(x = rmc_rsq, y = rsd_rsq)),
            xlab = paste0("Ratio ", toupper(i), " - Image"),
            ylab = paste0("Ratio ", toupper(i), " - Filter Settings"),
            col.main = 1
        )
        mytitle(substitute(list(phantom(x), y), list(x = rmc_rsq, y = rsd_rsq)), col.main = 6)

        if (i == "l") {
            mylegend(legendloc,
                c("Mean Cone Exc.", "Sta. Dev. Cone Exc."),
                pch = c(2, 6),
                col = c(1, 6)
            )
        }
    }
}

plotLMSMean <- function(pop.data, obs.data, meanstatname, neutral) {
    for (i in c("l", "m", "s")) {
        pop_x <- pop.data[[paste0("mean_cone_", i)]][, 1]
        pop_y <- pop.data[[paste0(meanstatname, i)]][, 1]
        pop_sd <- pop.data[[paste0(meanstatname, i)]][, 2]

        obs_x <- obs.data[[paste0("mean_cone_", i)]]
        obs_y <- obs.data[[paste0(meanstatname, i)]]

        minlim <- 0
        maxlim <- max(c(obs_x, obs_y)) + 0.1
        i <- toupper(i)

        myplot(pop_x,
            pop_y,
            xlim = c(minlim, maxlim),
            ylim = c(minlim, maxlim),
            col = pop.data$color,
            pch = pop.data$pch
        )
        my_V_ebar(pop_x, pop_y, pop_sd, col = pop.data$color)
        abline(0, 1)
        abline(rg <- lm(obs_y ~ obs_x), lty = 2, col = "darkslategray")

        rsq <- round(summary(rg)$adj.r.squared, digits = 2)
        mytitle(paste0("Adjusted R Squared = ", rsq),
            xlab = paste0("Mean Image ", i),
            ylab = paste0("Mean Observer Match ", i)
        )

        if (i == "L") {
            myaxis(1, at = seq(0, 1, 5))
            myaxis(2, at = seq(0, 1, 5))
            if (neutral) {
                mylegend("bottomright", c("Neutral Illuminant"), pch = c(19, 6))
            } else {
                mylegend("bottomright", c("Blue Illuminant", "Yellow Illuminant"), pch = c(19, 6))
            }
        } else if (i == "a") {
            myaxis(1, at = seq(0, 1, 5))
            myaxis(2, at = seq(0, 1, 5))
        } else if (i == "b") {
            myaxis(1, at = seq(0, 1, 5))
            myaxis(2, at = seq(0, 1, 5))
        }
    }
}

plotDKLMean <- function(pop.data, obs.data, meanstatname, neutral) {
    for (i in c("ld", "rg", "yv")) {
        pop_x <- pop.data[[paste0("gw", i)]][, 1]
        pop_y <- pop.data[[paste0(meanstatname, i)]][, 1]
        pop_sd <- pop.data[[paste0(meanstatname, i)]][, 2]

        obs_x <- obs.data[[paste0("gw", i)]]
        obs_y <- obs.data[[paste0(meanstatname, i)]]

        minlim <- -1
        maxlim <- 1
        i <- toupper(i)

        myplot(pop_x,
            pop_y,
            xlim = c(minlim, maxlim),
            ylim = c(minlim, maxlim),
            col = pop.data$color,
            pch = pop.data$pch
        )
        my_V_ebar(pop_x, pop_y, pop_sd, col = pop.data$color)
        abline(0, 1)
        abline(rg <- lm(obs_y ~ obs_x), lty = 2, col = "darkslategray")

        rsq <- round(summary(rg)$adj.r.squared, digits = 2)
        mytitle(paste0("Adjusted R Squared = ", rsq),
            xlab = paste0("Mean Image ", i),
            ylab = paste0("Mean Observer Match ", i)
        )

        if (i == "L") {
            myaxis(1, at = seq(-1, 1, 8))
            myaxis(2, at = seq(-1, 1, 8))
            if (neutral) {
                mylegend("bottomright", c("Neutral Illuminant"), pch = c(19, 6))
            } else {
                mylegend("bottomright", c("Blue Illuminant", "Yellow Illuminant"), pch = c(19, 6))
            }
        } else if (i == "a") {
            myaxis(1, at = seq(-1, 1, 8))
            myaxis(2, at = seq(-1, 1, 8))
        } else if (i == "b") {
            myaxis(1, at = seq(-1, 1, 8))
            myaxis(2, at = seq(-1, 1, 8))
        }
    }
}

plotConfEllipse <- function(covm, avg) {
	eig <- eigen(covm)
	
	# Get the index of the largest eigenvector
	largest_eigenvec_ind_c <- which.max(eig$values)
	largest_eigenvec <- eig$vectors[, largest_eigenvec_ind_c]

	# Get the largest eigenvalue
	largest_eigenval <- max(eig$values)

	# Get the smallest eigenvector and eigenvalue
	if (largest_eigenvec_ind_c == 1) {
	    smallest_eigenval <- eig$values[[2]]
	    smallest_eigenvec <- eig$vectors[, 2]
	} else {
	    smallest_eigenval <- eig$values[[1]]
	    smallest_eigenvec <- eig$vectors[, 1]
	}

	# Calculate the angle between the x-axis and the largest eigenvector
	angle <- atan2(largest_eigenvec[2], largest_eigenvec[1])

	# This angle is between -pi and pi.
	# Let's shift it such that the angle is between 0 and 2pi
	if (angle < 0) {
	    angle <- angle + 2*pi
	}

	# Get the 95% confidence interval error ellipse
	chisquare_val <- 2.4477
	theta_grid <- seq(0, 2*pi, length.out = 100)
	X0 <- avg[1]
	Y0 <- avg[2]
	a <- chisquare_val*sqrt(largest_eigenval)
	b <- chisquare_val*sqrt(smallest_eigenval)

	# the ellipse in x and y coordinates
	ellipse_x_r <- a*cos(theta_grid)
	ellipse_y_r <- b*sin(theta_grid)

	# Define a rotation matrix
	R <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), nrow = 2, ncol = 2)

	# let's rotate the ellipse to some angle phi
	r_ellipse <- matrix(c(ellipse_x_r, ellipse_y_r), nrow = length(ellipse_x_r), ncol = 2) %*% R
}