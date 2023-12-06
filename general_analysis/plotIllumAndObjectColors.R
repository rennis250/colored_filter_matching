colormatchPlanck <- function(spect, wlns) {
	cmfData <- read.csv('C:\\Users\\me\\Documents\\source_code\\MATLAB\\matlabHyper\\ciexyz31_1.csv')

	wavelength.cmf <- cmfData[ , 1]

	x.bar <- cmfData[ , 2]
	y.bar <- cmfData[ , 3]
	z.bar <- cmfData[ , 4]

	x.cmf.fun <- splinefun(x = wavelength.cmf, y = x.bar)
	x.cmf <- x.cmf.fun(wlns)
	x.cmf[x.cmf < 0.0] <- 0.0
	y.cmf.fun <- splinefun(x = wavelength.cmf, y = y.bar)
	y.cmf <- y.cmf.fun(wlns)
	y.cmf[y.cmf < 0.0] <- 0.0
	z.cmf.fun <- splinefun(x = wavelength.cmf, y = z.bar)
	z.cmf <- z.cmf.fun(wlns)
	z.cmf[z.cmf < 0.0] <- 0.0

	XYZ <- c(0, 0, 0)

	stepSize <- c(diff(wlns), 0)
	x.corrCMF <- 683 * x.cmf * stepSize
	y.corrCMF <- 683 * y.cmf * stepSize
	z.corrCMF <- 683 * z.cmf * stepSize

	XYZ[1] <- spect %*% x.corrCMF
	XYZ[2] <- spect %*% y.corrCMF
	XYZ[3] <- spect %*% z.corrCMF

	# monxyY <- read.csv('OLEDxyY.csv')
	# RGB <- XYZ2RGB(XYZ, monxyY)

	return(XYZ)
}

XYZtoxyY <- function(XYZ) {
	s <- sum(XYZ)
	xyY <- c(XYZ[1]/s, XYZ[2]/s, XYZ[2])
	return(xyY)
}

# load distributions for illum and object colors
lightredrefl <- read.table('../base_stimuli/spectra/munsell_red_lightest.spd', sep=' ')
lightbluerefl <- read.table('../base_stimuli/spectra/munsell_blue_lightest.spd', sep=' ')
lightgreenrefl <- read.table('../base_stimuli/spectra/munsell_green_lightest.spd', sep=' ')
lightyellowrefl <- read.table('../base_stimuli/spectra/munsell_yellow_lightest.spd', sep=' ')

wlns <- lightredrefl[ , 1]

lightredxyy <- XYZtoxyY(colormatchPlanck(lightredrefl[ , 2], wlns))
lightgreenxyy <- XYZtoxyY(colormatchPlanck(lightgreenrefl[ , 2], wlns))
lightbluexyy <- XYZtoxyY(colormatchPlanck(lightbluerefl[ , 2], wlns))
lightyellowxyy <- XYZtoxyY(colormatchPlanck(lightyellowrefl[ , 2], wlns))

bluespectra <- read.table('../base_stimuli/spectra/de_brighter_illum_blue.spd', sep=' ')
yellowspectra <- read.table('../base_stimuli/spectra/de_brighter_illum_yellow.spd', sep=' ')

wlns <- bluespectra[ , 1]

bluesxyy <- XYZtoxyY(colormatchPlanck(bluespectra[ , 2], wlns))
yellowsxyy <- XYZtoxyY(colormatchPlanck(yellowspectra[ , 2], wlns))

# load distributions for filter matching

filterredrefl <- read.table('../base_stimuli/spectra/munsell_red_better.spd', sep=' ')
filterbluerefl <- read.table('../base_stimuli/spectra/munsell_blue_better.spd', sep=' ')
filtergreenrefl <- read.table('../base_stimuli/spectra/munsell_green_better.spd', sep=' ')
filteryellowrefl <- read.table('../base_stimuli/spectra/munsell_yellow_better.spd', sep=' ')

wlns <- lightredrefl[ , 1]

filterredxyy <- XYZtoxyY(colormatchPlanck(filterredrefl[ , 2], wlns))
filtergreenxyy <- XYZtoxyY(colormatchPlanck(filtergreenrefl[ , 2], wlns))
filterbluexyy <- XYZtoxyY(colormatchPlanck(filterbluerefl[ , 2], wlns))
filteryellowxyy <- XYZtoxyY(colormatchPlanck(filteryellowrefl[ , 2], wlns))

# grab pre-computed xyY coords of CIE lasso and planckian locus for plotting
lassocs <- read.csv('C:\\Users\\me\\Documents\\source_code\\MATLAB\\xyYcoords.csv')
planck <- read.csv('C:\\Users\\me\\Documents\\source_code\\MATLAB\\planck.csv')

mypalette()

pdf("../figures/refl_illum_coords.pdf", width = 9, height = 4)
mypar(mfrow=c(1, 2), pty="s", cex.axis=1, cex = 1)
myplot(NA, NA, xlim=c(0, 0.8), ylim=c(0, 0.9))

mylines(lassocs[ , 1], lassocs[ , 2])
mylines(planck[ , 1], planck[ , 2])

mypoints(lightredxyy[1], lightredxyy[2], col=4, pch=16, cex=1.1)
mypoints(lightgreenxyy[1], lightgreenxyy[2], col=5, pch=16, cex=1.1)
mypoints(lightbluexyy[1], lightbluexyy[2], col=6, pch=16, cex=1.1)
mypoints(lightyellowxyy[1], lightyellowxyy[2], col=7, pch=16, cex=1.1)

mypoints(filterredxyy[1], filterredxyy[2], col=4, pch=17, cex=1.1)
mypoints(filtergreenxyy[1], filtergreenxyy[2], col=5, pch=17, cex=1.1)
mypoints(filterbluexyy[1], filterbluexyy[2], col=6, pch=17, cex=1.1)
mypoints(filteryellowxyy[1], filteryellowxyy[2], col=7, pch=17, cex=1.1)

mypoints(bluesxyy[1], bluesxyy[2], col=6, pch=8, cex=0.8)
mypoints(yellowsxyy[1], yellowsxyy[2], col=1, pch=8, cex=0.8)

myaxis(1, at = seq(0, 0.8, 0.8/5))
myaxis(2, at = round(seq(0, 0.89, 0.89/5), digits=2))
mytitle(xlab = "x", ylab = "y")

myplot(NA, NA, xlim=c(380, 830), ylim=c(0, max(bluespectra[, 2])+1))
mylines(bluespectra[, 1], bluespectra[, 2], col=6)
mylines(yellowspectra[, 1], yellowspectra[, 2], col=1)

myaxis(1, at = round(seq(380, 830, (830 - 380)/5)))
myaxis(2, at = round(seq(min(bluespectra[, 2]), max(bluespectra[, 2]), (max(bluespectra[, 2]) - min(bluespectra[, 2]))/5)))
mytitle(xlab = "Wavelength (nm)", ylab = expression(paste("Radiance ", ~~ (frac(W, str %*% m^{2})))))
dev.off()
