d <- data.frame(a = obs_data[["mean_a"]], b = obs_data[["mean_b"]])

covm <- cov(data.frame(a = obs_data[["mean_a"]], b = obs_data[["mean_b"]]))

# Get the correlation matrix
P <- cov2cor(covm)

mu <- c(mean(d[["a"]]), mean(d[["b"]]))

# Plot the data
dataEllipse(d[["a"]], d[["b"]], xlim = c(-100, 100), ylim = c(-100, 100), levels = c(0.95), center.pch = FALSE)
points(d[["a"]], d[["b"]])

# Plot the ellipse
lines(ellipse::ellipse(covm, centre = mu), col='red')

evals <- eigen(covm)$values
evecs <- eigen(covm)$vectors

# Angles of a circle
a <- seq(0, 2*pi, len=100)

# Get critical value
c2 <- qchisq(0.95, 2)
c <- sqrt(c2)

# Get the distances
xT <- c * sqrt(evals[1]) * cos(a)
yT <- c * sqrt(evals[2]) * sin(a)

M <- cbind(xT, yT)

# Covert the coordinates
transM <- evecs %*% t(M)
transM <- t(transM)

lines(transM + mu)

set.seed(17)

# Set the covariance matrix
sigma2 <- covm
# Set the means
mu <- c(5,5)
# Get the correlation matrix
P <- cov2cor(sigma2)
# Generate the data
p <- rmvnorm(n=50, mean=mu, sigma=sqrt(sigma2))
# Plot the data
plot(p)
dataEllipse(p[, 1], p[, 2], xlim = c(-20, 20), ylim = c(-20, 20), levels = c(0.95), center.pch = FALSE)
points(p[, 1], p[, 2])
# Plot the ellipse
lines( ellipse::ellipse(sqrt(covm), centre = c(5,5)) , col='red')
evals <- eigen(P)$values
evecs <- eigen(P)$vectors
# Angles of a circle
a <- seq(0, 2*pi, len=100)
# Get critical value
c2 <- qchisq(0.95, 2)
c <- sqrt(c2)
# Get the distances
xT <- c * sqrt(evals[1]) * cos(a)
yT <- c * sqrt(evals[2]) * sin(a)
M <- cbind(xT, yT)
# Covert the coordinates
transM <- evecs %*% t(M)
transM <- t(transM)
lines(transM + mu)