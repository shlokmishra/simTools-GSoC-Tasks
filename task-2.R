library(MASS)
library(mvtnorm)

logLikelihood <- function(x, mu, sigma) {
  p <- length(mu)
  -0.5 * (p * log(2 * pi) + log(det(sigma)) + t(x - mu) %*% solve(sigma) %*% (x - mu))
}

metropolisHastings <- function(n, mu, sigma) {
  p <- length(mu)
  x <- rep(0, p)
  accept <- 0
  chain <- matrix(0, n, p)
  for (i in 1:n) {
    x.prop <- rnorm(p, mean = x, sd = 1)
    log.ratio <- logLikelihood(x.prop, mu, sigma) - logLikelihood(x, mu, sigma)
    if (log(runif(1)) < log.ratio) {
      x <- x.prop
      accept <- accept + 1
    }
    chain[i, ] <- x
  }
  list(chain = chain, accept.rate = accept/n)
}

# Example usage
set.seed(1)
n <- 1e4
mu <- c(1, 2)
sigma <- matrix(c(2, 0.5, 0.5, 3), nrow = 2)
result <- metropolisHastings(n, mu, sigma)
accept.rate <- result$accept.rate
chain <- result$chain

plot(mvrnorm(n,mu,sigma,tol = 1e-6, empirical = FALSE, EISPACK = FALSE ), pch = 4, ylab = "", xlab = "")
points(chain, col = "blue")
legend("topright",legend =c("Target", "Sampled"), pch = c(4,1), col = c("black", "blue"))
title(main = "Samples from bivariate Normal distribution", xlab = " ", ylab = " ")



x <- seq(-6, 8, length.out = 100)
y <- seq(-6, 8, length.out = 100)
xy <- expand.grid(x = x, y = y)
z <- dmvnorm(xy, mean = mu, sigma = sigma)

zmat <- matrix(z, nrow = length(x), ncol = length(y), byrow = TRUE)

persp(x, y, zmat, theta = 30, phi = 30, col = "lightblue", border = NA,
      xlab = "X", ylab = "Y", zlab = "Density")
