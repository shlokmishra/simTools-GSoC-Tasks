library(MASS)

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
mu <- c(1, 40, 13)
sigma <- matrix(c(2, 0.5, 0.5, 0.5, 3, 0.2, 0.5, 0.2, 1), nrow = 3)
result <- metropolisHastings(n, mu, sigma)
accept.rate <- result$accept.rate
chain <- result$chain

plot(mvrnorm(n,mu,sigma,tol = 1e-6, empirical = FALSE, EISPACK = FALSE ), pch = 16)
points(chain, col = "blue")
legend("bottomleft",legend =c("Target", "Sampled"), pch = c(16,1), col = c("black", "blue"))
