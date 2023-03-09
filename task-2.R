library(MASS)

# Define the log likelihood function for a multivariate normal distribution
logLikelihood <- function(x, mu, sigma) {
  p <- length(mu)
  -0.5 * (p * log(2 * pi) + log(det(sigma)) + t(x - mu) %*% solve(sigma) %*% (x - mu))
}

# Define the Metropolis-Hastings algorithm function
metropolisHastings <- function(n, mu, sigma) {
  p <- length(mu)
  # Initialize the chain
  x <- rep(0, p)
  # Initialize the acceptance counter
  accept <- 0
  # Initialize the chain vector
  chain <- matrix(0, n, p)
  for (i in 1:n) {
    # Propose a new point
    x.prop <- rnorm(p, mean = x, sd = 1)
    # Calculate the acceptance ratio
    log.ratio <- logLikelihood(x.prop, mu, sigma) - logLikelihood(x, mu, sigma)
    # Accept or reject the proposed point
    if (log(runif(1)) < log.ratio) {
      x <- x.prop
      accept <- accept + 1
    }
    # Save the current point to the chain
    chain[i, ] <- x
  }
  # Return the chain and the acceptance rate
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

plot(mvrnorm(n,mu,sigma,tol = 1e-6, empirical = FALSE, EISPACK = FALSE ))
points(chain, col = "blue")
legend()