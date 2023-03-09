# Define the simulation function for AR(1) model
ar1_sim <- function(n, phi, sigma) {
  x <- numeric(length = n)
  x[1] <- rnorm(1)
  for (i in 2:n) {
    x[i] <- phi * x[i - 1] + rnorm(1, mean = 0, sd = sigma)
  }
  return(x)
}

# Example usage
set.seed(1)
n <- 1e3
phi <- 0.8
sigma <- 1
ar1_series <- ar1_sim(n, phi, sigma)
ts.plot(ar1_series)

acf(ar1_series)
acf(ar1_series,plot=FALSE)[0:10]
pacf(ar1_series)
pacf(ar1_series,plot=FALSE)[0:5]
