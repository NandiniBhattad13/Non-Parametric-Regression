f <- function(x)
{
  return((exp(x))^2)
}
convergence <- function(N)
{
  errors <- rnorm(N)
  hn <- N^(-1/3)
  ctr <- 0
  ctr1 <- 0
  fun <- function(x)
  {
    ctr <- ctr+1
    return((x^ctr)^2)
  }
  y_i <- numeric(N)
  for(i in 1:N)
  {
    r <- integrate(fun,0,1)
    y_i[i] <- sqrt(r$value) + errors[i]
  }
  
  func <- function(x)
  {
    ctr1 <- ctr1 + 1
    return(((exp(x))-(x^ctr1))^2)
  }
  num <- numeric(N)
  denom <- numeric(N)
  for(i in 1:N)
  {
   result <- integrate(func,0,1)
   num[i] <- dnorm(sqrt(result$value)/hn)*y_i[i]
   denom[i] <- dnorm(sqrt(result$value)/hn)
  }
  mn_hat <- sum(num)/sum(denom)
  return(mn_hat)
}
real <- integrate(f,0,1)
real <- sqrt(real$value)
########################################################

N_values <- seq(10, 10000, by = 10) 
# Calculate mn_hat for each N value
mn_hat_values <- sapply(N_values, convergence)

# Plot the convergence behavior
plot(N_values, mn_hat_values, type = "l", 
     xlab = "Number of Observations (N)", ylab = "Estimated Mean (mn_hat)",
     main = "Convergence of mn_hat as N Increases", ylim = c(-2,2))
abline(h = real, col = "red", lty = 2)  # True value line
