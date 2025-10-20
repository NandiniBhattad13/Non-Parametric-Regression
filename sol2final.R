library(polynom)
library(orthopolynom)
m_f <- function(x) 
{
  return(sin(x)^2*sin(x)^2)
}
convergence <- function(N, x) 
{
  
  onb <- legendre.polynomials(N, normalized = TRUE)
  Mj_f <- numeric(N)
  errors <- rnorm(N)
  e <- numeric(N)
  y <- numeric(N)
  hn <- N^(-1/5)
  
  for (j in 1:N) 
  {
    integrand <- function(x) 
    {
      m_f(x) * predict(onb[[j]], x)
    }
    error <- function(x) {
      errors[j] * predict(onb[[j]], x)
    }
    Mj_f[j] <- integrate(integrand, lower = 0, upper = 1)$value
    e[j] <- integrate(error, lower = 0, upper = 1)$value
    y[j] <- Mj_f[j] + e[j]
  }
  
  Mhat_j <- numeric(N)
  
  for (j in 1:N) {
    num <- numeric(N)
    denom <- numeric(N)
    for (i in 1:N) {
      t <- function(x) {
        y[i] * predict(onb[[j]], x^2)
      }
      foo <- integrate(t, lower = 0, upper = 1)$value
      num[i] <- dnorm((foo - Mj_f[j]) / hn) * foo
      denom[i] <- dnorm((foo - Mj_f[j]) / hn)
    }
    Mhat_j[j] <- sum(num) / sum(denom)
  }
  
  sum_result <- 0
  
  for (j in 1:N)
  {
    onb[[j]] <- Mhat_j[j] * onb[[j]]
    sum_result <- sum_result + onb[[j]]
  }
  return(as.function(sum_result)(x))
}

x <- seq(0, 1, by = 0.01)
y <- convergence(10, x)  
plot(x, m_f(x), type = "l", col = "red", ylim = c(-2, 2), ylab = "y")
lines(x, y, col = "blue")
legend("topright", legend = c("Original function", "Approximation"), col = c("red", "blue"), lty = 1)
