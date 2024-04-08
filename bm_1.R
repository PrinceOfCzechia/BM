rm(list=ls())

library(MASS)
library(pracma)

N = 50
y = 59.6
s = 11.1


###
# 1)
###

# mean
c( y-qt(0.975,N)*s/sqrt(N), y+qt(0.975,N)*s/sqrt(N) )

# std deviation
c( sqrt((N-1)*s^2/qchisq(0.975,N-1)), sqrt((N-1)*s^2/qchisq(0.025,N-1)) )

# inverse variance
c( qchisq(0.025,N-1)/(49*11.6^2), qchisq(0.975,N-1)/(49*11.6^2) )


###
# 2)
###
density_mu = function(N, g, y, s, h, mu) {
  numerator = gamma(N/2 + g)
  denominator = (0.5 * N * s^2 + N * (y - mu)^2 + h)^(N/2 + g)
  return = ( numerator / denominator )
}

mu_values = seq(50, 70, length.out = 100)

g = 0.001
h = 0.001
expression_values = sapply(mu_values, function(mu) density_mu(N, g, y, s, h, mu))
expression_values = expression_values * 1e65

df = data.frame(mu = mu_values, expression = expression_values)
ggplot(df, aes(x = mu, y = expression)) +
  geom_line() +
  labs(x = "Mean (mu)", y = "Expression Value", title = "Expression vs. Mean")

g = 1
h = 0.005
expression_values = sapply(mu_values, function(mu) density_mu(N, g, y, s, h, mu))
expression_values = expression_values * 1e65

df = data.frame(mu = mu_values, expression = expression_values)
ggplot(df, aes(x = mu, y = expression)) +
  geom_line() +
  labs(x = "Mean (mu)", y = "Expression Value", title = "Expression vs. Mean")

g = 0
h = 0
expression_values = sapply(mu_values, function(mu) density_mu(N, g, y, s, h, mu))
expression_values = expression_values * 1e65

df = data.frame(mu = mu_values, expression = expression_values)
ggplot(df, aes(x = mu, y = expression)) +
  geom_line() +
  labs(x = "Mean (mu)", y = "Expression Value", title = "Expression vs. Mean")

# Compute integral from -infinity to infinity
i_full = integrate(density_mu, lower = -100, upper = 100,
                           N = 50, g = 1, y = 59.6, s = 11.1, h = 0.005)

# Compute integral from 55 to infinity
i_55 = integrate(denisty_mu, lower = 55, upper = 100,
                             N = 50, g = 1, y = 59.6, s = 11.1, h = 0.005)

# Print results
print(i_full)
print(i_55)
(i_55$value/i_full$value)


###
# 3)
###
density_tau = function(N, g, s, h, tau)
{
  return = ( tau^(N/2 + g - 1) * exp(-0.5 * tau * (N * s^2 + h)) * sqrt(2 * pi / (N * tau)) )
}

g = 0.001
h = 0.001
tau_values = seq(0.001, 0.020, length.out = 100)
density_values = density_tau(N, g, s, h, tau_values) * 1e60

plot(tau_values, density_values, type = "l", xlab = "Tau", ylab = "density_tau",
     main = "density of tau")

density_sigma = function(N, g, s, h, sigma) # using trafo of density_tau
{
  result = density_tau(N, g, s, h, (sigma^2))*2*sigma
  return(result)
}

sigma_values = seq(0.05, 0.20, length.out = 100)
density_values = density_sigma(N, g, s, h, sigma_values) * 0.5e61

plot(sigma_values, density_values, type = "l", xlab = "sigma", ylab = "density_sigma",
     main = "density of sigma")
