rm(list=ls())

library(MASS)
library(pracma)
library(coda)

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
  return = ( 1e63 * numerator / denominator )
}

mu_values = seq(50, 70, length.out = 100)

g = 0.001
h = 0.001
expression_values = sapply(mu_values, function(mu) density_mu(N, g, y, s, h, mu))
expression_values = expression_values * 1e63

df = data.frame(mu = mu_values, expression = expression_values)
ggplot(df, aes(x = mu, y = expression)) +
  geom_line() +
  labs(x = "Mean (mu)", y = "Expression Value", title = "Expression vs. Mean")

g = 1
h = 0.005
expression_values = sapply(mu_values, function(mu) density_mu(N, g, y, s, h, mu))
expression_values = expression_values

df = data.frame(mu = mu_values, expression = expression_values)
ggplot(df, aes(x = mu, y = expression)) +
  geom_line() +
  labs(x = "Mean (mu)", y = "Expression Value", title = "Expression vs. Mean")

g = 0
h = 0
expression_values = sapply(mu_values, function(mu) density_mu(N, g, y, s, h, mu))
expression_values = expression_values

df = data.frame(mu = mu_values, expression = expression_values)
ggplot(df, aes(x = mu, y = expression)) +
  geom_line() +
  labs(x = "Mean (mu)", y = "Expression Value", title = "Expression vs. Mean")

# Compute integral from -infinity to infinity
i_full = integrate(density_mu, lower = -100, upper = 100,
                           N = 50, g = 1, y = 59.6, s = 11.1, h = 0.005)

# Compute integral from 55 to infinity
i_55 = integrate(density_mu, lower = 55, upper = 100,
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

g = 0
h = 0
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


###
# 6)
###

# ET for mu
CDF_mu = function(mu)
{ 
  CDF = integrate(function(mu) 1e63 * density_mu(N, g, y, s, h, mu), lower= - 100, upper = mu)$value
  norm = integrate(function(mu) 1e63 *density_mu(N, g, y, s, h, mu), lower= - 100, upper = 100)$value
  return = CDF/norm
}
qL_mu = uniroot(function(mu) {CDF_mu(mu) - 0.025}, interval = c(55, 60))$root
qU_mu = uniroot(function(mu) {CDF_mu(mu) - 0.975}, interval = c(60, 65))$root
ET_mu = c(qL_mu,qU_mu)

# ET for tau
CDF_tau = function(tau)
{ 
  CDF = integrate(function(tau) 1e60 * density_tau(N, g, s, h, tau), lower= 0.001, upper = tau)$value
  norm = integrate(function(tau) 1e60 * density_tau(N, g, s, h, tau), lower= 0.001, upper = 0.012)$value
  return = CDF/norm
}
qL_tau = uniroot(function(tau) {CDF_tau(tau) - 0.025}, interval = c(0.001, 0.007))$root
qU_tau = uniroot(function(tau) {CDF_tau(tau) - 0.975}, interval = c(0.009, 0.015))$root
ET_tau = c(qL_tau,qU_tau)

# ET for sigma
CDF_sigma = function(sigma)
{ 
  CDF = integrate(function(sigma) 0.5e61 * density_sigma(N, g, s, h, sigma), lower= 0.05, upper = sigma)$value
  norm = integrate(function(sigma) 0.5e61 * density_sigma(N, g, s, h, sigma), lower= 0.05, upper = 0.15)$value
  return = CDF/norm
}
qL_sigma = uniroot(function(sigma) {CDF_sigma(sigma) - 0.025}, interval = c(0.05, 0.09))$root
qU_sigma = uniroot(function(sigma) {CDF_sigma(sigma) - 0.975}, interval = c(0.09, 0.12))$root
ET_sigma = c(qL_sigma,qU_sigma)


###
# 7)
###

# HPD for mu
C = 0.95
mu_values = seq(50, 70, by = 0.001)
density_values = density_mu(N, g, y, s, h, mu_values)
sorted = sort(density_values, decreasing = TRUE)
cum_sum = cumsum(sorted)
index = which(cum_sum >= C * sum(sorted))[1]
hpd_idx = c(index, index + 1)
sorted[hpd_idx]
HPD = c(mu_values[which(density_values == sorted[hpd_idx][1])],
        mu_values[which(density_values == sorted[hpd_idx][2])])
HPD = c(min(HPD),max(HPD))
cat(HPD)

# HPD for tau
C = 0.95
tau_values = seq(0.001, 0.015, by = 5e-6)
density_values = density_tau(N, g, s, h, tau_values) * 1e60
sorted = sort(density_values, decreasing = TRUE)
cum_sum = cumsum(sorted)
index = which(cum_sum >= C * sum(sorted))[1]
hpd_idx = c(index, index + 1)
sorted[hpd_idx]
HPD = c(tau_values[which(density_values == sorted[hpd_idx][1])],
        tau_values[which(density_values == sorted[hpd_idx][2])])
HPD = c(min(HPD),max(HPD))
cat(HPD)


# HPD for sigma
C = 0.95
sigma_values = seq(0.05, 0.15, by = 0.0001)
density_values = density_sigma(N, g, s, h, sigma_values) * 0.5e61
sorted = sort(density_values, decreasing = TRUE)
cum_sum = cumsum(sorted)
index = which(cum_sum >= C * sum(sorted))[1]
hpd_idx = c(index, index + 1)
sorted[hpd_idx]
HPD = c(sigma_values[which(density_values == sorted[hpd_idx][1])],
        sigma_values[which(density_values == sorted[hpd_idx][2])])
HPD = c(min(HPD),max(HPD))
cat(HPD)
