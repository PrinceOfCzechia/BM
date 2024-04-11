rm(list=ls())

library(MASS)
library(pracma)
library(coda)

N = 50
y = 59.6
s = 11.1

g1 = 0.001
h1 = 0.001
g2 = 1
h2 = 0.005
g3 = 0
h3 = 0

colors = c("red", "forestgreen", "dodgerblue2")


###
# 1)
###

# mean
c( y-qt(0.975,N-1)*s/sqrt(N), y+qt(0.975,N-1)*s/sqrt(N) )

# std deviation
c( sqrt((N-1)*s^2/qchisq(0.975,N-1)), sqrt((N-1)*s^2/qchisq(0.025,N-1)) )

# inverse variance
c( qchisq(0.025,N-1)/((N-1)*11.6^2), qchisq(0.975,N-1)/((N-1)*11.6^2) )


###
# 2)
###

density_mu = function(mu, g=g1, h=h1)
{
  return = (1/2*((n-1)*s^2+n*(y-mu)^2+2*h))^(-n/2-g)*gamma(n/2+g)*((n-1)*s^2/2+h)^((n-1)/2+g)/gamma((n-1)/2+g)*sqrt(n/(2*pi))
}

mu_values = seq(50, 70, by=0.02)
f_mu_1 = density_mu(mu_values, g1, h1)
f_mu_2 = density_mu(mu_values, g2, h2)
f_mu_3 = density_mu(mu_values, g3, h3)

plot(f_mu_2 ~ mu_values, col = colors[2], lwd = 3, type = 'l',
     xlab = 'mu', ylab = 'f(mu,g,h)', main = 'Density of mu for different parameters')
lines(f_mu_1 ~ mu_values, col = colors[1], lwd = 3)
lines(f_mu_3 ~ mu_values, col = colors[3], lwd = 3, lty = 'dashed')
legend("topright", legend = c("g=0.001, h=0.001", "g=1, h=0.005", "g=0, h=0"),
       col = colors,
       lty = c(1, 1, 2), lwd = c(3, 3, 3),
       title = "Parameters")

i_55 = integrate(density_mu, lower = 55, upper = 100, g = g1, h = h1)
cat(i_55$value)
i_55 = integrate(density_mu, lower = 55, upper = 100, g = g2, h = h2)
cat(i_55$value)
i_55 = integrate(density_mu, lower = 55, upper = 100, g = g3, h = h3)
cat(i_55$value)


###
# 3)
###

# tau
density_tau = function(tau, g=g1, h=h1)
{
  return (tau^((N-1)/2+g-1)*exp(-tau*((N-1)/2*s^2+h))*((N-1)/2*s^2+h)^((N-1)/2+g)/gamma((N-1)/2+g))
}

tau_values = seq(0.001, 0.020, length.out = 200)
f_tau_1 = density_tau(tau_values, g1, h1)
f_tau_2 = density_tau(tau_values, g2, h2)
f_tau_3 = density_tau(tau_values, g3, h3)

plot(f_tau_2 ~ tau_values, col = colors[2], lwd = 3, type = 'l',
     xlab = 'tau', ylab = 'f(tau,g,h)', main = 'Density of tau for different parameters')
lines(f_tau_1 ~ tau_values, col = colors[1], lwd = 3)
lines(f_tau_3 ~ tau_values, col = colors[3], lwd = 3, lty = 'dashed')
legend("topright", legend = c("g=0.001, h=0.001", "g=1, h=0.005", "g=0, h=0"),
       col = colors,
       lty = c(1, 1, 2), lwd = c(3, 3, 3),
       title = "Parameters")

# sigma
density_sigma = function(sigma, g=g1, h=h1) # using trafo of density_tau
{
  return (density_tau(1/(sigma^2), g, h)*2/(sigma^3))
}

sigma_values = seq(5, 20, length.out = 200)
f_sigma_1 = sapply(sigma_values, density_sigma, g = g1, h = h1)
f_sigma_2 = sapply(sigma_values, density_sigma, g = g2, h = h2)
f_sigma_3 = sapply(sigma_values, density_sigma, g = g3, h = h3)

plot(f_sigma_2 ~ sigma_values, col = colors[2], lwd = 3, type = 'l',
     xlab = 'sigma', ylab = 'f(sigma,g,h)', main = 'Density of sigma for different parameters')
lines(f_sigma_1 ~ sigma_values, col = colors[1], lwd = 3)
lines(f_sigma_3 ~ sigma_values, col = colors[3], lwd = 3, lty = 'dashed')
legend("topright", legend = c("g=0.001, h=0.001", "g=1, h=0.005", "g=0, h=0"),
       col = colors,
       lty = c(1, 1, 2), lwd = c(3, 3, 3),
       title = "Parameters")


###
# 6)
###

# quantile function is needed for ET

# mu
CDF_mu = function(mu, g=g1, h=h1)
{ 
  return (integrate(function(m) density_mu(m, g, h), lower= 0, upper = mu)$value)
}
q_mu = function(m, g=g1, h=h1)
{
  ifelse( m>0.5,
          uniroot(function(mu) {CDF_mu(mu, g, h) - m}, interval = c(59.6, 100))$root,
          uniroot(function(mu) {CDF_mu(mu, g, h) - m}, interval = c(-100, 59.6))$root
        )
}
ET_mu_1 = c(q_mu(0.025, g1, h1),q_mu(0.975, g1, h1))
ET_mu_2 = c(q_mu(0.025, g2, h2),q_mu(0.975, g2, h2))
ET_mu_3 = c(q_mu(0.025, g3, h3),q_mu(0.975, g3, h3))

# tau
CDF_tau = function(tau, g=g1, h=h1)
{ 
  return (integrate(function(t) density_tau(t, g, h), lower= 0.001, upper = tau)$value)
}
q_tau = function(m, g=g1, h=h1)
{
  ifelse( m<0.5,
          uniroot(function(tau) {CDF_tau(tau, g, h) - m}, interval = c(0.001, 0.008007))$root,
          uniroot(function(tau) {CDF_tau(tau, g, h) - m}, interval = c(0.008007, 0.015))$root
  )
}
ET_tau_1 = c(q_tau(0.025),q_tau(0.975))
ET_tau_2 = c(q_tau(0.025, g2, h2),q_tau(0.975, g2, h2))
ET_tau_3 = c(q_tau(0.025, g3, h3),q_tau(0.975, g3, h3))

# sigma
CDF_sigma = function(sigma, g=g1, h=h1)
{ 
  return (integrate(function(s) density_sigma(s, g, h), lower= 0, upper = sigma)$value)
}
q_sigma = function(m, g=g1, h=h1)
{
  ifelse( m<0.5,
          uniroot(function(sigma) {CDF_sigma(sigma) - m}, interval = c(5, 10.9))$root,
          uniroot(function(sigma) {CDF_sigma(sigma) - m}, interval = c(10.9, 20))$root
  )
}
ET_sigma_1 = c(q_sigma(0.025),q_sigma(0.975))
ET_sigma_2 = c(q_sigma(0.025, g2, h2),q_sigma(0.975, g2, h2))
ET_sigma_3 = c(q_sigma(0.025, g3, h3),q_sigma(0.975, g3, h3))


###
# 7)
###

# HPD for mu
HPD_mu = function(g=g1, h=h1)
{
C = 0.95
mu_values = seq(50, 70, by = 0.0001)
density_values = density_mu(mu_values, g, h)
sorted = sort(density_values, decreasing = TRUE)
cum_sum = cumsum(sorted)
index = which(cum_sum >= C * sum(sorted))[1]
hpd_idx = c(index, index + 1)
sorted[hpd_idx]
HPD = c(mu_values[which(density_values == sorted[hpd_idx][1])],
        mu_values[which(density_values == sorted[hpd_idx][2])])
return = c(min(HPD),max(HPD))
}
HPD_mu_1 = HPD_mu()
HPD_mu_2 = HPD_mu(g=g2, h=h2)
HPD_mu_3 = HPD_mu(g=g3, h=h3)

# HPD for tau
HPD_tau = function(g=g1, h=h1)
{
  C = 0.95
  tau_values = seq(0.001, 0.015, by = 0.00001)
  density_values = density_tau(tau_values, g, h)
  sorted = sort(density_values, decreasing = TRUE)
  cum_sum = cumsum(sorted)
  index = which(cum_sum >= C * sum(sorted))[1]
  hpd_idx = c(index, index + 1)
  sorted[hpd_idx]
  HPD = c(tau_values[which(density_values == sorted[hpd_idx][1])],
          tau_values[which(density_values == sorted[hpd_idx][2])])
  return = c(min(HPD),max(HPD))
}
HPD_tau_1 = HPD_tau()
HPD_tau_2 = HPD_tau(g=g2, h=h2)
HPD_tau_3 = HPD_tau(g=g3, h=h3)

# HPD for sigma
HPD_sigma = function(g=g1, h=h1)
{
  C = 0.95
  sigma_values = seq(3, 18, by = 0.00001)
  density_values = density_sigma(sigma_values, g, h)
  sorted = sort(density_values, decreasing = TRUE)
  cum_sum = cumsum(sorted)
  index = which(cum_sum >= C * sum(sorted))[1]
  hpd_idx = c(index, index + 1)
  sorted[hpd_idx]
  HPD = c(sigma_values[which(density_values == sorted[hpd_idx][1])],
          sigma_values[which(density_values == sorted[hpd_idx][2])])
  return = c(min(HPD),max(HPD))
}
HPD_sigma_1 = HPD_sigma()
HPD_sigma_2 = HPD_sigma(g=g2, h=h2)
HPD_sigma_3 = HPD_sigma(g=g3, h=h3)

sigma_values = seq(5, 20, length.out = 200)
f_sigma_3 = sapply(sigma_values, density_sigma, g = g3, h = h3)
plot(f_sigma_3 ~ sigma_values, col = colors[3], lwd = 3, type = 'l',
     xlab = 'sigma', ylab = 'f(sigma,g,h)',
     main = 'ET and HPD for sigma')
abline(v=ET_sigma_1, col='red', lwd=2)
abline(v=HPD_sigma_1, col='magenta', lwd=2)
legend("topright", legend = c("density for g=0, h=0", "ET_sigma_3", "HPD_sigma_3"),
       col = c('dodgerblue2', 'red', 'magenta'),
       lty = c(1, 1, 1),
       lwd = c(3, 2, 2))


###
# 8)
###

gen_posterior = function(N_samples, g=g1, h=h1)
{
  N_samples = 1000
  U1 = runif(N_samples)
  U2 = runif(N_samples)
  
  mu_sample = numeric(N_samples)
  tau_sample = numeric(N_samples)
  
  for (i in 1:N_samples)
  {
    mu_sample[i] = q_mu(U1[i], g, h)
    tau_sample[i] = q_tau(U2[i], g, h)
  }
  result_list = list(mu_sample = mu_sample, tau_sample = tau_sample)
  return(result_list)
}

MC_sample = gen_posterior( 100000 ) # occasionally gives error, then just retry
hist( MC_sample$mu_sample )
hist( MC_sample$tau_sample )

# sample CI from the Monte Carlo
CI_mu = c( quantile( MC_sample$mu_sample, 0.025 ), quantile( MC_sample$mu_sample, 0.975 ) )
CI_tau = c( quantile( MC_sample$tau_sample, 0.025 ), quantile( MC_sample$tau_sample, 0.975 ) )
