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

expression_to_integrate = function(x, N, g, y_bar, s, h) {
  numerator = gamma(N/2 + g)
  denominator = (0.5 * N * s^2 + N * (y_bar - x)^2 + h)^(N/2 + g)
  return(numerator / denominator * 1e65)
}

# Compute integral from -infinity to infinity
i_full = integrate(expression_to_integrate, lower = -100, upper = 100,
                           N = 50, g = 1, y_bar = 59.6, s = 11.1, h = 0.005)

# Compute integral from 55 to infinity
i_55 = integrate(expression_to_integrate, lower = 55, upper = 100,
                             N = 50, g = 1, y_bar = 59.6, s = 11.1, h = 0.005)

# Print results
print(i_full)
print(i_55)
(i_55$value/i_full$value)


###
# 3)
###
density_tau = function(N, g, s, h, tau)
{
  return = ( 1e60 * tau^(N/2 + g - 1) * exp(-0.5 * tau * (N * s^2 + h)) * sqrt(2 * pi / (N * tau)) )
}

g = 0
h = 0
tau_values = seq(0.001, 0.025, length.out = 100)
density_values = density_tau(N, g, s, h, tau_values)

# Plot
plot(tau_values, density_values, type = "l", xlab = "Tau", ylab = "density_tau",
     main = "Expression vs. Tau")

###
# 4)
###

joint_posterior = function(mu, tau, N, y, s, g, h) {
  exponent = -0.5 * tau * (N * s^2 + N * (y - mu)^2 + h)
  density = tau^((N / 2) + g - 1) * exp(exponent)
  return(density)
}

plot_joint_posterior_contour = function(N, y, s, g, h)
{
  # Define the range of mu and tau values
  mu_values = seq(50, 70, length.out = 100)
  tau_values = seq(0.01, 0.2, length.out = 100)
  
  # Create a grid of mu and tau values
  grid = expand.grid(mu = mu_values, tau = tau_values)
  
  # Calculate the unnormalized joint posterior density for each mu and tau
  grid$density = with(grid, joint_posterior(mu, tau, N, y, s, g, h))
  
  # Reshape the data for contour plotting
  density_matrix = matrix(grid$density, nrow = length(mu_values), ncol = length(tau_values), byrow = TRUE)
  
  # Plot the contour
  contour(mu_values, tau_values, density_matrix, main = "Joint Posterior Density Contour Plot", 
          xlab = expression(mu), ylab = expression(tau))
}
plot_joint_posterior_contour(N, y, s, 0.001, 0.001)
plot_joint_posterior_contour(N, y, s, 1, 0.005)
plot_joint_posterior_contour(N, y, s, 0, 0)
