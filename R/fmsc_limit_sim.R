rfmsc <- function(n, tau, bias_coef, tau_sd, efficient_sd){
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  tau_hat <- tau + tau_sd * z1
  efficient <- bias_coef * tau + efficient_sd * z2
  valid <- efficient_sd * z2 - bias_coef * tau_sd * z1
  return(data.frame(efficient, valid, tau_hat))
}

pfmsc <- function(x, tau, bias_coef, tau_sd, efficient_sd){
  tau_var <- tau_sd^2
  K_G <- pnorm(sqrt(2) - tau/tau_sd) - pnorm(-sqrt(2) - tau/tau_sd)
  G <- K_G * pnorm((x - bias_coef * tau) / efficient_sd)
  S <- matrix(c(bias_coef^2 * tau_var + efficient_sd^2, -bias_coef * tau_var,
                -bias_coef * tau_var, tau_var),
              byrow = TRUE, nrow = 2, ncol = 2)
  H_lower <- mvtnorm::pmvnorm(sigma = S,
                              lower = c(-Inf, -Inf),
                              upper = c(x, -tau_sd * sqrt(2) - tau))
  H_upper <- mvtnorm::pmvnorm(sigma = S,
                              lower = c(-Inf, tau_sd * sqrt(2) - tau),
                              upper = c(x, Inf))
  out <- G + H_upper + H_lower
  attr(out, which = "error") <- NULL
  attr(out, which = "msg") <- NULL
  return(out)
}

qfmsc <- function(p, tau, bias_coef, tau_sd, efficient_sd){
  u_min <- 1.1 * tau_sd * sqrt(2)
  #Find x1 such that pfmsc(x1, ...) - p > 0
  a1 <- qnorm(1 - (1 - p)/4, 0, tau_sd)
  u1 <- max(abs(tau - a1), abs(tau + a1), u_min)
  m1 <- bias_coef * tau - min(-bias_coef * u1, bias_coef * u1)
  x1 <- qnorm(1 - (1 - p)/4, m1, efficient_sd)
  #Find x2 such that pfmsc(x2, ...) - p < 0
  a2 <- qnorm(1 - p/4, 0, tau_sd)
  u2 <- max(abs(tau - a2), abs(tau + a2), u_min)
  m2 <- bias_coef * tau - max(-bias_coef * u2, bias_coef * u2)
  x2 <- qnorm(p/4, m2, efficient_sd)
  x_lower <- min(x1, x2)
  x_upper <- max(x1, x2)
  #Additional bounds, one of which is typically much tighter
  #although the other can be trivial. Substantially speeds up
  #the root finder.
  Kplus <- pnorm(sqrt(2) - tau/tau_sd)
  Kminus <- pnorm(-sqrt(2) - tau/tau_sd)
  z1 <- bias_coef * tau +
    efficient_sd * qnorm(fmscr::clip(p / (Kplus - Kminus), 0, 1))
  z2 <- bias_coef * tau + (bias_coef < 0) * bias_coef * tau_sd * sqrt(2) +
    efficient_sd * qnorm(fmscr::clip(1 + (p - 1) / Kplus, 0, 1))
  z_lower <- min(z1, z2)
  z_upper <- max(z1, z2)
  lower <- max(x_lower, z_lower)
  upper <- min(x_upper, z_upper)
  uniroot(function(x) pfmsc(x, tau, bias_coef, tau_sd, efficient_sd) - p,
          interval = c(lower, upper))$root
}
