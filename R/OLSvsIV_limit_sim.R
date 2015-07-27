pfmscOLSvsIV <- function(x, tau, pi_sq){
  tau_var <- (1 - pi_sq) / pi_sq
  tau_sd <- sqrt(tau_var)
  c <- tau / tau_sd
  sqrt2 <- 1.414214
  G <- pnorm(x - tau) * (pnorm(sqrt2 - c) - pnorm(-sqrt2 - c))
  S <- matrix(NA, 2, 2)
  S[1,1] <- 1 / pi_sq
  S[1,2] <- -tau_var
  S[2,1] <- -tau_var
  S[2,2] <- tau_var
  a <- sqrt2 * tau_sd
  H_upper <- mvtnorm::pmvnorm(mean = c(0, tau), sigma = S,
                              lower = c(-Inf, a), upper = c(x, Inf))
  H_lower <- mvtnorm::pmvnorm(mean = c(0, tau), sigma = S,
                              lower = c(-Inf, -Inf), upper = c(x, -a))
  out <- G + H_upper + H_lower
  attr(out, which = "error") <- NULL
  attr(out, which = "msg") <- NULL
  return(out)
}

qfmscOLSvsIV <- function(p, tau, pi_sq){
  sqrt2 <- 1.414214
  tau_sd <- sqrt((1 - pi_sq) / pi_sq)
  u_min <- 1.1 * tau_sd * sqrt2
  #Find x1 such that pfmscOLSvsIV(x1, tau, pi_sq) - p > 0
  a1 <- qnorm(1 - (1 - p)/4, 0, tau_sd)
  u1 <- max(abs(tau - a1), abs(tau + a1), u_min)
  x1 <- qnorm(1 - (1 - p)/4, tau + u1, 1)
  #Find x2 such that pfmscOLSvsIV(x2, tau, pi_sq) - p < 0
  a2 <- qnorm(1 - p/4, 0, tau_sd)
  u2 <- max(abs(tau - a2), abs(tau + a2), u_min)
  x2 <- qnorm(p/4, tau - u2, 1)
  x_lower <- min(x1, x2)
  x_upper <- max(x1, x2)
  #Potentially tighter bounds depending on parameters
  c <- tau / tau_sd
  Kplus <- pnorm(sqrt2 - c)
  Kminus <- pnorm(-sqrt2 - c)
  z1 <- tau + qnorm(fmscr::clip(p / (Kplus - Kminus), 0, 1))
  z2 <- tau + qnorm(fmscr::clip(1 + (p - 1) / Kplus, 0, 1))
  z_lower <- min(z1, z2)
  z_upper <- max(z1, z2)
  lower <- max(x_lower, z_lower)
  upper <- min(x_upper, z_upper)
  uniroot(function(x) pfmscOLSvsIV(x, tau, pi_sq) - p,
          interval = c(lower, upper))$root
}
