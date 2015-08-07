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
    efficient_sd * qnorm(fmscr::myclip(p / (Kplus - Kminus), 0, 1))
  z2 <- bias_coef * tau + (bias_coef < 0) * bias_coef * tau_sd * sqrt(2) +
    efficient_sd * qnorm(fmscr::myclip(1 + (p - 1) / Kplus, 0, 1))
  z_lower <- min(z1, z2)
  z_upper <- max(z1, z2)
  lower <- max(x_lower, z_lower)
  upper <- min(x_upper, z_upper)
  uniroot(function(x) pfmsc(x, tau, bias_coef, tau_sd, efficient_sd) - p,
          interval = c(lower, upper))$root
}


cover_naive <- function(alpha, tau, bias_coef, tau_sd, efficient_sd){
  tau_var <- tau_sd^2
  z <- qnorm(1 - alpha/2)
  shift <- bias_coef * tau / efficient_sd
  cover_efficient <- pnorm(z - shift) - pnorm(-z - shift)
  prob_efficient <- pnorm(sqrt(2) - tau/tau_sd) - pnorm(-sqrt(2) - tau/tau_sd)
  p1 <- prob_efficient * cover_efficient
  S <- matrix(c(bias_coef^2 * tau_var + efficient_sd^2, -bias_coef * tau_var,
                -bias_coef * tau_var, tau_var),
              byrow = TRUE, nrow = 2, ncol = 2)
  ell <- z * sqrt(efficient_sd^2 + bias_coef^2 * tau_var)
  p2lower <- mvtnorm::pmvnorm(sigma = S,
                              lower = c(-ell, -Inf),
                              upper = c(ell, -tau_sd * sqrt(2) - tau))
  p2upper <- mvtnorm::pmvnorm(sigma = S,
                              lower = c(-ell, tau_sd * sqrt(2) - tau),
                              upper = c(ell, Inf))
  out <- p1 + p2lower + p2upper
  attr(out, which = "error") <- NULL
  attr(out, which = "msg") <- NULL
  return(out)
}

shortestCI_fmsc <- function(alpha, tau, bias_coef, tau_sd, efficient_sd){
  g <- function(a_lower){
    lower <- qfmsc(a_lower, tau, bias_coef, tau_sd, efficient_sd)
    upper <- qfmsc(1 - (alpha - a_lower), tau, bias_coef, tau_sd, efficient_sd)
    upper - lower
  }
  a_star <- optimize(g, c(0.001, alpha - 0.001))$minimum
  lower <- qfmsc(a_star, tau, bias_coef, tau_sd, efficient_sd)
  upper <- qfmsc(1 - (alpha - a_star), tau, bias_coef, tau_sd, efficient_sd)
  return(c(lower, upper))
}

rel_width_FMSCinfeas <- function(alpha, tau, bias_coef, tau_sd, efficient_sd,
                                 equal.tailed = TRUE){
  valid_sd <- sqrt(efficient_sd^2 + bias_coef^2 * tau_sd^2)
  valid_width <- 2 * qnorm(1 - alpha/2) * valid_sd
  if(equal.tailed){
    lower <- qfmsc(alpha/2, tau, bias_coef, tau_sd, efficient_sd)
    upper <- qfmsc(1 - alpha/2, tau, bias_coef, tau_sd, efficient_sd)
  }else{
    CI <- shortestCI_fmsc(alpha, tau, bias_coef, tau_sd, efficient_sd)
    lower <- CI[1]
    upper <- CI[2]
  }
  return((upper - lower) / valid_width)
}

expect_rel_width <- function(tau, bias_coef, tau_sd, efficient_sd){
  pA <- pnorm(sqrt(2) - tau/tau_sd) - pnorm(-sqrt(2) - tau/tau_sd)
  rel_width <- sqrt(efficient_sd^2 / (efficient_sd^2 + bias_coef^2 * tau_sd^2))
  return(rel_width * pA + (1 - pA))
}

onestepCI <- function(alpha, tau, bias_coef, tau_sd, efficient_sd,
                     equal.tailed = TRUE, ndraws = 500){
  tau_draws <- rnorm(ndraws, tau, tau_sd)
  if(equal.tailed){
  qlower <- sapply(tau_draws, function(x)
    qfmsc(p = alpha/2, tau = x, bias_coef, tau_sd, efficient_sd))
  qupper <- sapply(tau_draws, function(x)
    qfmsc(p = 1 - alpha/2, tau = x, bias_coef, tau_sd, efficient_sd))
  }else{
    CIs <- sapply(tau_draws, function(x)
      shortestCI_fmsc(alpha, tau = x, bias_coef, tau_sd, efficient_sd))
    qlower <- CIs[1,]
    qupper <- CIs[2,]
  }
  widths <- abs(qupper - qlower)
  pupper <- sapply(qupper, function(q)
    pfmsc(q, tau, bias_coef, tau_sd, efficient_sd))
  plower <- sapply(qlower, function(q)
    pfmsc(q, tau, bias_coef, tau_sd, efficient_sd))
  cprob <- pupper - plower
  return(data.frame(avgwidth = mean(widths), coverage = mean(cprob)))
}

twostepCI <- function(alpha, tau, bias_coef, tau_sd, efficient_sd,
                      a1 = alpha/2, ndraws = 500){
  tau_draws <- rnorm(ndraws, tau, tau_sd)
  ME <- qnorm(1 - a1/2) * tau_sd
  a2 <- alpha - a1
  get_lower <- function(tau_hat){
    optimize(function(x) qfmsc(a2/2, tau = x, bias_coef, tau_sd, efficient_sd),
           lower = tau_hat - ME, upper = tau_hat + ME, tol = 0.01)$objective
  }
  get_upper <- function(tau_hat){
    optimize(function(x) qfmsc(1 - a2/2, tau = x, bias_coef, tau_sd, efficient_sd),
           lower = tau_hat - ME, upper = tau_hat + ME, tol = 0.01,
           maximum = TRUE)$objective
  }
  qlower <- sapply(tau_draws, get_lower)
  qupper <- sapply(tau_draws, get_upper)
  widths <- abs(qupper - qlower)
  pupper <- sapply(qupper, function(q)
    pfmsc(q, tau, bias_coef, tau_sd, efficient_sd))
  plower <- sapply(qlower, function(q)
    pfmsc(q, tau, bias_coef, tau_sd, efficient_sd))
  cprob <- pupper - plower
  return(data.frame(avgwidth = mean(widths), coverage = mean(cprob)))
}


get_twostepCI <- function(alpha, tau_hat, bias_coef, tau_sd, efficient_sd,
                      a1 = alpha/2){
  ME <- qnorm(1 - a1/2) * tau_sd
  a2 <- alpha - a1
  lower <- optimize(function(x)
    qfmsc(a2/2, tau = x, bias_coef, tau_sd, efficient_sd),
    lower = tau_hat - ME, upper = tau_hat + ME, tol = 0.01)$objective

  upper <- optimize(function(x)
    qfmsc(1 - a2/2, tau = x, bias_coef, tau_sd, efficient_sd),
    lower = tau_hat - ME, upper = tau_hat + ME, tol = 0.01,
    maximum = TRUE)$objective

  return(c(lower, upper))
}













