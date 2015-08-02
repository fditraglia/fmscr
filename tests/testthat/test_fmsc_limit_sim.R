context("Functions for FMSC Limit Simulation")

test_that("rfmsc", {
  set.seed(1234)
  pi_sq <- 0.3
  M1 <- matrix(c(1, 1, 0,
                 1, 1/pi_sq, -(1-pi_sq)/pi_sq,
                 0, -(1-pi_sq)/pi_sq, (1-pi_sq)/pi_sq),
               byrow = TRUE, 3, 3)
  sim1 <- rfmsc(n = 1e6, tau = 2, bias_coef = 1,
                tau_sd = sqrt((1 - pi_sq)/pi_sq) , efficient_sd = 1)
  M1sim <- cov(sim1)
  expect_less_than(sum(abs(M1sim - M1)), 0.1)
  expect_less_than(sum(abs(colMeans(sim1) - c(2, 0, 2))), 0.01)
  gamma <- 0.2
  M2 <- matrix(c(1/(gamma^2 + 1/9), 1/(gamma^2 + 1/9), 0,
                 1/(gamma^2 + 1/9), 9, -9 * gamma,
                 0, -9 * gamma, 1 + 9 * gamma^2),
                 byrow = TRUE, 3, 3)
  sim2 <- rfmsc(n = 1e6, tau = 2, bias_coef = gamma/(gamma^2 + 1/9),
                tau_sd = sqrt(1 + 9 * gamma^2),
                efficient_sd = sqrt(1/(gamma^2 + 1/9)))
  M2sim <- cov(sim2)
  expect_less_than(sum(abs(M2sim - M2)), 0.1)
  m2 <- c(2*gamma/(gamma^2 + 1/9), 0, 2)
  m2sim <- colMeans(sim2)
  expect_less_than(sum(abs(m2sim - m2)), 0.01)
})

