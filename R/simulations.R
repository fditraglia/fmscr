rho <- 0.2
pi_sq <- 0.2
N <- 100
alpha <- 0.05
n_reps <- 10

sim_data <- as.data.frame(sim_OLSvsIV(rho, pi_sq, N, n_reps))
sim_data$tau_sd <- sqrt(sim_data$tau_var)
sim_data$bias_coef <- 1 / sim_data$s_x_sq
sim_data$efficient_sd <- with(sim_data, sqrt(s_e_sq_tsls / g_sq))
sim_data$fmsc_ols <- with(sim_data, abs(tau_hat) < sqrt(2) * tau_sd)
sim_data$b_fmsc <- with(sim_data, ifelse(fmsc_ols, b_ols, b_tsls))


#======================= OLS Intervals


#======================= TSLS Intervals

#======================= Naive Intervals
sim_data$SE_naive <- with(sim_data, ifelse(fmsc_ols, SE_ols, SE_tsls))
qz <- qnorm(1 - alpha / 2)
l_naive <- with(sim_data, b_fmsc - qz * SE_naive)
u_naive <- with(sim_data, b_fmsc + qz * SE_naive)
naive <- cbind(l_naive, u_naive)

cnaive <- coverage_prob(naive, 0.5)
wnaive <- mean(abs(naive[,2] - naive[,1]))

#======================= Equal-tailed 1-step Intervals
l_onestep_equal <- unlist(Map(function(tau_hat, bias_coef, tau_sd, efficient_sd)
  qfmsc(alpha/2, tau_hat, bias_coef, tau_sd, efficient_sd),
  sim_data$tau_hat, sim_data$bias_coef, sim_data$tau_sd, sim_data$efficient_sd))
u_onestep_equal <- unlist(Map(function(tau_hat, bias_coef, tau_sd, efficient_sd)
  qfmsc(1 - alpha/2, tau_hat, bias_coef, tau_sd, efficient_sd),
  sim_data$tau_hat, sim_data$bias_coef, sim_data$tau_sd, sim_data$efficient_sd))
onestep_equal <- cbind(l_onestep_equal, u_onestep_equal)
onestep_equal <- t(apply(sim_data$b_fmsc - onestep_equal / sqrt(N), 1, sort))

c1equal <- coverage_prob(onestep_equal, 0.5)
w1equal <- mean(abs(onestep_equal[,2] - onestep_equal[,1]))

#======================= Shortest 1-step Intervals
onestep_short <- Map(function(tau_hat, bias_coef, tau_sd, efficient_sd)
  shortestCI_fmsc(alpha, tau_hat, bias_coef, tau_sd, efficient_sd),
  sim_data$tau_hat, sim_data$bias_coef, sim_data$tau_sd, sim_data$efficient_sd)
onestep_short <- do.call(rbind, onestep_short)
onestep_short <- t(apply(sim_data$b_fmsc - onestep_short / sqrt(N), 1, sort))

c1short <- coverage_prob(onestep_short, 0.5)
w1short <- mean(abs(onestep_short[,2] - onestep_short[,1]))

#======================= 2-step Intervals (a1 = a2 = 0.5 * alpha)
twostep_equal <- Map(function(tau_hat, bias_coef, tau_sd, efficient_sd)
  get_twostepCI(alpha, tau_hat, bias_coef, tau_sd, efficient_sd,
                a1 = 0.5 * alpha), sim_data$tau_hat, sim_data$bias_coef,
  sim_data$tau_sd, sim_data$efficient_sd)
twostep_equal <- do.call(rbind, twostep_equal)
twostep_equal <- t(apply(sim_data$b_fmsc - twostep_equal / sqrt(N), 1, sort))

c2equal <- coverage_prob(onestep_short, 0.5)
w2equal <- mean(abs(twostep_equal[,2] - twostep_equal[,1]))

#======================= 2-step Intervals (a1 = 0.25 * alpha, a2 = 0.75 * alpha)
twostep_tau_wide <- Map(function(tau_hat, bias_coef, tau_sd, efficient_sd)
  get_twostepCI(alpha, tau_hat, bias_coef, tau_sd, efficient_sd,
                a1 = 0.25 * alpha), sim_data$tau_hat, sim_data$bias_coef,
  sim_data$tau_sd, sim_data$efficient_sd)
twostep_tau_wide <- do.call(rbind, twostep_tau_wide)
twostep_tau_wide <- t(apply(sim_data$b_fmsc - twostep_tau_wide / sqrt(N),
                            1, sort))

#======================= 2-step Intervals (a1 = 0.75 * alpha, a2 = 0.25 * alpha)
twostep_tau_narrow <- Map(function(tau_hat, bias_coef, tau_sd, efficient_sd)
  get_twostepCI(alpha, tau_hat, bias_coef, tau_sd, efficient_sd,
                a1 = 0.75 * alpha), sim_data$tau_hat, sim_data$bias_coef,
  sim_data$tau_sd, sim_data$efficient_sd)
twostep_tau_narrow <- do.call(rbind, twostep_tau_narrow)
twostep_tau_narrow <- t(apply(sim_data$b_fmsc - twostep_tau_narrow / sqrt(N),
                              1, sort))
