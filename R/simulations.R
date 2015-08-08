set.seed(4324)
rho <- 0.1
pi_sq <- 0.2
N <- 100
alpha <- 0.2
n_reps <- 200

sim_data <- as.data.frame(sim_OLSvsIV(rho, pi_sq, N, n_reps))
sim_data$tau_sd <- sqrt(sim_data$tau_var)
sim_data$bias_coef <- 1 / sim_data$s_x_sq
sim_data$efficient_sd <- with(sim_data, sqrt(s_e_sq_tsls / s_x_sq))
sim_data$fmsc_ols <- with(sim_data, abs(tau_hat) < sqrt(2) * tau_sd)
sim_data$b_fmsc <- with(sim_data, ifelse(fmsc_ols, b_ols, b_tsls))
qz <- qnorm(1 - alpha / 2)

#======================= OLS Intervals
l_ols <- with(sim_data, b_ols - qz * SE_ols)
u_ols <- with(sim_data, b_ols + qz * SE_ols)
c_ols <- coverage_prob(cbind(l_ols, u_ols), 0.5)
w_ols <- mean(abs(u_ols - l_ols))

#======================= TSLS Intervals
l_tsls <- with(sim_data, b_tsls - qz * SE_tsls)
u_tsls <- with(sim_data, b_tsls + qz * SE_tsls)
c_tsls <- coverage_prob(cbind(l_tsls, u_tsls), 0.5)
w_tsls <- mean(abs(u_tsls - l_tsls))

#======================= Naive Intervals
sim_data$SE_naive <- with(sim_data, ifelse(fmsc_ols, SE_ols, SE_tsls))
l_naive <- with(sim_data, b_fmsc - qz * SE_naive)
u_naive <- with(sim_data, b_fmsc + qz * SE_naive)
c_naive <- coverage_prob(cbind(l_naive, u_naive), 0.5)
w_naive <- mean(abs(u_naive - l_naive))

#======================= Equal-tailed 1-step Intervals
l_onestep_equal <- unlist(Map(function(tau_hat, bias_coef, tau_sd, efficient_sd)
  qfmsc(alpha/2, tau_hat, bias_coef, tau_sd, efficient_sd),
  sim_data$tau_hat, sim_data$bias_coef, sim_data$tau_sd, sim_data$efficient_sd))
u_onestep_equal <- unlist(Map(function(tau_hat, bias_coef, tau_sd, efficient_sd)
  qfmsc(1 - alpha/2, tau_hat, bias_coef, tau_sd, efficient_sd),
  sim_data$tau_hat, sim_data$bias_coef, sim_data$tau_sd, sim_data$efficient_sd))
onestep_equal <- cbind(l_onestep_equal, u_onestep_equal)
onestep_equal <- t(apply(sim_data$b_fmsc - onestep_equal / sqrt(N), 1, sort))

c_1equal <- coverage_prob(onestep_equal, 0.5)
w_1equal <- mean(abs(onestep_equal[,2] - onestep_equal[,1]))

#======================= Shortest 1-step Intervals
onestep_short <- Map(function(tau_hat, bias_coef, tau_sd, efficient_sd)
  shortestCI_fmsc(alpha, tau_hat, bias_coef, tau_sd, efficient_sd),
  sim_data$tau_hat, sim_data$bias_coef, sim_data$tau_sd, sim_data$efficient_sd)
onestep_short <- do.call(rbind, onestep_short)
onestep_short <- t(apply(sim_data$b_fmsc - onestep_short / sqrt(N), 1, sort))

c_1short <- coverage_prob(onestep_short, 0.5)
w_1short <- mean(abs(onestep_short[,2] - onestep_short[,1]))

#======================= 2-step Intervals (a1 = a2 = 0.5 * alpha)
twostep_equal <- Map(function(tau_hat, bias_coef, tau_sd, efficient_sd)
  get_twostepCI(alpha, tau_hat, bias_coef, tau_sd, efficient_sd,
                a1 = 0.5 * alpha), sim_data$tau_hat, sim_data$bias_coef,
  sim_data$tau_sd, sim_data$efficient_sd)
twostep_equal <- do.call(rbind, twostep_equal)
twostep_equal <- t(apply(sim_data$b_fmsc - twostep_equal / sqrt(N), 1, sort))

c_2equal <- coverage_prob(twostep_equal, 0.5)
w_2equal <- mean(abs(twostep_equal[,2] - twostep_equal[,1]))

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
