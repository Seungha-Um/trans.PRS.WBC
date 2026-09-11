#' Fit a BART Model with Study Site-Specific Random Intercepts
#'
#' Fits a Bayesian Additive Regression Trees (BART) model with random
#' intercepts for study sites. This function is based on the SoftBart R
#' package and the BART framework described by Linero and Yang (2018).
#'
#' @param df_train A data matrix containing the training predictors.
#' @param Y_scaled A numeric vector containing the scaled outcome values.
#' @param num_tree An integer specifying the number of trees in the BART ensemble.
#' @param sim_num An integer specifying the number of MCMC iterations.
#' @param df_test A data matrix containing the test predictors.
#' @param num_site An integer specifying the number of study sites.
#' @param factor_site A vector indicating the study site for each observation
#'   in the training data.
#'
#' @return A list containing predicted values for the training and test data
#'  and variable importance measures.
#'
#' @export
#'

################### random intercept BART
BART_site <- function(df_train, Y_scaled, num_tree = 200, sim_num = 5000, df_test,
                      num_site, factor_site){

  hypers <- Hypers(X = df_train, Y = Y_scaled, normalize_Y = FALSE)
  hypers$num_tree <- num_tree
  opts <- Opts(update_s = FALSE)
  forest <- MakeForest(hypers, opts, warn = FALSE)

  NN <- dim(X_mat)[1]
  N_test <- dim(df_test)[1]

  mu_hat <- matrix(NA, sim_num, NN)
  mu_test <- matrix(NA, sim_num, N_test)
  alpha_mat <- matrix(NA, sim_num, num_site)

  #mu_test_hat <- matrix(NA, sim_num/step_size, df_train)
  ave_mu_test <- list()
  ave_mu_test_all <- list()

  beta_mat <- rep(1, num_site)
  ran_int <- model.matrix(~ factor_site - 1) %*% beta_mat

  Sigma_vec <- rep(NA, sim_num)
  Z_beta <- list()
  var_count <- matrix(NA, sim_num, dim(df_train)[2])

  for(s in 1:sim_num){
    R = Y_scaled - ran_int
    mu_hat[s,] <- forest$do_gibbs(df_train, R, df_train, 1)
    mu_test[s,] <- forest$do_predict(df_test)
    Sigma <- forest$get_sigma()
    var_count[s,] <- forest$get_counts()
    Sigma_vec[s] <- Sigma
    delta <- Y_scaled - mu_hat[s,]
    alpha_vec <- rep(NA, num_site)
    ran_int_new <- rep(NA, NN)

    for(g in name_site){
      ind_sites <- X_mat_site == g
      num_obs <- sum(ind_sites)
      Cov <- Sigma^2/num_obs
      mu_vec <- sum(delta[ind_sites])/num_obs
      alpha_g <- rnorm(1, mu_vec, Cov)
      alpha_vec[which(g == name_site)] <- alpha_g
      ran_int_new[ind_sites] <- alpha_g
    }
    alpha_mat[s,] <- alpha_vec
    ran_int <- ran_int_new
    if(s %% 100 == 0) cat("MCMC iter :", s, "\n")
  }

  bin_mat <- model.matrix(~ factor_site - 1)
  rand_ind_mat <- alpha_mat %*% t(bin_mat)
  temp_f <- mu_hat + rand_ind_mat
  temp_f <- temp_f * scale_Y + center_Y

  mu_test_hat <- mu_test * scale_Y + center_Y

  var_mat <- var_count[(sim_num/2):sim_num, ]
  var_im <- var_mat / rowSums(var_mat)

  return(list(mu_train_hat = temp_f,
              mu_test_hat,
              var_im))
}
