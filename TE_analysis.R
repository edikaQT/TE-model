#########################################################################################
#########################################################################################
#
#      DOCUMENTATION OF THE MAJOR FUNCTIONS (useful to those who program)
#      (skip to line 1111 for where a non-programmer can adjust the inputs)
#
#########################################################################################
#
# True and Error model estimation
# 
# Function 'TE_estimation_data” fits the TE model to a dataset by 
# minimizing a chi-squared statistic. It estimates the probabilities 
# of eight true preference patterns denoted by 
# p111, p112, p121, p122, p211, p212, p221, and p222, 
# along with three errors e1, e2 and e3, for the three choice problems. 
# The estimation assumes strict positive probabilities and errors, which 
# are constrained to be no greater than a half. Syntax is as follows:
# 
# TE_estimation_data(data)
# 
# Where "data" is a data frame with 17 columns describing the frequencies 
# of responses in G choices and in both, G and F, choices for patterns 
# 111, 112, 121, 122, 211, 212, 221, and 222. 
# Each row contains data for a different group or individual. Also, columns should be named as:  
# "ID" (the participant ID), "G111", "G112", "G121", "G122", "G211", "G212", "G221", "G222", 
# "both111", "both112", "both121", "both122", "both211", "both212", "both221" and "both222". 
# 
# The program generates the following components are estimated for each dataset: 
# probability estimates, error estimates, a Chi-squared statistic, a p value, 
# observed and predicted proportions of response patterns, and observed and predicted 
# frequencies of response patterns. These components can be obtained by typing:
# 
# ... [Case Number]$parameters
# ... [Case Number]$errors
# ... [Case Number]$chisq
# ... [Case Number]$p.value
# ... [Case Number]$contrast.of.proportions
# ... [Case Number]$contrast.of.frequencies
# in which, the "..." is the name of the list storing the results.
# 
#
# Monte Carlo simulation of True and Error model
# 
# Function "TE_MonteCarlo_data" performs simulations of the TE model by 
# fitting samples created under the null of the model (i.e., the predictions 
# of the model estimated from the empirical data). Chi-squared statistics are computed 
# for each simulated sample and are then used to determine a Monte Carlo p value 
# for the model. 
# Two types of simulations are performed for each participant: 
# Re-fitted Monte Carlo simulations and Conservative Monte Carlo simulations. 
#   With Re-fitted Monte Carlo simulations, TE fits new parameters in each sample and calculates a Chi squared statistic.
#   With Conservative Monte Carlo simulations, each Chi squared statistic compares of each new sample with the original predictions 
#   of the model (i.e., no re-fitting is performed).
# 
# Syntax is as follows:
# 
# TE_MonteCarlo_data(data, number_samples)
# 
# A list containing the following components is constructed for each case: 
# the parameter estimates of each sample (for Re-fitted Monte Carlo), 
# two Chi-squared statistics in each sample (For Re-fitted and Conservative Monte Carlo), 
# the model p values (for Re-fitted and Conservative Monte Carlo), 
# the simulated samples, and the original Chi-squared statistic. These components are specified as follows:
# 
# ... [Case Number]$parameters_Re.fitted.MC
# ... [Case  Number]$errors_Re.fitted.MC
# ... [Case  Number]$chisq_Re.fitted.MC
# ... [Case  Number]$p.value_Re.fitted.MC
# ... [Case  Number]$chisq_Conservative.MC
# ... [Case  Number]$p.value_Conservative.MC
# ... [Case  Number]$samples
# ... [Case  Number]$original.chisq
# 
# in which, the "..." is the name of the list storing the results.
# 
# Bootstrapping simulation of True and Error model
# 
# Function "TE_Bootstrapping.CI_data" performs simulations of the TE model 
# by fitting random samples resampled from the original data with replacement. 
#  Parameters are estimated for each simulated sample and are then used to 
# estimate 95% confidence intervals using the percentile method. Syntax:
# 
# TE_Bootstrapping.CI_data (data, number_samples)
# 
# This function generates a list containing the simulations performed for each case. 
# The following components are computed: probability estimates for each sample, 
# error estimates for each sample, simulated samples, and confidence intervals 
# for each parameter estimate. These components can be obtained by typing:
# 
# ... [Case Number]$boot.parameters
# ... [Case Number]$boot.errors
# ... [Case Number]$boot.samples
# ... [Case Number]$confidence.intervals
# 
# Independence model estimation
# 
# Function "Independence_estimation_data" fits the Independence model by 
# minimizing a chi-squared index. This model assumes that the predicted probability 
# of observing, for instance, the response pattern 112 is given by: 
# P112=(p1xx) * (px1x) * (1-pxx1).
# It estimates three marginal choice probabilities denoted by p1xx, px1x and pxx1. 
# These probabilities relate to proportions for response 1 in choice problems 
# G+ versus G0, G0 versus G-, and G+ versus G-, respectively. 
# The estimation assumes probabilities between 0 and 1. Syntax:
# 
# Independence_estimation_data(data)
# 
# The function generates a list with following components: 
# marginal choice probability estimates, a Chi-squared statistic, a p value, 
# observed and predicted proportions of response patterns, and observed and predicted frequencies 
# of response patterns. These components can be obtained by typing:
# 
# ... [Case Number]$parameters
# ... [Case Number]$chisq
# ... [Case Number]$p.value
# ... [Case Number]$contrast.of.proportions
# ... [Case Number]$contrast.of.frequencies
# 
# Monte Carlo simulation of Independence model
# 
# Function "Independence_MonteCarlo_data" perform simulations of the 
# Independence model by fitting samples created under the null of that model 
# (i.e., the predictions of the model given the empirical data). 
# Chi-squared statistics are computed for each simulated sample and are 
# then used to determine a Monte Carlo p value for the model. 
# Two types of simulations are performed for each participant: 
# Re-fitted Monte Carlo simulations and Conservative Monte Carlo simulations. 
#   With Re-fitted Monte Carlo simulations,  new parameters are estimated in each sample for the Chi squared statistic.
#   With Conservative Monte Carlo simulations, the original parameters are used for all samples in computing statistic.
# 
# Independence_MonteCarlo_data (data, number_samples)
# 
# A list containing the following components is constructed for each case: 
# the parameter estimates of each sample (for Re-fitted Monte Carlo), 
# two Chi-squared statistics in each sample (For Re-fitted and Conservative Monte Carlo), 
# the model p values (for Re-fitted and Conservative Monte Carlo), 
# the simulated samples, and the original Chi-squared statistic. These components are specified as follows:
# 
# ... [Case Number]$parameters_Re.fitted.MC
# ... [Case  Number]$errors_Re.fitted.MC
# ... [Case  Number]$chisq_Re.fitted.MC
# ... [Case  Number]$p.value_Re.fitted.MC
# ... [Case  Number]$chisq_Conservative.MC
# ... [Case  Number]$p.value_Conservative.MC
# ... [Case  Number]$samples
# ... [Case  Number]$original.chisq
# 
# in which, the "..." is the name of the list storing the results.
#######################################################################################
#
##########################################################################################
#
# FUNCTIONS APPLIED IN THE ESTIMATION OF TE MODEL
#
# - To estimate TE parameters
#    +++ true_error: sets the model equations
#    +++ wrapper: sets the parameters to be used in the optimization algorithm
#    +++ TE_estimation: estimates parameters using responses of an individual
#    +++ TE_estimation_data: estimates parameters for all the individuals presented 
#                            in the researcher's data file
#    +++ aggregate_estimates: accumulates results obtained by "TE_estimation_data"
#                             in a data frame.
#
# - To perform Monte Carlo simulations
#    +++ MonteCarlo_TE: performs a number of Monte Carlo simulations given a vector of 
#                       predicted frequencies of the model (i.e., the null of the model that is
#                       used to create the samples)
#    +++ MonteCarlo_simulation_TE: provides the vector of predicted frequencies to
#                                  function "MonteCarlo_TE"
#    +++ TE_MonteCarlo_data: computes simulations for all individuals presented 
#                            in the reseaercher's data file
#    +++ aggregate_MonteCarlo_p.value: accumulates all Monte Carlo p values estimated by
#                                      "TE_MonteCarlo_data". 
#
# - To perform Bootstrapping simulations to obtain confidence intervals
#    +++ CI_estimates_bootstrap: defines the model to be used in the simulation
#    +++ CI_Bootstrap_TE: performs bootstrapping simulations using responses of an individual
#    +++ TE_Bootstrapping.CI_data: performs bootstrapping simulations for all the individuals
#                                  preesented in the researcher's data file
#    +++ aggregate_boot.confidence.intervals: aggregates the results (confidence intervals)
#                                             in a data frame
#
###########################################################################################


# FUNCTIONS TO ESTIMATE TE PARAMETERS
# -------------------------------------

# +++ Function "true_error"

true_error<-function( p111, p112,  p121, p122,  p211,   p212,   p221, p222, 
                      error1, error2, error3, G, G_and_F ){
  
  # We ensure that the parameters used by the model are strictly positive
  p111 <- ifelse(p111<0.00001,0.00001,p111)
  p112 <- ifelse(p112<0.00001,0.00001,p112)
  p121 <- ifelse(p121<0.00001,0.00001,p121)
  p122 <- ifelse(p122<0.00001,0.00001,p122)
  p211 <- ifelse(p211<0.00001,0.00001,p211)
  p212 <- ifelse(p212<0.00001,0.00001,p212)
  p221 <- ifelse(p221<0.00001,0.00001,p221)
  p222 <- ifelse(p222<0.00001,0.00001,p222)
  error1 <- ifelse(error1<0.00001,0.00001,error1)
  error2 <- ifelse(error2<0.00001,0.00001,error2)
  error3 <- ifelse(error3<0.00001,0.00001,error3)
  
  sum <- sum(c(p111, p112, p121, p122, p211, p212, p221,p222))
  true_prob <- c(p111, p112, p121, p122, p211, p212, p221,p222)/sum
  G_exclusive <- G-G_and_F
  Total_G <- sum(G) #Number of participants
  
  # Creating a vector with all possible combinations of errors in three choices
  error <- c((1-error1)*(1-error2)*(error3), 
             (1-error1)*(1-error2)*(1-error3), 
             (1-error1)*(error2)*(error3),
             (1-error1)*(error2)*(1-error3),
             (error1)*(1-error2)*(error3),
             (error1)*(1-error2)*(1-error3),
             (error1)*(error2)*(error3), 
             (error1)*(error2)*(1-error3))
  
  order_111 <- c(2,1,4,3,6,5,8,7)
  order_112 <- c(1,2,3,4,5,6,7,8)
  order_121 <- c(4,3,2,1,8,7,6,5)
  order_122 <- c(3,4,1,2,7,8,5,6)
  order_211 <- c(6,5,8,7,2,1,4,3)
  order_212 <- c(5,6,7,8,1,2,3,4)
  order_221 <- c(8,7,6,5,4,3,2,1)
  order_222 <- c(7,8,5,6,3,4,1,2)
  
  # Ordering of errors for each of the eight equations estimating true preferences
  # for each pattern
  
  error_m <- rbind(error[order_111], 
                   error[order_112],
                   error[order_121], 
                   error[order_122],
                   error[order_211], 
                   error[order_212],
                   error[order_221], 
                   error[order_222])
  
  
  # Errors are squared when subjects chose twice the same pattern              
  error2_m <- error_m*error_m 
  
  # Probability of observing a pattern in G choices
  prob_observe_G <- rowSums(t(t(error_m) * true_prob)) 
  
  # Probability of observing a pattern in both G and F choices
  prob_observe_GF <- rowSums(t(t(error2_m) * true_prob)) 
  
  # Probability of observing a pattern in G choices but not in F choices
  prob_observe_only_G <- prob_observe_G-prob_observe_GF 
  
  # Computing Chi-squared values
  Chi_GF_vector <- ((G_and_F-prob_observe_GF*Total_G)^2)/(prob_observe_GF*Total_G)
  Chi_GF <- sum( Chi_GF_vector)
  
  Chi_only_G_vector <- ((G_exclusive - prob_observe_only_G*Total_G)^2)/
    (prob_observe_only_G*Total_G)
  Chi_only_G <- sum(Chi_only_G_vector)
  
  Chi_total = Chi_GF+Chi_only_G
  
  list (Chi_total,   c(prob_observe_only_G,prob_observe_GF))
}

# +++ Function "wrapper"

wrapper <- function(p, G, G_and_F ) {
  
  true_error(p[1], p[2],p[3], p[4], p[5],p[6],p[7], p[8],p[9],p[10], p[11],
             G, G_and_F )[[1]]
}


# +++ Function "TE_estimation"

TE_estimation <- function(G, G_and_F){
  
  sol <- optim(c(.125,.125,.125,.125,.125,.125,.125,.125,.2,.2,.2), wrapper, control=list(fnscale = 1),
               method="L-BFGS-B", 
               lower=c(0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,
                       0.00001, 0.00001, 0.00001),
               upper=c(1,1,1,1,1,1,1,1,.5,.5,.5),
               G=G, G_and_F=G_and_F)
  
  
  # Parameters estimated for each pattern
  parameters<-data.frame(sol$par[1:8]/sum(sol$par[1:8])) 
  colnames(parameters) <- c("Estimates")
  rownames(parameters) <- c("p111", "p112", "p121", "p122", "p211", "p212", "p221", "p222")
  
  # Errors estimated for each choice
  estimated_errors <- data.frame(sol$par[9:11])
  colnames(estimated_errors) <- c("Estimates")
  rownames(estimated_errors) <- c("e1", "e2", "e3")
  
  # Chi-squared statistic
  Chi2_test_statistic <- as.data.frame(sol$value) 
  colnames(Chi2_test_statistic) <- c("Chi2 statistic")
  
  # P value from Chi-squared table
  Chi2_table_p_value <- as.data.frame(1-pchisq(sol$value,df=5)) 
  colnames(Chi2_table_p_value) <- c("p value (DF=5)")
  
  
  # Contrast between observed and predicted proportions of responses
  G_exclusive <- G-G_and_F
  
  # Observed proportions
  obs_prop <- c(G_exclusive, G_and_F)/sum(G_and_F,G_exclusive ) 
  
  # Predicted proportions
  pred_prop <- true_error( sol$par[1], sol$par[2], sol$par[3],sol$par[4], sol$par[5],
                           sol$par[6], sol$par[7], sol$par[8],sol$par[9], sol$par[10],
                           sol$par[11],G=G, G_and_F= G_and_F)[[2]]
  
  contrast_prop<-cbind(obs_prop,pred_prop)
  rownames(contrast_prop) <- c("Only G 111", "Only G 112", "Only G 121", "Only G 122", 
                               "Only G 211", "Only G 212", "Only G 221", "Only G 222",
                               "Both (G and F) 111", "Both (G and F) 112", 
                               "Both (G and F) 121", "Both (G and F) 122", 
                               "Both (G and F) 211", "Both (G and F) 212", 
                               "Both (G and F) 221", "Both (G and F) 222")
  colnames(contrast_prop) <- c("Observed proportions", "Predicted proportions")
  
  # Contrast between observed and predicted frequencies of responses
  # Observed frequencies
  obs_freq <- c(G_exclusive, G_and_F) 
  size <- sum(obs_freq) #number of people
  
  # Predicted frequencies
  pred_freq <- pred_prop*size 
  contrast_freq <- cbind(obs_freq,pred_freq )
  rownames(contrast_freq) <- c("Only G 111", "Only G 112", "Only G 121", "Only G 122",
                               "Only G 211", "Only G 212", "Only G 221", "Only G 222",
                               "Both (G and F) 111", "Both (G and F) 112", 
                               "Both (G and F) 121", "Both (G and F) 122", 
                               "Both (G and F) 211", "Both (G and F) 212", 
                               "Both (G and F) 221", "Both (G and F) 222")
  colnames(contrast_freq) <- c("Observed frequencies", "Predicted frequencies")
  
  list(parameters=parameters, errors=estimated_errors, chisq=Chi2_test_statistic,
       p.value=Chi2_table_p_value,contrast.of.proportions=contrast_prop,
       contrast.of.frequencies=contrast_freq)
}


# +++ Function "TE_estimation_data"

TE_estimation_data <- function(data){
  
  data_list<-apply(data,1,list)
  number=nrow(data)
  
  results_list <-list() 
  for(i in 1:number) {
    G=data_list[[i]][[1]] [c(2:9)]
    G_and_F=data_list[[i]][[1]][c(10:17)]
    results_list[[i]]=TE_estimation(G,G_and_F)
  }
  
  results_list
  
}

# +++ Function "aggregate_estimates"

aggregate_estimates <- function (results, model) {
  n_elements <- length(results)
  
  parameters_data <- data.frame()
  
  if (model=="TE") {
    errors_data <- data.frame()
  }
  
  chisq_data <- data.frame()
  p.value_data <- data.frame()
  contrast_frequencies_data_Obs <- data.frame()
  contrast_frequencies_data_Pred <- data.frame()
  
  for (i in 1:n_elements) {
    
    parameters_results <- as.data.frame(t(results[[i]]$parameters))
    parameters_data <- rbind(parameters_data, parameters_results)
    
    if (model=="TE") {
      errors_results <- as.data.frame(t(results[[i]]$errors))
      errors_data <- rbind(errors_data, errors_results)
    }
    
    chisq_results <- as.data.frame(t(results[[i]]$chisq))
    chisq_data <- rbind(chisq_data, chisq_results)
    
    p.value_results <- as.data.frame(t(results[[i]]$p.value))
    p.value_data <- rbind(p.value_data, p.value_results)
    
    contrast_frequencies_results_Obs <- as.data.frame(t(results[[i]]$contrast.of.frequencies[,1]))
    contrast_frequencies_results_Pred <- as.data.frame(t(results[[i]]$contrast.of.frequencies[,2]))
    
    contrast_frequencies_data_Obs <- rbind(contrast_frequencies_data_Obs, 
                                           contrast_frequencies_results_Obs)
    contrast_frequencies_data_Pred <- rbind(contrast_frequencies_data_Pred, 
                                            contrast_frequencies_results_Pred)
    
    
  }
  
  if (model=="TE") {
    all_results <- cbind(contrast_frequencies_data_Obs, contrast_frequencies_data_Pred, 
                         parameters_data, errors_data, chisq_data, p.value_data)  
    
    names(all_results)[17:32] <- c( "Prediction Only G 111", "Prediction Only G 112", 
                                    "Prediction Only G 121", "Prediction Only G 122", 
                                    "Prediction Only G 211", "Prediction Only G 212",
                                    "Prediction Only G 221", "Prediction Only G 222",
                                    "Prediction Both (G and F) 111", 
                                    "Prediction Both (G and F) 112",
                                    "Prediction Both (G and F) 121",
                                    "Prediction Both (G and F) 122",
                                    "Prediction Both (G and F) 211",
                                    "Prediction Both (G and F) 212",
                                    "Prediction Both (G and F) 221",
                                    "Prediction Both (G and F) 222")
    names(all_results)[44:45] <- c("Chi-squared statistic", "P value (DF=5)")
    
  }
  
  if (model=="Independence") {
    all_results <- cbind(contrast_frequencies_data_Obs, contrast_frequencies_data_Pred, 
                         parameters_data, chisq_data, p.value_data)  
    
    names(all_results)[17:32] <- c( "Prediction Only G 111", "Prediction Only G 112", 
                                    "Prediction Only G 121", "Prediction Only G 122",
                                    "Prediction Only G 211", "Prediction Only G 212",
                                    "Prediction Only G 221", "Prediction Only G 222",
                                    "Prediction Both (G and F) 111", 
                                    "Prediction Both (G and F) 112",
                                    "Prediction Both (G and F) 121", 
                                    "Prediction Both (G and F) 122",
                                    "Prediction Both (G and F) 211", 
                                    "Prediction Both (G and F) 212",
                                    "Prediction Both (G and F) 221",
                                    "Prediction Both (G and F) 222")
    names(all_results)[36:37] <- c("Chi-squared statistic", "P value (DF=12)")
    
  }
  
  all_results
  
}




# FUNCTIONS TO PERFORM MONTE CARLO SIMULATIONS
# ----------------------------------------------


# +++ Function " MonteCarlo_TE"

MonteCarlo_TE <- function(number_samples, 
                          population_size, 
                          expected_frequencies_under_null, chisq.statistic) {
  
  # Random samples are generated according to the null of TE (i.e., predicted 
  # probabilities of the model). 
  # For the Re-Fitted Monte Carlo simulation: Each sample is then fit with TE, generating  
  # a set of 11 parameters (8 true probabilities and 3 errors) and the corresponding
  # Chi-squared statistics.
  # For the Conservative Monte Carlo simulation: Each sample is compared to the 
  # predicted frequencies using the original parameter estimates of the models.

  samples <- rmultinom(number_samples, population_size, expected_frequencies_under_null) 
  
  # In "samples", each column is a sample from a distribution (according to the null of 
  # TE). Each column has 16 values representing 2 sets of mutually exclusive 
  # frequencies: 8 values for patterns observed in G choices only; and other 8 values 
  # for patters observed in both, G and F, choices.
  
  # Re-Fitted Monte Carlo
  # ---------------------
  
  G_minus_both_montecarlo <- samples[1:8,]
  G_both_montecarlo <- samples[9:16,]
  G_montecarlo <- G_minus_both_montecarlo+G_both_montecarlo
  
  chi_all <- matrix(nrow=number_samples,ncol=1)
  sol_all <- matrix(nrow=number_samples,ncol=11)
  
  for(i in 1:number_samples) {
    sol <- optim(c(.125,.125,.125,.125,.125,.125,.125,.125,.2,.2,.2), wrapper, control=list(fnscale = 1),
                 method="L-BFGS-B", 
                 lower=c(0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,
                         0.00001, 0.00001, 0.00001, 0.00001),
                 upper=c(1,1,1,1,1,1,1,1,.5,.5,.5),
                 G=G_montecarlo[,i], G_and_F=G_both_montecarlo[,i])
    sol_all[i,] <- sol$par 
    chi_all[i,1] <- sol$value
    chi_values <- chi_all
  }
  
  parameters_all <- data.frame(sol_all[,1:8]/rowSums(sol_all[,1:8])) 
  estimated_errors_all <- sol_all[,9:11]
  colnames(parameters_all) <- c("p111", "p112", "p121", "p122", "p211", "p212", "p221",
                                "p222")
  colnames(estimated_errors_all) <- c("e1","e2","e3")
  colnames(chi_values) <- c("Chi squared statistics of Re-fitted MC")
  rownames(samples) <- c("Only G 111", "Only G 112", "Only G 121", "Only G 122",
                         "Only G 211", "Only G 212", "Only G 221", "Only G 222",
                         "Both (G and F) 111", "Both (G and F) 112", "Both (G and F) 121", 
                         "Both (G and F) 122", "Both (G and F) 211", "Both (G and F) 212", 
                         "Both (G and F) 221", "Both (G and F) 222")
  
  statistic <- chisq.statistic
  MonteCarlo_p_value <- (length(chi_values[chi_values>=as.numeric(statistic)   ]    ))/length(chi_values)
  
  # Conservative Monte Carlo
  # ------------------------
  
  
  p=expected_frequencies_under_null
  chi_all_conservative<-matrix(nrow=number_samples,ncol=1)
  
  for(i in 1:number_samples) {
    
    x=samples[,i]
    
    chi_sample <- (x - p)^2/p
    chi<-sum(chi_sample)
    chi_all_conservative[i,1]<-chi
  }
  
  
  p.value.conservative <- (length(chi_all_conservative[chi_all_conservative >= as.numeric(statistic) ]))/
                           length(chi_all_conservative)
  
  colnames(chi_all_conservative) <- c("Chi squared statistics of Conservative MC")
  
  
  
    list(parameters_Re.fitted.MC=parameters_all, errors_Re.fitted.MC= estimated_errors_all, 
       chisq_Re.fitted.MC=chi_values, p.value_Re.fitted.MC=MonteCarlo_p_value, 
       chisq_Conservative.MC=chi_all_conservative, p.value_Conservative.MC=p.value.conservative,
       samples=samples, original.chisq=statistic)
  
}



# +++ Function "MonteCarlo_simulation_TE"

MonteCarlo_simulation_TE <- function(G, G_and_F, number_samples) {
  
  TE_results <- TE_estimation(G=G, G_and_F=G_and_F)
  pred_freq <- TE_results$contrast.of.frequencies[,2]
  size <- sum(G)  
  statistic <- TE_results$chisq
  
  simulation_MonteCarlo_TE <- MonteCarlo_TE(number_samples=number_samples,
                                            population_size=size,
                                            expected_frequencies_under_null=pred_freq, 
                                            chisq.statistic = statistic)
  simulation_MonteCarlo_TE 
}


# +++ Function "TE_MonteCarlo_data"

TE_MonteCarlo_data <- function (data, number_samples) {
  data_list<-apply(data,1,list)
  number=nrow(data)
  
  results_list <-list() 
  for(i in 1:number) {
    G=data_list[[i]][[1]] [c(2:9)]
    G_and_F=data_list[[i]][[1]][c(10:17)]
    results_list[[i]]=MonteCarlo_simulation_TE(G=G, G_and_F=G_and_F, 
                                               number_samples=number_samples)
  }
  
  results_list
  
}

# +++ Function "aggregate_MonteCarlo_p.value"


aggregate_MonteCarlo_p.value <- function(results) {
  
  n_elements <- length(results)
  
  MonteCarlo.Re.fitted.p.value_data <- data.frame()
  MonteCarlo.Conservative.p.value_data <- data.frame()
  
  for (i in 1:n_elements) {

    MonteCarlo.Re.fitted.p.value_results <- as.data.frame(results[[i]]$p.value_Re.fitted.MC)
    MonteCarlo.Re.fitted.p.value_data <- rbind(MonteCarlo.Re.fitted.p.value_data, MonteCarlo.Re.fitted.p.value_results)
    
    MonteCarlo.Conservative.p.value_results <- as.data.frame(results[[i]]$p.value_Conservative.MC)
    MonteCarlo.Conservative.p.value_data <- rbind(MonteCarlo.Conservative.p.value_data, MonteCarlo.Conservative.p.value_results)
    
  }
  
  
  MonteCarlo.p.value_data_all<- cbind(MonteCarlo.Conservative.p.value_data, MonteCarlo.Re.fitted.p.value_data )
  
  colnames(MonteCarlo.p.value_data_all) <- c("Conservative MC - p value", "Re-fitted MC - p value")
  
  
  MonteCarlo.p.value_data_all
}





# FUNCTIONS TO PERFORM BOOTSTRAPPING SIMULATIONS TO OBTAIN CONFIDENCE INTERVALS
# ------------------------------------------------------------------------------


# +++ Function "CI_estimates_bootstrap"

CI_estimates_bootstrap <- function(obs_patterns_raw_data, d) {
  # obs_patterns_raw_data is a 1xNumber.of.participants vector in which each
  # cell represent an individual. 
  # Numbers on each cell will indicate the response pattern showed by the individual.
  # Numbers 1-8 represent response patterns 111, 112, 121, 122, 211, 212, 221 and 222, 
  # respectively, presented in G choices only (and not in F choices)
  # Numbers 9-16 represent response patterns presented in both, G and F, choices. 
  # "d" has indexes that will be used to build each sample.
  # "obs_sample" has the frequencies of each pattern in a sample. 
  # There are 8 cells for patterns observed in G choices only, and 8 cells for patterns 
  # observed in both, G and F, choices.
  
  samples <- obs_patterns_raw_data[d]
  obs_sample <- vector()
  obs_sample[1] <- sum(samples == 1)
  obs_sample[2] <- sum(samples == 2)
  obs_sample[3] <- sum(samples == 3)
  obs_sample[4] <- sum(samples == 4)
  obs_sample[5] <- sum(samples == 5)
  obs_sample[6] <- sum(samples == 6) 
  obs_sample[7] <- sum(samples == 7)
  obs_sample[8] <- sum(samples == 8)
  obs_sample[9] <- sum(samples == 9)
  obs_sample[10] <- sum(samples == 10)
  obs_sample[11] <- sum(samples == 11)
  obs_sample[12] <- sum(samples == 12)
  obs_sample[13] <- sum(samples == 13)
  obs_sample[14] <- sum(samples == 14)
  obs_sample[15] <- sum(samples == 15)
  obs_sample[16] <- sum(samples == 16)
  
  G_less_both_boot <- obs_sample[1:8]
  G_both_boot <- obs_sample[9:16]
  G_boot <- G_less_both_boot+G_both_boot
  
  
  sol <- optim(c(.125,.125,.125,.125,.125,.125,.125,.125,.2,.2,.2), wrapper, control=list(fnscale = 1),
               method="L-BFGS-B", 
               lower=c(0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,
                       0.00001,0.00001, 0.00001, 0.00001, 0.00001),
               upper=c(1,1,1,1,1,1,1,1,.5,.5,.5),
               G=G_boot, G_and_F=G_both_boot)
  sol <- sol$par 
  
  parameters <- sol[1:8]/sum(sol[1:8])
  estimated_errors<-sol[9:11]
  
  c(parameters, estimated_errors, obs_sample)
  
}


# +++ Function "CI_Bootstrap_TE"

CI_Bootstrap_TE <- function(G, G_and_F, number_samples) {
  
  G_exclusive <- G-G_and_F
  obs_freq <- c(G_exclusive, G_and_F) #observed frequencies
  names(obs_freq) <- c("Only G 111", "Only G 112", "Only G 121", "Only G 122", "Only G 211", 
                       "Only G 212", "Only G 221", "Only G 222", "Both (G and F) 111", 
                       "Both (G and F) 112", "Both (G and F) 121", "Both (G and F) 122", 
                       "Both (G and F) 211", "Both (G and F) 212", "Both (G and F) 221", 
                       "Both (G and F) 222")
  
  # Creating a vector in which each cell represents the responses of an individual
  obs_patterns_raw_data <- rep(1:16, obs_freq) 
  
  # Performing bootstrapping simulations
  result <- boot(obs_patterns_raw_data, CI_estimates_bootstrap, R=number_samples)
  
  # Estimating confidence intervals (percentile method)
  all_CI <- matrix(nrow=11,ncol=2)
  
  if (!is.null(boot.ci(result, type ="perc",index=1)[[4]][4:5])) { 
    all_CI[1,] <-  boot.ci(result, type ="perc",index=1)[[4]][4:5]}
  
  if (!is.null(boot.ci(result, type ="perc",index=2)[[4]][4:5])) { 
    all_CI[2,] <-  boot.ci(result, type ="perc",index=2)[[4]][4:5]}
  
  if (!is.null(boot.ci(result, type ="perc",index=3)[[4]][4:5])) { 
    all_CI[3,] <-  boot.ci(result, type ="perc",index=3)[[4]][4:5]}
  
  if (!is.null(boot.ci(result, type ="perc",index=4)[[4]][4:5])) { 
    all_CI[4,] <-  boot.ci(result, type ="perc",index=4)[[4]][4:5]}
  
  if (!is.null(boot.ci(result, type ="perc",index=5)[[4]][4:5])) { 
    all_CI[5,] <-  boot.ci(result, type ="perc",index=5)[[4]][4:5]}
  
  if (!is.null(boot.ci(result, type ="perc",index=6)[[4]][4:5])) { 
    all_CI[6,] <-  boot.ci(result, type ="perc",index=6)[[4]][4:5]}
  
  if (!is.null(boot.ci(result, type ="perc",index=7)[[4]][4:5])) { 
    all_CI[7,] <-  boot.ci(result, type ="perc",index=7)[[4]][4:5]}
  
  if (!is.null(boot.ci(result, type ="perc",index=8)[[4]][4:5])) { 
    all_CI[8,] <-  boot.ci(result, type ="perc",index=8)[[4]][4:5]}
  
  if (!is.null(boot.ci(result, type ="perc",index=9)[[4]][4:5])) { 
    all_CI[9,] <-  boot.ci(result, type ="perc",index=9)[[4]][4:5]}
  
  if (!is.null(boot.ci(result, type ="perc",index=10)[[4]][4:5])) { 
    all_CI[10,] <-  boot.ci(result, type ="perc",index=10)[[4]][4:5]}
  
  if (!is.null(boot.ci(result, type ="perc",index=11)[[4]][4:5])) { 
    all_CI[11,] <-  boot.ci(result, type ="perc",index=11)[[4]][4:5]}
  
  
  table_CI <- cbind(result$t0[1:11], all_CI)
  parnames = c("p111","p112","p121","p122","p211","p212","p221","p222","error1", 
               "error2", "error3")
  rownames(table_CI) = parnames
  columnnames = c("Estimates", "Lower bound", "Upper bound")
  colnames(table_CI) = columnnames
  
  boot_parameters <- result$t[,1:8]
  boot_errors <- result$t[,9:11]
  boot_samples <- result$t[,12:27]
  
  colnames(boot_parameters) <- c("p111", "p112", "p121", "p122", "p211", "p212", "p221",
                                 "p222")
  colnames(boot_errors) <- c("e1", "e2", "e3")
  colnames(boot_samples) <- c("Only G 111", "Only G 112", "Only G 121", "Only G 122",
                              "Only G 211", "Only G 212", "Only G 221", "Only G 222",
                              "Both (G and F) 111", "Both (G and F) 112", 
                              "Both (G and F) 121", "Both (G and F) 122", 
                              "Both (G and F) 211", "Both (G and F) 212",
                              "Both (G and F) 221", "Both (G and F) 222")
  
  
  list(boot.parameters=boot_parameters, boot.errors=boot_errors, boot.samples=boot_samples, 
       confidence.intervals=table_CI)
}

# +++ Function "TE_Bootstrapping.CI_data"


TE_Bootstrapping.CI_data <- function (data, number_samples) {
  data_list<-apply(data,1,list)
  number=nrow(data)
  
  results_list <-list() 
  for(i in 1:number) {
    G=data_list[[i]][[1]] [c(2:9)]
    G_and_F=data_list[[i]][[1]][c(10:17)]
    results_list[[i]]=CI_Bootstrap_TE(G=G, G_and_F=G_and_F, number_samples=number_samples)
  }
  results_list
  
}


# +++ Function "aggregate_boot.confidence.intervals"

aggregate_boot.confidence.intervals <- function (results) {
  results <- results.Boot 
  
  n_elements <- length(results)
  
  confidence.intervals_data <- data.frame()
  
  
  for (i in 1:n_elements) {
    
    confidence.intervals_results.estimates <- as.data.frame(t(results[[i]]$confidence.intervals[,1] ))
    confidence.intervals_results.lower.bound <- as.data.frame(t(results[[i]]$confidence.intervals[,2] ))
    confidence.intervals_results.upper.bound <- as.data.frame(t(results[[i]]$confidence.intervals[,3] ))
    confidence.intervals_results <- cbind( confidence.intervals_results.estimates, 
                                           confidence.intervals_results.lower.bound, 
                                           confidence.intervals_results.upper.bound)
    
    confidence.intervals_data<-rbind(confidence.intervals_data, confidence.intervals_results)
    
  }
  
  
  names(confidence.intervals_data)[12:22] <- c( "LowerBound p111", "LowerBound p112", 
                                                "LowerBound p121", "LowerBound p122", 
                                                "LowerBound p211", "LowerBound p212",
                                                "LowerBound p221", "LowerBound p222",
                                                "LowerBound error1", "LowerBound error2",
                                                "LowerBound error3")
  names(confidence.intervals_data)[23:33] <- c( "UpperBound p111", "UpperBound p112", 
                                                "UpperBound p121", "UpperBound p122", 
                                                "UpperBound p211", "UpperBound p212",
                                                "UpperBound p221", "UpperBound p222",
                                                "UpperBound error1", "UpperBound error2",
                                                "UpperBound error3")
  
  confidence.intervals_data

}



##########################################################################################
#
# FUNCTIONS APPLIED IN THE ESTIMATION OF INDEPENDENCE MODEL
#
# - To estimate Independence model parameters
#    +++ independence: sets the model equations
#    +++ wrapper_independence: sets the parameters to be used in the optimization algorithm
#    +++ Independence_estimation: estimates parameters using responses of an individual
#    +++ Independence_estimation_data: estimates parameters for all the individuals presented 
#                            in the researcher's data file
#    +++ aggregate_estimates: accumulates results obtained by "Independence_estimation_data"
#                             in a data frame.
#
# - To perform Monte Carlo simulations
#    +++ MonteCarlo_Independence: performs a number of Monte Carlo simulations given a vector of 
#                       predicted frequencies of the model (i.e., the null of the model that is
#                       used to create the samples)
#    +++ MonteCarlo_simulation_Independence: provides the vector of predicted frequencies to
#                                  function "MonteCarlo_Independence"
#    +++ Independence_MonteCarlo_data: computes simulations for all individuals presented 
#                            in the reseaercher's data file
#    +++ aggregate_MonteCarlo_p.value: accumulates all Monte Carlo p values estimated by
#                                      "Independence_MonteCarlo_data". 
#
###########################################################################################

# FUNCTIONS TO ESTIMATE INDEPENDENCE MODEL PARAMETERS
# -------------------------------------

# +++ Function "independence"

independence <- function( p1XX, pX1X,  pXX1, G, G_and_F ){
  
  G_exclusive <- G-G_and_F
  Total_G<-sum(G) #Number of participants
  
  # Defining probabilities of observing a particular pattern in G choices
  p111 <- p1XX*pX1X*pXX1
  p112 <- p1XX*pX1X*(1-pXX1)
  p121 <- p1XX*(1-pX1X)*pXX1
  p122 <- p1XX*(1-pX1X)*(1-pXX1)
  p211 <- (1-p1XX)*pX1X*pXX1
  p212 <- (1-p1XX)*pX1X*(1-pXX1)
  p221 <- (1-p1XX)*(1-pX1X)*pXX1
  p222 <- (1-p1XX)*(1-pX1X)*(1-pXX1)
  
  # Probability of observing a pattern in G choices
  prob_observe_G <- c(p111, p112, p121, p122, p211, p212, p221,p222)
  
  # Probability of observing a pattern in both G and F choices
  prob_observe_GF <- prob_observe_G^2
  
  # Probability of observing a pattern in G choices but not in F choices
  prob_observe_only_G <- prob_observe_G*(1-prob_observe_G)
  
  # Computing Chi-squared values
  Chi_GF_vector <- ((G_and_F-prob_observe_GF*Total_G)^2)/(prob_observe_GF*Total_G)
  Chi_GF<-sum( Chi_GF_vector)
  
  Chi_only_G_vector <- ((G_exclusive- prob_observe_only_G*Total_G)^2)/
    (prob_observe_only_G*Total_G)
  Chi_only_G<-sum(Chi_only_G_vector)
  
  Chi_total=Chi_GF+Chi_only_G
  
  list (Chi_total, c(prob_observe_only_G,prob_observe_GF))
}


# +++ Function "wrapper_independence"

wrapper_independence <- function(p, G, G_and_F ) {
  
  independence(p[1], p[2],p[3], G, G_and_F )[[1]]
}


# +++ Function "Independence_estimation"

Independence_estimation<-function(G, G_and_F) {
  
  
  sol <- optim(c(.5,.5,.5), wrapper_independence, control=list(fnscale = 1),
               method="L-BFGS-B", 
               lower=c(0.00001,0.00001,0.00001),
               upper=c(.99999,.99999,.99999),
               G=G, G_and_F=G_and_F)
  
  parameters<-data.frame(sol$par)
  colnames(parameters) <- c("Estimates")
  rownames(parameters) <- c("p1XX", "pX1X",  "pXX1")
  
  Chi2_test_statistic <- as.data.frame(sol$value) #Chi-squared statistic
  Chi2_table_p_value <- as.data.frame(1-pchisq(sol$value,df=12)) 
  colnames(Chi2_test_statistic) <- c("Chi2 statistic")
  colnames(Chi2_table_p_value) <- c("p value (DF=12)")
  
  G_exclusive <- G-G_and_F
  
  # Observed proportions
  obs_prop <- c(G_exclusive, G_and_F)/sum(G_and_F,G_exclusive ) 
  
  # Predicted proportions
  pred_prop <- independence( sol$par[1], sol$par[2], sol$par[3],G=G, 
                             G_and_F= G_and_F)[[2]] 
  
  contrast_prop<-cbind(obs_prop,pred_prop)
  rownames(contrast_prop) <- c("Only G 111", "Only G 112", "Only G 121", "Only G 122", 
                               "Only G 211", "Only G 212", "Only G 221", "Only G 222",
                               "Both (G and F) 111", "Both (G and F) 112", 
                               "Both (G and F) 121", "Both (G and F) 122", 
                               "Both (G and F) 211", "Both (G and F) 212",
                               "Both (G and F) 221", "Both (G and F) 222")
  colnames(contrast_prop) <- c("Observed proportions", "Predicted proportions")
  
  # Observed frequencies
  obs_freq <- c(G_exclusive, G_and_F) 
  size <- sum(obs_freq) #number of people
  
  # Predicted frequencies
  pred_freq <- pred_prop*size 
  contrast_freq <- cbind(obs_freq,pred_freq )
  rownames(contrast_freq) <- c("Only G 111", "Only G 112", "Only G 121", "Only G 122", 
                             "Only G 211", "Only G 212", "Only G 221", "Only G 222",
                             "Both (G and F) 111", "Both (G and F) 112", 
                             "Both (G and F) 121", "Both (G and F) 122", 
                             "Both (G and F) 211", "Both (G and F) 212", 
                             "Both (G and F) 221", "Both (G and F) 222")
  colnames(contrast_freq) <- c("Observed frequencies", "Predicted frequencies")
  
  list(parameters=parameters, chisq=Chi2_test_statistic, p.value=Chi2_table_p_value,
       contrast.of.proportions=contrast_prop, contrast.of.frequencies=contrast_freq)
}


# +++ Function "Independence_estimation_data"

Independence_estimation_data <- function(data){
  
  data_list <- apply(data,1,list)
  number=nrow(data)
  
  results_list <-list() 
  for(i in 1:number) {
    G=data_list[[i]][[1]] [c(2:9)]
    G_and_F=data_list[[i]][[1]][c(10:17)]
    results_list[[i]]=Independence_estimation(G,G_and_F)
  }
  
  results_list
  
}

# +++ Function "aggregate_estimates" (already defined in functions for TE estimations)


# FUNCTIONS TO PERFORM MONTE CARLO SIMULATIONS
# ----------------------------------------------


# +++ Function "MonteCarlo_Independence"

MonteCarlo_Independence <- function(number_samples, 
                                    population_size, 
                                    expected_frequencies_under_null, chisq.statistic) {
  
  samples <- rmultinom(number_samples, population_size, expected_frequencies_under_null) 
  
  # Re-fitted Monte Carlo
  # ---------------------
  G_minus_both_montecarlo <- samples[1:8,]
  G_both_montecarlo <- samples[9:16,]
  G_montecarlo <- G_minus_both_montecarlo+G_both_montecarlo
  
  chi_all <- matrix(nrow=number_samples,ncol=1)
  sol_all <- matrix(nrow=number_samples,ncol=3)
  
  for(i in 1:number_samples) {
    
    sol <- optim(c(.5,.5,.5), wrapper_independence, control=list(fnscale = 1),
                 method="L-BFGS-B", 
                 lower=c(0.00001,0.00001,0.00001),
                 upper=c(.99999,.99999,.99999),
                 G=G_montecarlo[,i], G_and_F=G_both_montecarlo[,i])  
    
    sol_all[i,] <- sol$par 
    chi_all[i,1] <- sol$value
    chi_values <- chi_all
  }
  
  parameters_all <- data.frame(sol_all)
  colnames(parameters_all) <- c("p1xx", "px1x", "pxx1")
  colnames(chi_values) <- c("Chi squared statistics of Re-fitted MC")
  rownames(samples) <- c("Only G 111", "Only G 112", "Only G 121", "Only G 122", 
                         "Only G 211", "Only G 212", "Only G 221", "Only G 222",
                         "Both (G and F) 111", "Both (G and F) 112",
                         "Both (G and F) 121", "Both (G and F) 122", 
                         "Both (G and F) 211", "Both (G and F) 212",
                         "Both (G and F) 221", "Both (G and F) 222")
  
  
  statistic <-chisq.statistic
  MonteCarlo_p_value <- (length(chi_values[chi_values>=as.numeric(statistic)]))/
    length(chi_values)
  
  
  
  
  # Conservative Monte Carlo
  # ------------------------
  
  
  p=expected_frequencies_under_null
  chi_all_conservative<-matrix(nrow=number_samples,ncol=1)
  
  for(i in 1:number_samples) {
    
    x=samples[,i]
    
    chi_sample <- (x - p)^2/p
    chi<-sum(chi_sample)
    chi_all_conservative[i,1]<-chi
  }
  
  
  p.value.conservative <- (length(chi_all_conservative[chi_all_conservative >= as.numeric(statistic) ]))/
    length(chi_all_conservative)
  
  colnames(chi_all_conservative) <- c("Chi squared statistics of Conservative MC")
  
  
  
  list(parameters_Re.fitted.MC=parameters_all,  
       chisq_Re.fitted.MC=chi_values, p.value_Re.fitted.MC=MonteCarlo_p_value, 
       chisq_Conservative.MC=chi_all_conservative, p.value_Conservative.MC=p.value.conservative,
       samples=samples, original.chisq=statistic)
  
}


# +++ Function "MonteCarlo_simulation_Independence"

MonteCarlo_simulation_Independence <- function(G, G_and_F, number_samples) {
  
  Independence_results <- Independence_estimation(G=G, G_and_F=G_and_F)
  pred_freq<-Independence_results$contrast.of.frequencies[,2]
  size <- sum(G)
  
  
  statistic<-Independence_results$chisq
  
  simulation_MonteCarlo_Independence <- MonteCarlo_Independence(number_samples=number_samples,
                                                                population_size=size,
                                                                expected_frequencies_under_null=pred_freq, 
                                                                chisq.statistic = statistic)
  
  simulation_MonteCarlo_Independence
}


# +++ Function "Independence_MonteCarlo_data"

Independence_MonteCarlo_data <- function (data, number_samples) {
  data_list<-apply(data,1,list)
  number=nrow(data)
  
  results_list <-list() 
  for(i in 1:number) {
    G=data_list[[i]][[1]] [c(2:9)]
    G_and_F=data_list[[i]][[1]][c(10:17)]
    results_list[[i]]=MonteCarlo_simulation_Independence(G=G, G_and_F=G_and_F, 
                                                         number_samples=number_samples)
  }
  
  results_list
  
}

# +++ Function "aggregate_MonteCarlo_p.value" (already defined in functions for TE estimations)
#################################################################################################
#
#################################################################################################
#
#
#############################################################################################
#       EXAMPLE SET-UP:  THE USER CAN REVISE THIS SECTION
#       This section illustrates use of the functions
#       (Additional Documentation of the main functions at the top of this listing)
#############################################################################################
data_file<-"example.txt"            # Input data are in a tab-delimited text file with header
output_file<-"output_TE.txt"        #  Output file is a text file
number_runs<-1000                  # Change the number of runs to 10000 for more accuracy
seed_number<- 18743                 # This random number seed can be changed.
#                                     Keep the seed number fixed to reproduce results
#   
#                  READ IN THE DATA
response_patterns <- read.table(file = data_file, header=T)
# The file "example.txt" contains the data.
# The operator " <- " assign values to variables.
#   The next line puts back out the first few lines of the input, so you can check
head(response_patterns)
#                  
#                  ESTIMATE PARAMETERS AND FIT THE MODELS
# The next lines fit the TE and Independence models to the data and store results
results_TE <- TE_estimation_data(data=response_patterns)
results_Independence <-Independence_estimation_data(data=response_patterns)
#
label_TE=c("True and Error Model Results")
label_Indep=c("Independence Model Results")
#   Output is directed to the file designated below
sink(output_file,append=TRUE)
print(label_TE)
print(results_TE)
print(label_Indep)
print(results_Independence)
sink()
#
#  Here results of TE and Indep Models put in tables and saved as CSV files
all_estimates_TE <-aggregate_estimates(results=results_TE, model="TE")
write.csv(all_estimates_TE, file="parameter estimates of TE.csv")
#   The first argument of write.csv() should be a data frame (i.e., a set of
#   vectors of equal length, resembling a table).
#   Thus, a list, such as "results_TE" or "results_Independence" needs to be
#   converted into a data frame before being exported into a CSV file.
all_estimates_Independence <- aggregate_estimates(results=results_Independence, model="Independence")
write.csv(all_estimates_Independence, file="parameter estimates of Independence.csv")
#
#                MONTE CARLO SECTION
#
#   Observe that we declare the seed for the random generator process prior
#   to the simulations in order to create results that can be reproduced,
#   which is done with command set.seed(<seed>)".
#
#     Monte Carlo of TE Model
set.seed(seed_number)   # We can use any other seed number for other pseudo random samples                                                                                                                                                                                                                                                       
results_TE_MC <- TE_MonteCarlo_data(data=response_patterns,number_samples = number_runs )
#   It may take 20 min. per case to simulate and analyze 10,000 samples via MC
#   Here is a printout of the first few samples, parameters for case [1]
sink(output_file,append=TRUE)
print("First few samples parameters for Case 1")
print(head(results_TE_MC[[1]]$parameters))
sink()
all_p.values_MC_TE <- aggregate_MonteCarlo_p.value(results=results_TE_MC)
#
#       Monte Carlo of Independence Model
set.seed(seed_number)  # We could use any other seed number  
results_Independence_MC <- Independence_MonteCarlo_data(data=response_patterns,number_samples = number_runs)
#   Aggregating the p values of all of our participants in data frame "all_p.values_MC_Independence".
all_p.values_MC_Independence <-aggregate_MonteCarlo_p.value(results=results_Independence_MC)
#
#       Write out the p-values into files for both Models
write.csv(all_p.values_MC_TE, file= "Monte Carlo p values TE.csv")
write.csv(all_p.values_MC_Independence, file="Monte Carlo p values Ind.csv")
#
#      Draw Histograms of the Fit of the TE Model
library(scales)
#   Draw Histogram for Case #1 TE Model Refit to each sample
pdf("MonteCarlo_TE_Refit_Case1.pdf")
MonteCarlo_statistic_histiogram <- hist(results_TE_MC[[1]]$chisq_Re.fitted.MC,col=scales::alpha('skyblue',.5),breaks=50,xlim = c(0, 30),border=F,main="Monte Carlo TE Model",xlab="Chi2 Re-Fit in each Sample")
abline(v=results_TE_MC[[1]]$original.chisq,col="red")
dev.off()
pdf("MonteCarlo_TE_Conservative_Case1.pdf")
MonteCarlo_statistic_histiogram <- hist(results_TE_MC[[1]]$chisq_Conservative.MC,col=scales::alpha('skyblue',.5),breaks=50,xlim = c(0, 30),border=F,main="Monte Carlo TE Model (Conservative)",xlab="Chi2 Conservative")
abline(v=results_TE_MC[[1]]$original.chisq,col="red")
dev.off()
#
#   Draw Histograms for Case #1 Independence Model
pdf("MonteCarlo_Independence_Refit_Case1.pdf")
MonteCarlo_statistic_histiogram <- hist(results_Independence_MC[[1]]$chisq_Re.fitted.MC,col=scales::alpha('skyblue',.5),breaks=50,xlim = c(0, 30),border=F,main="Monte Carlo Independence Model (Re-fit in each Sample)",xlab="Chi2 Re-Fit in each Sample")
abline(v=results_Independence_MC[[1]]$original.chisq,col="red")
dev.off()
#
pdf("MonteCarlo_Independence_Conservative_Case1.pdf")
MonteCarlo_statistic_histiogram <- hist(results_Independence_MC[[1]]$chisq_Conservative.MC,col=scales::alpha('skyblue',.5),breaks=50,xlim = c(0, 30),border=F,main="Monte Carlo Independence Model (Conservative)",xlab="Chi2 Fit Original Data")
abline(v=results_Independence_MC[[1]]$original.chisq,col="red")
dev.off()
#
##############################################################################
#                  BOOTSTRAPPING SECTION
#
library(boot)
set.seed(seed_number) # Sets the random number seed so you get the same results as we get. 
results.Boot <- TE_Bootstrapping.CI_data (data=response_patterns, number_samples=number_runs) 
all_boot_CI <- aggregate_boot.confidence.intervals(results=results.Boot)
write.csv(all_boot_CI, file="Bootstrapped confidence intervals TE.csv") 
#
#  The next section shows how to graph the results of the bootstrapping distributions
#  Here smoothed density distributions are drawn.  One could use histograms instead.
#
#            SAVE BOOTSTRAPPED DENSITY GRAPHS  
pdf("Bootstrapped Density of p111.pdf")
plot(density(results.Boot[[1]]$boot.parameters[,c("p111")]),xlim = c(0, 1))
dev.off()
pdf("Bootstrapped Density of p112.pdf")
plot(density(results.Boot[[1]]$boot.parameters[,c("p112")]),xlim = c(0, 1))
dev.off()
pdf("Bootstrapped Density of p121.pdf")
plot(density(results.Boot[[1]]$boot.parameters[,c("p121")]),xlim = c(0, 1))
dev.off()
pdf("Bootstrapped Density of p122.pdf")
plot(density(results.Boot[[1]]$boot.parameters[,c("p122")]),xlim = c(0, 1))
dev.off()
pdf("Bootstrapped Density of p211.pdf")
plot(density(results.Boot[[1]]$boot.parameters[,c("p211")]),xlim = c(0, 1))
dev.off()
pdf("Bootstrapped Density of p212.pdf")
plot(density(results.Boot[[1]]$boot.parameters[,c("p212")]),xlim = c(0, 1))
dev.off()
pdf("Bootstrapped Density of p221.pdf")
plot(density(results.Boot[[1]]$boot.parameters[,c("p221")]),xlim = c(0, 1))
dev.off()
pdf("Bootstrapped Density of p222.pdf")
plot(density(results.Boot[[1]]$boot.parameters[,c("p222")]),xlim = c(0, 1))
dev.off()
pdf("Bootstrapped Density of e1.pdf")
plot(density(results.Boot[[1]]$boot.errors[,c("e1")]),xlim = c(0, 1))
dev.off()
pdf("Bootstrapped Density of e2.pdf")
plot(density(results.Boot[[1]]$boot.errors[,c("e2")]),xlim = c(0, 1))
dev.off()
pdf("Bootstrapped Density of e3.pdf")
plot(density(results.Boot[[1]]$boot.errors[,c("e3")]),xlim = c(0, 1))
dev.off()
#
##########################################################################################
#   This example is set up to save only selected results to illustrate the program
#   For example, it saves only graphs for the first case.
#   You can easily revise it to save other information and make other graphs  
#   See documentation section.
##########################################################################################
