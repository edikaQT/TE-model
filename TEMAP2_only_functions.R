
#########################################################################################
#
#      DOCUMENTATION OF THE MAJOR FUNCTIONS (useful to those who program)
#      (skip to line 2325, where a non-programmer can adjust the inputs)
#
#########################################################################################
#
#                       _______________________________________
#
#                                EXPERIMENTAL DESIGN
#                       _______________________________________
#
# We have two choices (1) R vs S and (2) R' vs S', each choice with two replications.
#  
# Labelling of preferences
# --------------------------
# In the Choice between R and S (or R' and S'):
# Preference for R (or R') is coded as 0
# Preference for S (or S') is coded as 1.
#
# True Preferences
# ------------------
# Then, 4 patterns of true preferences are possible:
# a_00: true preference for R in Choice 1 and R' in Choice 2
# a_01: true preference for R in Choice 1 and S' in Choice 2
# a_10: true preference for S in Choice 1 and R' in Choice 2
# a_11: true preference for S in Choice 1 and S' in Choice 2
#  
# Errors in a choice
# ---------------------
# When you truly prefer R (or R') in a choice: 
#    error in choice (1) is labelled as "error_R1"
#    error in choice (2) is labelled as "error_R2"
#  
# When you truly prefer S (or S') in a choice: 
#    error in choice (1) is labelled as "error_S1"
#    error in choice (2) is labelled as "error_S2"
#  
# Patterns of responses
# ------------------------------
# Then, we have 16 pattern of responses PXYZW: P0000,  P0001,P0010,P0011, etc. 
# Where "XY" are responses in Replication 1 and "ZW" are responses in replication 2.
# For instance, P1001 is the pattern:
# [Replication 1, Choice1=S][Replication1, Choice2=R'][Replication 2, Choice1=R][Replication2, Choice2=S']
#
#
#     __________________________________________________________________
#
#             MAIN FUCTIONS FOR TRUE AND ERROR MODEL ESTIMATION
#     __________________________________________________________________
#
# True and Error model estimation
# 
# Function 'TE_estimation_data' fits the TE model to a dataset by 
# minimizing a chi-squared statistic. It estimates the probabilities 
# of four true preference patterns denoted by 
# a_00, a_01,a_10, a_11, 
# along with four errors error_R1, error_R2, error_S1, error_S2
# for two choice problems with two replications. 
# The estimation assumes strict positive probabilities and errors, which 
# are constrained to be no greater than a half. Syntax is as follows:
# 
# TE_estimation_data(data,restrictions, equal.prob.of.error.in.choice)
#
# Where "data" is a data frame with 17 columns describing the ID of participants 
# and the 16 possible frequencies of patterns "P0000", "P0001","P0010",
# "P0011", "P0100","P0101","P0110", "P0111", "P1000","P1001","P1010",
# "P1011", "P1100","P1101","P1110","P1111".
# Each row contains data for a different individual. 
#
# Where "restrictions" is a 4-element vector that indicates which if any
# probability of true preference is set to a fixed value.
# One or two parameters of [a_00, a_01, a_10, a_11] can be fixed. 
# For instance, to fix a_00 to .3, and a_11 to .5, set:
#
# user_restrictions <- c(.3, "FREE", .5, "FREE") 
#
# Or to let all parameters to be free: 
#
# user_restrictions <- c("FREE", "FREE", "FREE", "FREE") 
#
# There are 4 error terms in the model: "error_R1", error_R2", "error_S1", "error_S2".
# The user can constraint the errors by specifying:
#
# equal.prob.of.error.in.choice = "TE-4"
# (to estimate 4 error terms)
# equal.prob.of.error.in.choice = "TE-2"
# (to estimate 2 error terms, assuming errors in a given choice are the same, i.e., error_R1=error_S1 & error_R2=error_S2)
# equal.prob.of.error.in.choice = "TE-1"
# (to estimate 1 error term, assuming errors in any choice and replication are the same, i.e., error_R1=error_S1=error_R2=error_S2)
#
# The program generates the following components: 
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
# ... [Case Number]$matrix.contrast.frequencies
# in which, the "..." is the name of the list storing the results.
# 
#
# Monte Carlo simulation of True and Error model
# ----------------------------------------------
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
# TE_MonteCarlo_data(data, restrictions, equal.prob.of.error.in.choice, number_samples) 

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
# ----------------------------------------------
# 
# Function "TE_Bootstrapping.CI_data" performs simulations of the TE model 
# by fitting random samples resampled from the original data with replacement. 
# Parameters are estimated for each simulated sample and are then used to 
# estimate 95% confidence intervals using the percentile method. Syntax:
# 
# TE_Bootstrapping.CI_data (data, number_samples, restrictions, equal.prob.of.error.in.choice)
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
# ----------------------------------------------
# 
# Function "Independence_estimation_data" computes the chi-squared index for a model
# that assumes that the predicted probability of observing a response pattern is given
# by the product of the marginal choice probabilities.
# For instance, pattern P1001 has the following probability:
# P1001 = p(S)p(R')P(R)p(S')
# where the first two components, p(S) and p(R'), are the marginal choice probabilities for Replication 1
# of the choices R vs S and R' vs S',
# and the last two components, P(R) and p(S'), are are the marginal choice probabilities for Replication 2
# of the same choices.
#
# The function generates a list with following components: 
# marginal choice probability estimates, a Chi-squared statistic, 
# observed and predicted proportions of response patterns, and observed and predicted frequencies 
# of response patterns. These components can be obtained by typing:
# 
# ... [Case Number]$Marginal.prob.Independence
# ... [Case Number]$chisq
# ... [Case Number]$contrast.of.frequencies
# ... [Case Number]$matrix.contrast.frequencies
# in which, the "..." is the name of the list storing the results.
#
# Tests of Chi2 differences to compare Restricted and Unrestricted models
# -----------------------------------------------------------------
# 
# Function "test_chi2_diff_data" computes the standard test of differences in chi-squared values
# for a restricted model (in which the user fix some parameters) and an unrestricted model (in which
# all parameters are FREE).
# Function "MonteCarlo_test_chi2_diff_data" computes Monte Carlo simulations for the same test.
# samples are created under the null of the Restricted model.
#
# If the chi2-diff-value (from Chi-squared table or from Monte Carlo simulations) is significant (p<.05), the 
# Unrestricted model fits the data better than the Restricted model in which a number of parameters are 
# fixed to certain value.
# Otherwise, both models fit equally well and we can accept the Restricted model and fix some parameters.
#
# The function generates a list with following components: 
# Standard Chi-squared difference test and Monte Carlo Re-fitted Chi-squared difference test.
# These components can be obtained by typing:
# 
# ... [Case Number]$[TE Model]$Chi2_diff_standard
# ... [Case Number]$[TE Model]$Chi2_diff_Re.fitted.MC
#
# in which, the "..." is the name of the list storing the results, and the TE Model is either TE-4,
# TE-2 or TE-1 (e.g., ...[Case Number]$TE_1_Chi2_diff_test$Chi2_diff_standard)
#
#######################################################################################
#
#
##########################################################################################
#
# FUNCTIONS APPLIED IN THE ESTIMATION OF TE MODEL
#
# - To estimate TE parameters
#    +++ true_error: sets the model equations
#    +++ wrapper: sets the parameters to be used in the optimization algorithm
#    +++ optim_TE: set the algorithm for optimization
#    +++ optim_TE_N_times: to avoid local minima, the optimization is solved "N" times
#                          from the best fitted value in time "N-1". 
#                          Currently N is set to 3 times.
#    +++ TE_estimation: estimates parameters using responses of an individual
#    +++ TE_estimation_data: estimates parameters for all the individuals presented 
#                            in the researcher's data file
#    +++ aggregate_estimates: accumulates results obtained by "TE_estimation_data"
#                             in a data frame
#
# - To perform Monte Carlo simulations
#    +++ extract_results: extract chi_values, parameters and errors from each simulation
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
#
#///////////////////////////////////////////////////////////////////////////////////////
# ___________________________________________________
#
#    FUNCTIONS TO ESTIMATE THE TRUE AND ERROR MODEL
# ___________________________________________________
#
#///////////////////////////////////////////////////////////////////////////////////////

# +++ Function "true_error"



true_error <- function(prob, responses, restrictions, fit_method){
  
  # "prob" is an 8-element vector containing the parameters of the model
  # "responses" has the participant choice patterns
  # "restrictions" is a 4-element vector indicating which parameters have fixed
  # values
  
  # We ensure that the parameters used by the model are strictly positive
  # Parameters for "true" probability of preferences
  a_00<-ifelse(prob[1]<0.000001,0.000001,prob[1])
  a_01<-ifelse(prob[2]<0.000001,0.000001,prob[2])
  a_10<-ifelse(prob[3]<0.000001,0.000001,prob[3])
  a_11<-ifelse(prob[4]<0.000001,0.000001,prob[4])
  # Parameters for errors
  error_R1 <- ifelse(prob[5]<0.000001,0.000001,prob[5])
  error_R2 <- ifelse(prob[6]<0.000001,0.000001,prob[6])
  error_S1 <- ifelse(prob[7]<0.000001,0.000001,prob[7])
  error_S2 <- ifelse(prob[8]<0.000001,0.000001,prob[8])
  
  
  # "fixed" indicates which "true" probability preferences are restricted to certain values
  fixed <- as.numeric(ifelse(restrictions=="FREE", 0, 1))
  
  parameters <- c(a_00, a_01,a_10, a_11)
  
  # Adjusting "true" probabilities to sum 1
  parameters_logical <-data.frame(cbind(logical=fixed, parameters))
  parameters_logical$order.names <- c(1:4 )
  
  parameters.fixed <- parameters_logical[which(parameters_logical$logical==1),2:3] 
  parameters.no.fixed <- parameters_logical[which(parameters_logical$logical!=1),2:3]
  
  sum.parameters.fixed    <- sum(parameters.fixed$parameters) 
  sum.parameters.no.fixed <- sum(parameters.no.fixed$parameters)
  
  true.prob.parameters.no.fixed <- (parameters.no.fixed$parameters / sum.parameters.no.fixed)*(1-sum.parameters.fixed)
  
  true_prob <- data.frame(cbind(c(parameters.fixed$parameters, true.prob.parameters.no.fixed ), 
                                order=c(parameters.fixed$order.names, parameters.no.fixed$order.names)))
  
  true_prob <- true_prob[order(true_prob$order),] [,1] 
  
  # Number of participants
  Total_responses <- sum(responses) 
  
  # Creating a matrix with all possible combinations of errors in two choices 
  # with two replications
  
  V1 <- c(1-error_R1,	1-error_R1,	1-error_R1,	1-error_R1,
          1-error_R1,	1-error_R1,	1-error_R1,	1-error_R1,
          error_R1,		error_R1,		error_R1,		error_R1,
          error_R1,		error_R1,		error_R1,		error_R1)
  
  V2 <- c(1-error_R2,	1-error_R2,	1-error_R2,	1-error_R2,
          error_R2,		error_R2,		error_R2,		error_R2,
          1-error_R2,	1-error_R2,	1-error_R2,	1-error_R2,
          error_R2,		error_R2,		error_R2,		error_R2)
  
  V3 <- c(1-error_R1,	1-error_R1,		error_R1,		error_R1,
          1-error_R1,	1-error_R1,		error_R1,		error_R1,
          1-error_R1,	1-error_R1,		error_R1,		error_R1,
          1-error_R1,	1-error_R1,		error_R1,		error_R1)
  
  V4 <- c(1-error_R2,	error_R2,		1-error_R2,		error_R2,
          1-error_R2,	error_R2,		1-error_R2,		error_R2,
          1-error_R2,	error_R2,		1-error_R2,		error_R2,
          1-error_R2,	error_R2,		1-error_R2,		error_R2)
  
  V5 <- c(error_S2,		error_S2,		error_S2,		error_S2,
          1-error_S2,	1-error_S2,	1-error_S2,	1-error_S2,
          error_S2,		error_S2,		error_S2,		error_S2,
          1-error_S2,	1-error_S2,	1-error_S2,	1-error_S2)
  
  V6 <- c(error_S2,		1-error_S2,	error_S2,		1-error_S2,
          error_S2,		1-error_S2,	error_S2,		1-error_S2,
          error_S2,		1-error_S2,	error_S2,		1-error_S2,
          error_S2,		1-error_S2,	error_S2,		1-error_S2)
  
  V7 <- c(error_S1,		error_S1,		error_S1,		  error_S1,
          error_S1,		error_S1,		error_S1,		  error_S1,
          1-error_S1,	1-error_S1,	1-error_S1,	  1-error_S1,
          1-error_S1,	1-error_S1,	1-error_S1,	  1-error_S1)
  
  V8 <- c(error_S1,		error_S1,		1-error_S1,		1-error_S1,
          error_S1,		error_S1,		1-error_S1,		1-error_S1,
          error_S1,		error_S1,		1-error_S1,		1-error_S1,
          error_S1,		error_S1,		1-error_S1,		1-error_S1)
  
  # sum(a_00* V1*V2*V3*V4 + a_01*V1*V5*V3*V6 + a_10*V7*V2*V8*V4 + a_11*V7*V5*V8*V6)
  # sum to 1
  a_00<- true_prob[1]
  a_00<- true_prob[1]
  a_01<- true_prob[2]
  a_10<- true_prob[3]
  a_11<- true_prob[4]
  
  prob.responses  <-c(a_00* V1*V2*V3*V4 + a_01*V1*V5*V3*V6 +
                        a_10*V7*V2*V8*V4 + a_11*V7*V5*V8*V6)
  
  
  if (fit_method=="CHI2"){
    
    # Computing Chi-squared values
    
    Chi_vector <- ((responses - prob.responses*Total_responses)^2)/
      (prob.responses*Total_responses)
    Chi_total <- sum(Chi_vector)
  } else {
    
    # Computing the G2 index or deviance statistic (default method)
    
    G_data <-  as.data.table(cbind(responses , G= log(responses / (prob.responses*Total_responses) ) ))
    G_data[responses==0, G:=0 ]
    G_data[, G_vector:=responses*G]  
    G_index <- 2*sum(G_data$G_vector)
    Chi_total <- G_index
  }
  list (Chi_total,  prob.responses , true_prob) #errors never are an output here
  #we call "Chi_total" to the index of fit selected by the user, which can be  either the Chi-square or the G2 index
}

#///////////////////////////////////////////////////////////////////////////////////////



# +++ Function "wrapper"

wrapper <- function(p, responses, restrictions, fit_method  ) {
  
  if (length(p)==5) {
    true_error(prob=c(p[1], p[2],p[3], p[4], p[5],p[5],p[5], p[5]),
               responses=responses, restrictions=restrictions, fit_method=fit_method )[[1]]
  } else if (length(p)==6){
    true_error(prob=c(p[1], p[2],p[3], p[4], p[5],p[6],p[5], p[6]),
               responses=responses, restrictions=restrictions, fit_method=fit_method )[[1]]
  } else {
    true_error(prob=c(p[1], p[2],p[3], p[4], p[5],p[6],p[7], p[8]),
               responses=responses, restrictions=restrictions, fit_method=fit_method )[[1]]			 
  }
}

#///////////////////////////////////////////////////////////////////////////////////////

# +++ Function "optim_TE"

optim_TE <- function(initial_values, lower_values, upper_values, 
                     responses, restrictions, equal.prob.of.error.in.choice, fit_method) {
  
  sol <- optim(c(initial_values), 
               wrapper, control=list(fnscale = 1),
               method="L-BFGS-B", 
               lower=lower_values,
               upper=upper_values,
               responses=responses, 
               restrictions=restrictions, 
               fit_method=fit_method)
  
  
  if (equal.prob.of.error.in.choice=="TE-1") {
    
    parameters <- data.frame(true_error( prob=c(sol$par[1], sol$par[2], sol$par[3],sol$par[4], sol$par[5],
                                                sol$par[5], sol$par[5], sol$par[5]),
                                         responses=responses, restrictions=restrictions, fit_method=fit_method)[[3]])
    #these are the true probabilitites, 4 values
    
    estimated_errors <- data.frame(sol$par[5])
    initial_values <-c(parameters[,1],estimated_errors[,1])									 
    
  } else  if (equal.prob.of.error.in.choice=="TE-2") {
    
    
    #parameters: a_00, a_01,   a_10, a_11, error_R1 , error_R2,  error_S1, error_S2 
    #In TE2, error_R1=error_S1 & error_R2=error_S2
    
    parameters <- data.frame(true_error( prob=c(sol$par[1], sol$par[2], sol$par[3],sol$par[4], sol$par[5],
                                                sol$par[6], sol$par[5], sol$par[6]),
                                         responses=responses, restrictions=restrictions, fit_method=fit_method)[[3]])
    
    estimated_errors <- data.frame(sol$par[5:6])
    initial_values <-c(parameters[,1],estimated_errors[,1])									 
    
  } else { 
    parameters <- data.frame(true_error(prob= c(sol$par[1], sol$par[2], sol$par[3],sol$par[4], sol$par[5],
                                                sol$par[6], sol$par[7], sol$par[8]),
                                        responses=responses, restrictions=restrictions, fit_method = fit_method)[[3]])
    
    estimated_errors <- data.frame(sol$par[5:8])
    initial_values <-c(parameters[,1],estimated_errors[,1])
  }
  list(parameters, estimated_errors, initial_values, sol) 
}

#///////////////////////////////////////////////////////////////////////////////////////

# +++ Function "optim_TE_N_times"

optim_TE_N_times <- function (initial_values, 
                              lower_values, 
                              upper_values, 
                              responses, 
                              restrictions,
                              equal.prob.of.error.in.choice, 
                              fit_method,
                              N ) {
  for(i in 1:N) {
    
    R <- optim_TE(initial_values=initial_values, 
                  lower_values=lower_values, 
                  upper_values=upper_values, 
                  responses=responses, 
                  restrictions=restrictions, 
                  equal.prob.of.error.in.choice=equal.prob.of.error.in.choice,
                  fit_method=fit_method)
    
    initial_values <-R[[3]]
  }
  R # contains a list of parameters, estimated_errors, initial_values and the optimization solution
}





#///////////////////////////////////////////////////////////////////////////////////////

# +++ Function "TE_estimation"

TE_estimation <-  function (responses, restrictions, equal.prob.of.error.in.choice, fit_method) {
  
  lower_values <- c(0.000001,0.000001,0.000001,0.000001)
  upper_values <- c(1,1,1,1,.5,.5,.5,.5)
  
  # "restrictions" can constrain any of 4 true probabilities: a_00, a_01,a_10,a_11
  # "FREE" indicates parameters unconstrained
  suppressWarnings(new_upper_values <- ifelse(restrictions=="FREE", 
                                              upper_values, as.numeric(restrictions)+0.000001))
  suppressWarnings(new_lower_values <- ifelse(restrictions=="FREE", 
                                              lower_values, as.numeric(restrictions)))
  
  new_lower_values <- ifelse(new_lower_values==0,new_lower_values+0.000001, new_lower_values)
  new_upper_values <- ifelse(new_upper_values==new_lower_values, new_upper_values+0.000001,new_upper_values)
  
  initial_values <- (new_upper_values-new_lower_values ) / 2
  initial_values <- c(initial_values,.25,.25,.25,.25)
  upper_values <-c(new_upper_values,.5,.5,.5,.5)
  lower_values <-c(new_lower_values,0.000001,0.000001,0.000001,0.000001)
  
  if (equal.prob.of.error.in.choice=="TE-2") {
    initial_values<-initial_values[c(1:5,6)]
    lower_values<-lower_values[c(1:5,6)]
    upper_values<-upper_values[c(1:5,6)]
  } else if (equal.prob.of.error.in.choice=="TE-1") {
    initial_values<-initial_values[c(1:5)]
    lower_values<-lower_values[c(1:5)]
    upper_values<-upper_values[c(1:5)]
  }
  
  # Running optimization
  R<- optim_TE_N_times(initial_values=initial_values,
                       lower_values=lower_values, 
                       upper_values=upper_values, 
                       responses=responses, 
                       restrictions=restrictions, 
                       equal.prob.of.error.in.choice=equal.prob.of.error.in.choice,
                       fit_method=fit_method,
                       N=3) 
  
  parameters <- R[[1]]
  estimated_errors <- R[[2]]
  sol <- R[[4]]
  
  colnames(parameters) <- c("Estimates")
  rownames(parameters) <- c("a_00", "a_01","a_10", "a_11")
  colnames(estimated_errors) <- c("Estimates")
  
  if (nrow(estimated_errors)==1) {
    rownames(estimated_errors) <- c("error")
  } else  if (nrow(estimated_errors)==2) {
    rownames(estimated_errors) <- c("error_1", "error_2")
  }   else  {
    rownames(estimated_errors) <- c("error_R1", "error_R2", "error_S1","error_S2")
  }
  
  # Chi-squared statistic
  Chi2_test_statistic <- as.data.frame(sol$value) 
  
  
  if (fit_method=="CHI2") {
    colnames(Chi2_test_statistic) <- c("Chi2 statistic")
  } else {
    colnames(Chi2_test_statistic) <- c("G2 deviance")
    
  }
  # P value from Chi-squared table & Degrees of freedom
  # -- We have 15 degrees of freedom in the data
  # -- We use "DF1 -1" degrees of freedom to estimate true probabilities
  # -- We use "DF2" degrees of freedom to estimate error terms
  
  DF1 <- sum(ifelse(restrictions=="FREE", 1, 0))
  
  ##########################
  if (equal.prob.of.error.in.choice=="TE-2") {
    DF2=2} else if (equal.prob.of.error.in.choice=="TE-1") {
      DF2=1} else {DF2=4}
  
  df=(15-(DF1-1)-DF2) # total degrees of freedom
  
  Chi2_table_p_value <- as.data.frame(1-pchisq(sol$value,df=df ) ) 
  colnames(Chi2_table_p_value ) <-  c(paste("p value, DF=", df))
  
  # Contrast between observed and predicted proportions of responses
  
  # Observed proportions
  obs_prop <- responses/sum(responses) 
  
  # Predicted proportions
  
  if (nrow(estimated_errors)==1) {
    
    pred_prop <-true_error(prob=c(parameters[,1],estimated_errors[1,1],estimated_errors[1,1],
                                  estimated_errors[1,1],estimated_errors[1,1]),    
                           responses=responses, 
                           restrictions=restrictions,
                           fit_method = fit_method)[[2]]
  } else   if (nrow(estimated_errors)==2) {
    
    pred_prop <-true_error(prob=c(parameters[,1],estimated_errors[1,1],estimated_errors[2,1],
                                  estimated_errors[1,1],estimated_errors[2,1]),    
                           responses=responses, 
                           restrictions=restrictions,
                           fit_method = fit_method)[[2]]
  } else {
    pred_prop <-true_error(prob=c(parameters[,1],estimated_errors[,1]),
                           responses=responses, 
                           restrictions=restrictions,
                           fit_method = fit_method)[[2]]  
  }
  contrast_prop<-cbind(obs_prop,pred_prop)
  
  rownames(contrast_prop) <- c("P0000", "P0001","P0010","P0011",
                               "P0100","P0101","P0110","P0111",
                               "P1000","P1001","P1010","P1011",
                               "P1100","P1101","P1110","P1111")
  colnames(contrast_prop) <- c("Observed proportions", "Predicted TE proportions")
  
  # Contrast between observed and predicted frequencies of responses
  
  # Predicted frequencies
  size <- sum(responses) #number of people
  pred_freq <- pred_prop*size
  
  contrast_freq <- cbind(responses,pred_freq )
  
  
  contrast_freq
  
  rownames(contrast_freq) <- c("P0000", "P0001","P0010","P0011",
                               "P0100","P0101","P0110","P0111",
                               "P1000","P1001","P1010","P1011",
                               "P1100","P1101","P1110","P1111")
  colnames(contrast_freq) <- c("Observed frequencies", "Predicted TE frequencies")
  
  # Predicted frequencies in matrix form
  Replication1_RRp<-responses[1:4]
  Replication1_RSp<-responses[5:8]
  Replication1_SRp<-responses[9:12]
  Replication1_SSp<-responses[13:16]
  
  matrix_frequencies <-rbind(Replication1_RRp,Replication1_RSp,  Replication1_SRp,Replication1_SSp)
  rownames(matrix_frequencies) <- c("Replication1_RR'", "Replication1_RS'","Replication1_SR'","Replication1_SS'")
  colnames(matrix_frequencies) <- c("Replication2_RR'", "Replication2_RS'","Replication2_SR'","Replication2_SS'")
  
  Replication1_RRp<-pred_freq[1:4]
  Replication1_RSp<-pred_freq[5:8]
  Replication1_SRp<-pred_freq[9:12]
  Replication1_SSp<-pred_freq[13:16]
  
  matrix_TE_predictions <-rbind(Replication1_RRp,Replication1_RSp,  Replication1_SRp,Replication1_SSp)
  rownames(matrix_TE_predictions) <- c("Replication1_RR'", "Replication1_RS'","Replication1_SR'","Replication1_SS'")
  colnames(matrix_TE_predictions) <- c("Replication2_RR'", "Replication2_RS'","Replication2_SR'","Replication2_SS'")
  
  list(parameters=parameters, errors=estimated_errors, chisq=Chi2_test_statistic,
       p.value=Chi2_table_p_value,contrast.of.proportions=contrast_prop,
       contrast.of.frequencies=contrast_freq, 
       matrix.contrast.frequencies=list(matrix.of.frequencies=matrix_frequencies,
                                        matrix.of.TE.predictions=matrix_TE_predictions))
}

#///////////////////////////////////////////////////////////////////////////////////////

# +++ Function "TE_estimation_data"


TE_estimation_data <- function(data,restrictions, equal.prob.of.error.in.choice,
                               fit_method) {
  
  data_list<-apply(data,1,list)
  number=nrow(data)
  
  results_list <-list() 
  
  # Fitting TE to each participant
  for(i in 1:number) {
    responses=data_list[[i]][[1]] [c(2:17)]
    results_list[[i]]=TE_estimation(responses=responses, 
                                    restrictions=restrictions, 
                                    equal.prob.of.error.in.choice=equal.prob.of.error.in.choice,
                                    fit_method=fit_method)
    
    print(paste0("Fitting True and Error Model ", equal.prob.of.error.in.choice ," for case:  ", i))  
  }
  
  results_list
  
}



#///////////////////////////////////////////////////////////////////////////////////////

# +++ Function "aggregate_estimates"

aggregate_estimates <- function (results, model, fit_method) {
  n_elements <- length(results)
  fit_method 
  if (model=="TE") {
    
    parameters_data <- data.frame()
    errors_data <- data.frame()
    chisq_data <- data.frame()
    p.value_data <- data.frame()
    contrast_frequencies_data_Obs <- data.frame()
    contrast_frequencies_data_Pred <- data.frame()
    
    for (i in 1:n_elements) {
      parameters_results <- as.data.frame(t(results[[i]]$parameters))
      parameters_data <- rbind(parameters_data, parameters_results)
      
      errors_results <- as.data.frame(t(results[[i]]$errors))
      errors_data <- rbind(errors_data, errors_results)
      
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
    
    all_results <- cbind(contrast_frequencies_data_Obs, contrast_frequencies_data_Pred, 
                         parameters_data, errors_data, chisq_data, p.value_data)  
    
    names(all_results)[17:32] <-  c("Prediction P0000", 
                                    "Prediction P0001",
                                    "Prediction P0010",
                                    "Prediction P0011",
                                    "Prediction P0100",
                                    "Prediction P0101",
                                    "Prediction P0110",
                                    "Prediction P0111",
                                    "Prediction P1000",
                                    "Prediction P1001",
                                    "Prediction P1010",
                                    "Prediction P1011", 
                                    "Prediction P1100", 
                                    "Prediction P1101", 
                                    "Prediction P1110", 
                                    "Prediction P1111")
    
    if (fit_method=="CHI2"){
      names(all_results)[ (ncol(all_results)-1 ):ncol(all_results)] <- c("Chi-squared statistic", colnames(results[[1]]$p.value))
    } else {
      names(all_results)[ (ncol(all_results)-1 ):ncol(all_results)] <- c("G2 deviance", colnames(results[[1]]$p.value))
    }
    
  }
  
  
  
  if (model=="Independence") {
    chisq_data <- data.frame()
    contrast_frequencies_data_Obs <- data.frame()
    contrast_frequencies_data_Pred <- data.frame()
    Marginal.prob.Independence_data <-data.frame()
    
    for (i in 1:n_elements) {
      chisq_results <- as.data.frame(t(results[[i]]$chisq))
      chisq_data <- rbind(chisq_data, chisq_results)
      
      contrast_frequencies_results_Obs <- as.data.frame(t(results[[i]]$contrast.of.frequencies[,1]))
      contrast_frequencies_results_Pred <- as.data.frame(t(results[[i]]$contrast.of.frequencies[,2]))
      marg.prob <- as.data.frame(t(results[[i]]$Marginal.prob.Independence[,1]  ))
      
      contrast_frequencies_data_Obs <- rbind(contrast_frequencies_data_Obs, 
                                             contrast_frequencies_results_Obs)
      contrast_frequencies_data_Pred <- rbind(contrast_frequencies_data_Pred, 
                                              contrast_frequencies_results_Pred)
      Marginal.prob.Independence_data <- rbind(Marginal.prob.Independence_data, 
                                               marg.prob) 
    }
    
    all_results <- cbind(Marginal.prob.Independence_data, contrast_frequencies_data_Obs, contrast_frequencies_data_Pred, 
                         chisq_data)  
    
    names(all_results)[25:40] <-  c("Prediction P0000", 
                                    "Prediction P0001",
                                    "Prediction P0010",
                                    "Prediction P0011",
                                    "Prediction P0100",
                                    "Prediction P0101",
                                    "Prediction P0110",
                                    "Prediction P0111",
                                    "Prediction P1000",
                                    "Prediction P1001",
                                    "Prediction P1010",
                                    "Prediction P1011", 
                                    "Prediction P1100", 
                                    "Prediction P1101", 
                                    "Prediction P1110", 
                                    "Prediction P1111")
    
    
    
    if (fit_method=="CHI2"){
      names(all_results)[ ncol(all_results)] <- c("Chi-squared statistic")
    } else {
      names(all_results)[ ncol(all_results)] <- c("G2 deviance")
    }
    
    
    
  }
  rownames(all_results) <- NULL
  all_results
}


#///////////////////////////////////////////////////////////////////////////////////////
# ____________________________________________________________
#
# FUNCTIONS TO PERFORM MONTE CARLO SIMULATIONS FOR TE MODEL
# ____________________________________________________________
#
#///////////////////////////////////////////////////////////////////////////////////////

# +++ Function "extract_results"


extract_results <- function(results_list_samples_in) {
  chi_values_restricted <- results_list_samples_in[[4]]$value
  parameters_all_restricted <- results_list_samples_in[[1]][,1]
  estimated_errors_all_restricted <- results_list_samples_in[[2]][,1]
  list(chi_values_restricted, parameters_all_restricted ,  estimated_errors_all_restricted)
}

#///////////////////////////////////////////////////////////////////////////////////////

# +++ Function "MonteCarlo_TE"

MonteCarlo_TE <- function(number_samples, 
                          population_size, 
                          expected_frequencies_under_null, 
                          chisq.statistic,
                          restrictions,
                          equal.prob.of.error.in.choice,
                          fit_method) {
  
  # Random samples are generated according to the null of TE (i.e., predicted 
  # probabilities of the model). 
  # For the Re-Fitted Monte Carlo simulation: Each sample is then fit with TE, generating  
  # a set of parameters (true probabilities and errors) and the corresponding
  # Chi-squared statistics (or G2).
  # For the Conservative Monte Carlo simulation: Each sample is compared to the 
  # predicted frequencies using the original parameter estimates of the models.
  
  samples <- rmultinom(number_samples, population_size, expected_frequencies_under_null) 
  
  # In "samples", each column is a sample from a distribution (according to the null of 
  # TE). Each column has 16 values representing the 16 possible pattern of responses
  
  # Re-Fitted Monte Carlo
  # ---------------------
  
  responses_montecarlo <- samples
  
  # Determining starting values of the optimization algorithm 
  # Determining bounds of parameters according to the restrictions impossed 
  # by the user. Restrictions can be in any true probability: a_00, a_01,a_10,a_11              
  
  lower_values <- c(0.000001,0.000001,0.000001,0.000001)
  upper_values<-c(1,1,1,1,.5,.5,.5,.5)
  
  suppressWarnings(new_upper_values <- ifelse(restrictions=="FREE", 
                                              upper_values, as.numeric(restrictions)+0.000001))
  suppressWarnings(new_lower_values <- ifelse(restrictions=="FREE", 
                                              lower_values, as.numeric(restrictions)))
  
  new_lower_values <- ifelse(new_lower_values==0,new_lower_values+0.000001, new_lower_values)
  new_upper_values <- ifelse(new_upper_values==new_lower_values, new_upper_values+0.000001,new_upper_values)
  
  initial_values <- (new_upper_values-new_lower_values ) / 2
  initial_values <- c(initial_values,.25,.25,.25,.25)
  upper_values <-c(new_upper_values,.5,.5,.5,.5)
  lower_values <-c(new_lower_values,0.000001,0.000001,0.000001,0.000001)
  ##CORRECTION
  # rbind(upper_values + 1, lower_values + 1)
  
  if (equal.prob.of.error.in.choice=="TE-2") {
    initial_values<-initial_values[c(1:5,6)]
    lower_values<-lower_values[c(1:5,6)]
    upper_values<-upper_values[c(1:5,6)]
  } else if (equal.prob.of.error.in.choice=="TE-1") {
    initial_values<-initial_values[c(1:5)]
    lower_values<-lower_values[c(1:5)]
    upper_values<-upper_values[c(1:5)]
  }
  
  # Performing the optimization for each sample 
  results_list_samples <- lapply(seq_along(1:number_samples), function(i) {
    
    optim_TE_N_times(initial_values=initial_values,
                     lower_values=lower_values, 
                     upper_values=upper_values, 
                     responses=responses_montecarlo[,i], 
                     restrictions=restrictions, 
                     equal.prob.of.error.in.choice=equal.prob.of.error.in.choice,
                     fit_method = fit_method,
                     N=1) 
  })
  
  results_extracted <- lapply(seq_along(1:number_samples), function(i) {
    extract_results(results_list_samples_in=results_list_samples[[i]])
  })
  
  # Extracting results from the optimization
  chi_values <-  matrix(unlist(    lapply(results_extracted , function(l) l[[1]])     ),  
                        ncol = 1, byrow = TRUE)
  parameters_all <- matrix( unlist( lapply(results_extracted , function(l) l[[2]])   ),
                            ncol = 4, byrow = TRUE)
  estimated_errors_all <- matrix(unlist( lapply(results_extracted , function(l) l[[3]]) ),
                                 ncol=(length(initial_values)-4), byrow=TRUE)
  
  colnames(parameters_all) <- c("a_00", "a_01","a_10", "a_11")
  if (ncol(estimated_errors_all)==1) {
    colnames(estimated_errors_all)<- c("error")  
  } else if (ncol(estimated_errors_all)==2) {
    colnames(estimated_errors_all) <- c("error_1", "error_2")
  }  else   {
    colnames(estimated_errors_all)<- c("error_R1", "error_R2", "error_S1","error_S2")  
  }
  
  
  
  if (fit_method=="CHI2") {
    colnames(chi_values) <- c("Chi squared statistics of Re-fitted MC")
  } else {
    colnames(chi_values) <- c("G2 deviance of Re-fitted MC")
  }
  
  
  statistic <- chisq.statistic
  MonteCarlo_p_value <- (length(chi_values[chi_values>=as.numeric(statistic)   ]    ) )/length(chi_values)
  
  # Conservative Monte Carlo
  # ------------------------
  
  p=expected_frequencies_under_null
  chi_all_conservative <-matrix(nrow=number_samples,ncol=1)
  
  for(i in 1:number_samples) {
    
    x=samples[,i]
    
    chi_sample <- ((x - p)^2)/p
    chi<-sum(chi_sample)
    chi_all_conservative[i,1]<-chi
  }
  
  p.value.conservative <- (length(chi_all_conservative[chi_all_conservative >= as.numeric(statistic) ]))/
    length(chi_all_conservative)
  
  
  
  if (fit_method=="CHI2") {
    
    colnames(chi_all_conservative) <- c("Chi squared statistics of Conservative MC")
  } else {
    colnames(chi_all_conservative) <- c("G2 deviance of Conservative MC")
    
  }
  
  
  
  
  
  
  list(parameters_Re.fitted.MC=parameters_all, errors_Re.fitted.MC= estimated_errors_all, 
       chisq_Re.fitted.MC=chi_values, p.value_Re.fitted.MC=MonteCarlo_p_value, 
       chisq_Conservative.MC=chi_all_conservative, p.value_Conservative.MC=p.value.conservative,
       samples=samples, original.chisq=statistic)
  
}

#///////////////////////////////////////////////////////////////////////////////////////

# +++ Function "MonteCarlo_simulation_TE"

MonteCarlo_simulation_TE <- function(responses,
                                     restrictions,
                                     equal.prob.of.error.in.choice,
                                     number_samples,
                                     fit_method, i) {
  
  TE_results <- TE_estimation(responses=responses, restrictions=restrictions, 
                              equal.prob.of.error.in.choice=equal.prob.of.error.in.choice,
                              fit_method = fit_method)
  
  pred_freq <- TE_results$contrast.of.frequencies[,2]
  size <- sum(responses)  
  statistic <- TE_results$chisq
  
  simulation_MonteCarlo_TE <- MonteCarlo_TE(number_samples=number_samples,
                                            population_size=size,
                                            expected_frequencies_under_null=pred_freq, 
                                            chisq.statistic = statistic,
                                            restrictions=restrictions ,
                                            equal.prob.of.error.in.choice=equal.prob.of.error.in.choice,
                                            fit_method=fit_method)
  
  print(paste0("Performing Monte Carlo simulations of ", equal.prob.of.error.in.choice, " for case:  ", i)) 
  
  simulation_MonteCarlo_TE 
}

#///////////////////////////////////////////////////////////////////////////////////////

# +++ Function "TE_MonteCarlo_data"

TE_MonteCarlo_data <- function (data, 
                                restrictions,
                                equal.prob.of.error.in.choice,
                                number_samples,
                                fit_method=fit_method) {
  
  data_list<-apply(data,1,list)
  number=nrow(data)
  
  results_list <-list() 
  
  list_responses <- list()
  
  for(i in 1:number) {
    list_responses[[i]]=data_list[[i]][[1]] [c(2:17)]
  }
  
  results_list <- lapply(seq_along(1:number), function(i) {
    
    MonteCarlo_simulation_TE    (responses=list_responses[[i]],
                                 restrictions=restrictions,
                                 equal.prob.of.error.in.choice=equal.prob.of.error.in.choice,
                                 number_samples=number_samples, 
                                 fit_method=fit_method,
                                 i=i)
  })
  results_list
  
}

#///////////////////////////////////////////////////////////////////////////////////////

# +++ Function "aggregate_MonteCarlo_p.value"

aggregate_MonteCarlo_p.value <- function(results) {
  
  n_elements <- length(results)
  
  MonteCarlo.Re.fitted.p.value_data <- data.frame()
  MonteCarlo.Conservative.p.value_data <- data.frame()
  
  for (i in 1:n_elements) {
    #i=1
    MonteCarlo.Re.fitted.p.value_results <- as.data.frame(results[[i]]$p.value_Re.fitted.MC)
    MonteCarlo.Re.fitted.p.value_data <- rbind(MonteCarlo.Re.fitted.p.value_data, MonteCarlo.Re.fitted.p.value_results)
    
    MonteCarlo.Conservative.p.value_results <- as.data.frame(results[[i]]$p.value_Conservative.MC)
    MonteCarlo.Conservative.p.value_data <- rbind(MonteCarlo.Conservative.p.value_data, MonteCarlo.Conservative.p.value_results)
    
  }
  
  MonteCarlo.p.value_data_all<- cbind(MonteCarlo.Conservative.p.value_data, MonteCarlo.Re.fitted.p.value_data )
  
  colnames(MonteCarlo.p.value_data_all) <- c("Conservative MC - p value", "Re-fitted MC - p value")
  
  
  MonteCarlo.p.value_data_all
}

#///////////////////////////////////////////////////////////////////////////////////////
# ________________________________________________________________________________
#
# FUNCTIONS TO PERFORM BOOTSTRAPPING SIMULATIONS TO OBTAIN CONFIDENCE INTERVALS
# ________________________________________________________________________________
#
#///////////////////////////////////////////////////////////////////////////////////////
#
# +++ Function "CI_estimates_bootstrap"

CI_estimates_bootstrap <- function(obs_patterns_raw_data, 
                                   equal.prob.of.error.in.choice,
                                   restrictions,
                                   fit_method,
                                   d) {
  # obs_patterns_raw_data is a 1xNumber.of.participants vector in which each
  # cell represent an individual. 
  # Numbers on each cell will indicate the response pattern showed by the individual.
  # Numbers 1-16 represent response patterns: "P0000", "P0001","P0010","P0011",
  #  "P0100","P0101","P0110","P0111", "P1000","P1001","P1010","P1011",
  #  "P1100","P1101","P1110","P1111"
  # "d" has indexes that will be used to build each sample.
  # "obs_sample" has the frequencies of each pattern in a sample. 
  
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
  
  responses_boot <- obs_sample
  
  # Setting the bounds for the estimation
  lower_values <- c(0.000001,0.000001,0.000001,0.000001)
  upper_values<-c(1,1,1,1,.5,.5,.5,.5)
  
  
  # Restrictions can be set in any true probability: a_00, a_01,a_10,a_11              
  suppressWarnings(new_upper_values <- ifelse(restrictions=="FREE", 
                                              upper_values, as.numeric(restrictions)+0.000001))
  suppressWarnings(new_lower_values <- ifelse(restrictions=="FREE", 
                                              lower_values, as.numeric(restrictions)))
  
  new_lower_values <- ifelse(new_lower_values==0,new_lower_values+0.000001, new_lower_values)
  new_upper_values <- ifelse(new_upper_values==new_lower_values, new_upper_values+0.000001,new_upper_values)
  
  initial_values <- (new_upper_values-new_lower_values ) / 2
  initial_values <- c(initial_values,.25,.25,.25,.25)
  upper_values <-c(new_upper_values,.5,.5,.5,.5)
  lower_values <-c(new_lower_values,0.000001,0.000001,0.000001,0.000001)
  
  
  
  if (equal.prob.of.error.in.choice=="TE-2") {
    initial_values<-initial_values[c(1:5,6)]
    lower_values<-lower_values[c(1:5,6)]
    upper_values<-upper_values[c(1:5,6)]
  } else if (equal.prob.of.error.in.choice=="TE-1") {
    initial_values<-initial_values[c(1:5)]
    lower_values<-lower_values[c(1:5)]
    upper_values<-upper_values[c(1:5)]
  }
  
  
  
  
  # Running optimizations
  R<-   optim_TE_N_times(initial_values=initial_values,
                         lower_values=lower_values, 
                         upper_values=upper_values, 
                         responses=responses_boot, 
                         restrictions=restrictions, 
                         equal.prob.of.error.in.choice=equal.prob.of.error.in.choice,
                         fit_method=fit_method,
                         N=1) 
  
  parameters <- R[[1]][,1]
  estimated_errors <- R[[2]][,1]
  
  c(parameters, estimated_errors, obs_sample)
  
}

#///////////////////////////////////////////////////////////////////////////////////////

# +++ Function "CI_Bootstrap_TE"

CI_Bootstrap_TE <- function(responses, number_samples, 
                            restrictions,
                            equal.prob.of.error.in.choice,
                            fit_method,
                            i_number) {
  
  obs_freq <- responses #observed frequencies
  names(obs_freq) <-   c("P0000", "P0001","P0010","P0011",
                         "P0100","P0101","P0110","P0111",
                         "P1000","P1001","P1010","P1011",
                         "P1100","P1101","P1110","P1111")
  # Creating a vector in which each cell represents the responses of an individual
  obs_patterns_raw_data <- rep(1:16, obs_freq) 
  
  # Performing bootstrapping simulations
  result <- boot(obs_patterns_raw_data, 
                 CI_estimates_bootstrap, 
                 equal.prob.of.error.in.choice=equal.prob.of.error.in.choice,
                 restrictions=restrictions,
                 fit_method=fit_method,
                 R=number_samples)
  
  # Estimating confidence intervals (percentile method)
  
  number_estimates= (length(result$t0)-16)
  all_CI <- matrix(nrow= number_estimates,ncol=2)
  
  suppressWarnings(   for(i in 1:  number_estimates   ) {
    
    if (!is.null(  boot.ci(result, type ="perc",index=i)[[4]][4:5]  )) { 
      all_CI[i,] <-  boot.ci(result, type ="perc",index=i)[[4]][4:5]  }
  } )
  
  table_CI <- cbind(result$t0[1: number_estimates    ], all_CI)
  
  if (number_estimates==8) {
    parnames = c("a_00", "a_01","a_10", "a_11", 
                 "error_R1", "error_R2", "error_S1","error_S2")
  } else  if (number_estimates==6) {
    parnames = c("a_00", "a_01","a_10", "a_11", 
                 "error_1", "error_2")
  } else  {
    parnames = c("a_00", "a_01","a_10", "a_11", 
                 "error")
  }
  
  
  rownames(table_CI) = parnames
  columnnames = c("Estimates", "Lower bound", "Upper bound")
  colnames(table_CI) = columnnames
  
  boot_parameters <- result$t[,1:4]
  boot_errors <- result$t[,5:number_estimates]
  boot_samples <- result$t[,(number_estimates+1):(number_estimates+16)]
  
  colnames(boot_parameters) <- c("a_00", "a_01","a_10", "a_11")
  
  if (number_estimates==8) {
    colnames(boot_errors) <- c("error_R1", "error_R2", "error_S1","error_S2")
  } else  if (number_estimates==6) {
    colnames(boot_errors) <- c( "error_1", "error_2")
  } else  {
    boot_errors<-as.data.frame(boot_errors)
    colnames(boot_errors) <- c( "error")
  }
  
  colnames(boot_samples) <- c("P0000", "P0001","P0010","P0011",
                              "P0100","P0101","P0110","P0111",
                              "P1000","P1001","P1010","P1011",
                              "P1100","P1101","P1110","P1111")
  
  print(paste0("Performing Bootstrap simulations of ",equal.prob.of.error.in.choice ," for case:  ", i_number)) 
  
  list(boot.parameters=boot_parameters, boot.errors=boot_errors, boot.samples=boot_samples, 
       confidence.intervals=table_CI)
}


#///////////////////////////////////////////////////////////////////////////////////////

# +++ Function "TE_Bootstrapping.CI_data"

TE_Bootstrapping.CI_data <- function (data, 
                                      number_samples,
                                      restrictions,
                                      equal.prob.of.error.in.choice,
                                      fit_method) {
  data_list<-apply(data,1,list)
  number=nrow(data)
  
  responses_list <- list()
  
  for(i in 1:number) {
    responses_list[[i]]=data_list[[i]][[1]][c(2:17)]
  }
  
  results_list <- lapply(seq_along(1:number), function(i) {
    CI_Bootstrap_TE(responses=responses_list[[i]], 
                    number_samples=number_samples, 
                    restrictions=restrictions,
                    equal.prob.of.error.in.choice=equal.prob.of.error.in.choice,
                    fit_method=fit_method,
                    i_number = i)
  })
  results_list
}

#///////////////////////////////////////////////////////////////////////////////////////

# +++ Function "aggregate_boot.confidence.intervals"

aggregate_boot.confidence.intervals <- function (results) {
  
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
  
  if (  ncol(confidence.intervals_data)==24) {
    names(confidence.intervals_data)[9:16]  <-c("LowerBound a_00", 
                                                "LowerBound a_01",
                                                "LowerBound a_10", 
                                                "LowerBound a_11",
                                                "LowerBound error_R1", 
                                                "LowerBound error_R2", 
                                                "LowerBound error_S1",
                                                "LowerBound error_S2" )
    names(confidence.intervals_data)[17:24] <- c("UpperBound a_00", 
                                                 "UpperBound a_01",
                                                 "UpperBound a_10", 
                                                 "UpperBound a_11",
                                                 "UpperBound error_R1", 
                                                 "UpperBound error_R2", 
                                                 "UpperBound error_S1",
                                                 "UpperBound error_S2" )
  } else   if (  ncol(confidence.intervals_data)==18) {
    names(confidence.intervals_data)[7:12]  <-c("LowerBound a_00", 
                                                "LowerBound a_01",
                                                "LowerBound a_10", 
                                                "LowerBound a_11",
                                                "LowerBound error_1", 
                                                "LowerBound error_2")
    
    names(confidence.intervals_data)[13:18] <- c("UpperBound a_00", 
                                                 "UpperBound a_01",
                                                 "UpperBound a_10", 
                                                 "UpperBound a_11",
                                                 "UpperBound error_1", 
                                                 "UpperBound error_2")
  }  else  {
    names(confidence.intervals_data)[6:10]  <-c("LowerBound a_00", 
                                                "LowerBound a_01",
                                                "LowerBound a_10", 
                                                "LowerBound a_11",
                                                "LowerBound error")
    
    names(confidence.intervals_data)[11:15] <- c("UpperBound a_00", 
                                                 "UpperBound a_01",
                                                 "UpperBound a_10", 
                                                 "UpperBound a_11",
                                                 "UpperBound error")
  }
  confidence.intervals_data
}


##########################################################################################
#
# FUNCTIONS APPLIED IN THE ESTIMATION OF INDEPENDENCE MODEL
#
# - To estimate Independence model parameters
#    +++ Independence: computes prediction of the model using responses of an individual
#    +++ Independence_estimation_data: computes predictions for all the individuals presented 
#                            in the researcher's data file
#    +++ aggregate_estimates: accumulates results obtained by "Independence_estimation_data"
#                             in a data frame.
#
###########################################################################################


# +++ Function "Independence"

Independence <- function( responses, fit_method ){
  
  # Computing marginal probabilities from response patterns
  Marginal.prob.Replication1_R<-sum(responses[1:8])/sum(responses)
  Marginal.prob.Replication1_S<-sum(responses[9:16])/sum(responses)
  Marginal.prob.Replication1_Rp<-sum(responses[1:4]+responses[9:12])/sum(responses)
  Marginal.prob.Replication1_Sp<-sum(responses[5:8]+responses[13:16])/sum(responses)
  Marginal.prob.Replication2_R<-sum(responses[c(1, 5, 9, 13, 2, 6, 10, 14)])/sum(responses)
  Marginal.prob.Replication2_S<-sum(responses[c(3, 7, 11, 15, 4, 8, 12, 16)])/sum(responses)
  Marginal.prob.Replication2_Rp<-sum(responses[c(1, 5, 9, 13, 3, 7, 11, 15)])/sum(responses)
  Marginal.prob.Replication2_Sp<-sum(responses[c(2, 6, 10, 14, 4, 8, 12, 16)])/sum(responses)
  
  Marginal.prob.Independence<-rbind(Marginal.prob.Replication1_R,
                                    Marginal.prob.Replication1_S,
                                    Marginal.prob.Replication1_Rp,
                                    Marginal.prob.Replication1_Sp,
                                    Marginal.prob.Replication2_R,
                                    Marginal.prob.Replication2_S,
                                    Marginal.prob.Replication2_Rp,
                                    Marginal.prob.Replication2_Sp)
  
  rownames(Marginal.prob.Independence)<-c("Marginal.prob.Replication1_R",
                                          "Marginal.prob.Replication1_S",
                                          "Marginal.prob.Replication1_R'",
                                          "Marginal.prob.Replication1_S'",
                                          "Marginal.prob.Replication2_R",
                                          "Marginal.prob.Replication2_S",
                                          "Marginal.prob.Replication2_R'",
                                          "Marginal.prob.Replication2_S'")
  
  V1<- c(0,	0,	0,	0,	0,	0,	0,	0,	1,	1,	1,	1,	1,	1,	1,	1)
  V2<- c(0,	0,	0,	0,	1,	1,	1,	1,	0,	0,	0,	0,	1,	1,	1,	1)
  V3<- c(0,	0,	1,	1,	0,	0,	1,	1,	0,	0,	1,	1,	0,	0,	1,	1)
  V4<- c(0,	1,	0,	1,	0,	1,	0,	1,	0,	1,	0,	1,	0,	1,	0,	1)
  
  M1<-ifelse(V1==1,Marginal.prob.Replication1_S,1-Marginal.prob.Replication1_S)
  M2<-ifelse(V2==1,Marginal.prob.Replication1_Sp, 1-Marginal.prob.Replication1_Sp)
  M3<-ifelse(V3==1,Marginal.prob.Replication2_S,1-Marginal.prob.Replication2_S)
  M4<-ifelse(V4==1,Marginal.prob.Replication2_Sp,1-Marginal.prob.Replication2_Sp)
  
  # Predictions of the model
  Independence.predictions <- (M1*M2*M3*M4)* sum(responses)
  
  
  # Computing Chi-squared values
  
  if (fit_method=="CHI2"){
    Chi_vector <- ((responses - Independence.predictions)^2)/
      (Independence.predictions)
    Chi_total <- sum(Chi_vector)
  } else{
    
    # Computing the G2 index or deviance statistic (default method)
    
    G_data <-  as.data.table(cbind(responses , G= log(responses / (Independence.predictions) ) ))
    G_data[responses==0, G:=0 ]
    G_data[, G_vector:=responses*G]  
    G_index <- 2*sum(G_data$G_vector)
    Chi_total <- G_index
  }
  
  
  
  
  
  
  Replication1_RRp<-responses[1:4]
  Replication1_RSp<-responses[5:8]
  Replication1_SRp<-responses[9:12]
  Replication1_SSp<-responses[13:16]
  
  matrix_frequencies <-rbind(Replication1_RRp,Replication1_RSp,  Replication1_SRp,Replication1_SSp)
  rownames(matrix_frequencies) <- c("Replication1_RR'", "Replication1_RS'","Replication1_SR'","Replication1_SS'")
  colnames(matrix_frequencies) <- c("Replication2_RR'", "Replication2_RS'","Replication2_SR'","Replication2_SS'")
  
  Replication1_RRp<-Independence.predictions[1:4]
  Replication1_RSp<-Independence.predictions[5:8]
  Replication1_SRp<-Independence.predictions[9:12]
  Replication1_SSp<-Independence.predictions[13:16]
  
  matrix_Independence_predictions <-rbind(Replication1_RRp,Replication1_RSp,  Replication1_SRp,Replication1_SSp)
  rownames(matrix_Independence_predictions) <- c("Replication1_RR'", "Replication1_RS'","Replication1_SR'","Replication1_SS'")
  colnames(matrix_Independence_predictions) <- c("Replication2_RR'", "Replication2_RS'","Replication2_SR'","Replication2_SS'")
  
  # Constrans of observed and predicted frequencies
  contrast_freq <- cbind(responses,Independence.predictions)
  rownames(contrast_freq) <- c("P0000", "P0001","P0010","P0011",
                               "P0100","P0101","P0110","P0111",
                               "P1000","P1001","P1010","P1011",
                               "P1100","P1101","P1110","P1111")
  colnames(contrast_freq) <- c("Observed frequencies", "Predicted TE frequencies")
  
  list( chisq=Chi_total,Marginal.prob.Independence=Marginal.prob.Independence,
        contrast.of.frequencies=contrast_freq, 
        matrix.contrast.frequencies=list(matrix.of.frequencies=matrix_frequencies,
                                         matrix.of.Independence.predictions=matrix_Independence_predictions))
}

#///////////////////////////////////////////////////////////////////////////////////////

# +++ Function "Independence_estimation_data"

Independence_estimation_data <- function(data, fit_method){
  
  data_list <- apply(data,1,list)
  number=nrow(data)
  
  results_list <-list() 
  for(i in 1:number) {
    responses=data_list[[i]][[1]] [c(2:17)]
    results_list[[i]]=Independence(responses=responses, fit_method=fit_method)
  }
  
  results_list
  
}

#///////////////////////////////////////////////////////////////////////////////////////

# +++ Function "aggregate_estimates" (already defined in functions for TE estimations)
#

#///////////////////////////////////////////////////////////////////////////////////////
#
#


##########################################################################################
#
# FUNCTIONS APPLIED TO READ THE DATA, ESTIMATE SELECTED MODELS, AND COMPARE COMPETING MODELS
#
# - To perform Chi2 (or G2) difference test between restricted and unrestricted models
#    +++ test_chi2_diff_case: Gathers models' predictions for a participant. 
#                             Obtains predicted frequencies for the unrestricted model(s)
#                             and Chi2 value(s) 
#
#    +++ Chi2_diff: Computes Chi2 difference test for a given participant using 
#                   the predictions of the unrestricted model. 
#
#    +++ test_chi2_diff_data: Computes Chi2 difference test for all the participants
#
#
# - To perform Monte Carlo Chi2 (or G2) difference test between restricted and unrestricted models
#
#    +++ MonteCarlo_test_chi2_diff_case: Gathers models' predictions for a participant. 
#                             Obtains predicted frequencies for the unrestricted model(s)
#                             and Chi2 value(s) 
#
#    +++ MonteCarlo_Chi2_diff: Computes Monte Carlo Chi2 difference test for a given
#                              participant using the predictions of the unrestricted model. 
#
#    +++ MonteCarlo_Chi2_diff: Computes Monte Carlo Chi2 difference test for all the participants
#                              Restuls include:
#                              (1) standard Chi2 difference test
#                              (2) Monte Carlo Re-fitted Chi2 difference test
#                              (3) Monte Carlo Conservative Chi2 difference test
#
#
# - To read the data and estimate the models selected by the user
#    +++ TE_read_inputs_and_MC: reads data file and compute Monte Carlo and Bootstrap simulations 
#                               if required by the user (for the selected models). 
#                               Results are exported in csv files.
# 
#    +++ TE_READ_DATA: Runs funtion 'TE_read_inputs_and_MC' and adds to the results tests to 
#                      compare competing models if required by the user. 
###########################################################################################



#///////////////////////////////////////////////////////////////////////////////////////
# _______________________________________________________________________________________
#
#    FUNCTIONS TO PERFORM CHI2 (or G2) DIFFERENCE TEST BETWEEN RESTRICTED AND UNRESTRICTED MODELS
# _______________________________________________________________________________________
#
#///////////////////////////////////////////////////////////////////////////////////////

# +++ Function "test_chi2_diff_case"
#
test_chi2_diff_case <- function(simulation_results, case) {
  # Gathering results from simulations based on the restricted model(s). They are stored in 'simulation_results'  
  TE_1_estimation <- simulation_results$results_TE_1_errors[[case]]
  TE_2_estimation <- simulation_results$results_TE_2_errors[[case]]
  TE_4_estimation <- simulation_results$results_TE_4_errors[[case]]
  restrictions <- simulation_results$user_restrictions
  fit_method <- simulation_results$fit_method
  # Obtaining predicted frequencies for the unrestricted model(s) & Chi2 (or G2) value(s) 
  
  observed_responses <- NULL
  chisq_TE_1_unrestricted <- NULL
  chisq_TE_2_unrestricted <- NULL
  chisq_TE_4_unrestricted <- NULL
  TE_1_unrestricted_results <-NULL
  TE_2_unrestricted_results <-NULL
  TE_4_unrestricted_results <-NULL
  
  
  chisq_TE_1_restricted <- TE_1_estimation$chisq
  chisq_TE_2_restricted <- TE_2_estimation$chisq
  chisq_TE_4_restricted <- TE_4_estimation$chisq 
  
  
  
  if (is.null(TE_1_estimation)==0) {
    observed_responses <-TE_1_estimation$contrast.of.frequencies[,1]
    TE_1_unrestricted_results <- TE_estimation(responses=observed_responses, 
                                               restrictions=c("FREE", "FREE", "FREE", "FREE"),
                                               equal.prob.of.error.in.choice="TE-1",
                                               fit_method = fit_method)
    predicted_freq_TE_1_unrestricted <- TE_1_unrestricted_results$contrast.of.frequencies[,2]
    chisq_TE_1_unrestricted <-   TE_1_unrestricted_results$chisq
    
  } 
  if (is.null(TE_2_estimation)==0) {
    observed_responses <-TE_2_estimation$contrast.of.frequencies[,1]
    TE_2_unrestricted_results <- TE_estimation(responses=observed_responses, 
                                               restrictions=c("FREE", "FREE", "FREE", "FREE"),
                                               equal.prob.of.error.in.choice="TE-2",
                                               fit_method = fit_method)
    predicted_freq_TE_2_unrestricted <- TE_2_unrestricted_results$contrast.of.frequencies[,2]
    chisq_TE_2_unrestricted <-   TE_2_unrestricted_results$chisq
    
  } 
  if (is.null(TE_4_estimation)==0){
    observed_responses <-TE_4_estimation$contrast.of.frequencies[,1]
    TE_4_unrestricted_results <- TE_estimation(responses=observed_responses, 
                                               restrictions=c("FREE", "FREE", "FREE", "FREE"),
                                               equal.prob.of.error.in.choice="TE-4",
                                               fit_method = fit_method)
    predicted_freq_TE_4_unrestricted <- TE_4_unrestricted_results$contrast.of.frequencies[,2]
    chisq_TE_4_unrestricted <-   TE_4_unrestricted_results$chisq
    
  }
  
  
  # Running Chi2 difference test
  
  TE_1_Chi2_diff_test <- NULL
  TE_2_Chi2_diff_test <- NULL
  TE_4_Chi2_diff_test <- NULL
  
  if (is.null(TE_1_estimation)==0 ) {
    TE_1_Chi2_diff_test <- Chi2_diff(           predicted_freq_unrestricted=predicted_freq_TE_1_unrestricted,
                                                chisq_restricted=chisq_TE_1_restricted,
                                                chisq_unrestricted=chisq_TE_1_unrestricted,
                                                restrictions_tested=restrictions ) 
    
  }
  
  if (is.null(TE_2_estimation)==0 ) {
    TE_2_Chi2_diff_test <- Chi2_diff(           predicted_freq_unrestricted=predicted_freq_TE_2_unrestricted,
                                                chisq_restricted=chisq_TE_2_restricted,
                                                chisq_unrestricted=chisq_TE_2_unrestricted,
                                                restrictions_tested=restrictions) 
    
  }
  
  if (is.null(TE_4_estimation)==0 ) {
    TE_4_Chi2_diff_test <- Chi2_diff(           predicted_freq_unrestricted=predicted_freq_TE_4_unrestricted,
                                                chisq_restricted=chisq_TE_4_restricted,
                                                chisq_unrestricted=chisq_TE_4_unrestricted,
                                                restrictions_tested=restrictions ) 
    
  }
  
  print(paste0("Performing Chi-square (or G2) difference test of selected models for case:  ", case)) 
  
  
  list(  TE_1_Chi2_diff_test=TE_1_Chi2_diff_test,
         TE_2_Chi2_diff_test=TE_2_Chi2_diff_test,
         TE_4_Chi2_diff_test=TE_4_Chi2_diff_test,
         TE_1_unrestricted_results=TE_1_unrestricted_results,
         TE_2_unrestricted_results=TE_2_unrestricted_results,
         TE_4_unrestricted_results=TE_4_unrestricted_results,
         fit_method=fit_method)
  
}  

#///////////////////////////////////////////////////////////////////////////////////////

# +++ Function "Chi2_diff" 
#

Chi2_diff <- function(           predicted_freq_unrestricted,
                                 chisq_restricted,
                                 chisq_unrestricted,
                                 restrictions_tested)  {
  
  
  # ______________________________________________
  #
  #  Standard Chi2 (or G2) difference test
  # ______________________________________________
  
  Chi2.diff.data <-  chisq_restricted -  chisq_unrestricted
  
  val=names(chisq_restricted)
  # P value from Chi-squared table
  
  df <-   4-sum(restrictions_tested=="FREE")  #how many restrictions are imposed
  Chi2.diff.p.value <-  as.data.frame(1-pchisq(as.numeric(Chi2.diff.data)  ,df=df)) 
  colnames(Chi2.diff.p.value) <-  c(paste("p value, DF=", df))
  
  
  colnames(chisq_unrestricted) <- c(paste(val," - Unrestricted Model"))
  colnames(chisq_restricted) <- c(paste(val," - Restricted Model"))
  
  colnames(Chi2.diff.data) <-   c(paste("Difference in ", val, " of Restricted and Unrestricted TE"))
  
  
  
  
  
  Chi2_diff_results <- cbind(chisq_restricted,chisq_unrestricted,Chi2.diff.data, Chi2.diff.p.value )
  
  
  # If the chi2-diff-value (from Chi-squared table or from Monte Carlo simulations) is significant (p<.05), the 
  # Unrestricted model fits the data better than the Restricted model in which a number of parameters are 
  # fixed to certain value.
  # Otherwise, both models fit equally well and we can accept the Restricted model and fix some parameters.
  
  
  list( Chi2_diff_standard=Chi2_diff_results )
}





#///////////////////////////////////////////////////////////////////////////////////////

# +++ Function "test_chi2_diff_data"

test_chi2_diff_data <- function(simulation_results) {
  
  number <- simulation_results$number_cases
  results_list2 <- lapply(seq_along(1:number), function(i) {
    
    test_chi2_diff_case(simulation_results=simulation_results,
                        case=i)
    
  })
  
  results_list2
}



#/////////////////////////////////////////////////////////////////////////////////////////////////////
# ____________________________________________________________________________________________________
#
#    FUNCTIONS TO PERFORM MONTE CARLO CHI2 (OR G2) DIFFERENCE TEST BETWEEN RESTRICTED AND UNRESTRICTED MODELS
# ____________________________________________________________________________________________________
#
#/////////////////////////////////////////////////////////////////////////////////////////////////////

# +++ Function "MonteCarlo_test_chi2_diff_case"
#
MonteCarlo_test_chi2_diff_case <- function(simulation_results, case) {
  
  fit_method=simulation_results$fit_method
  # Gathering results from simulations based on the restricted model(s). They are stored in 'simulation_results'  
  TE_1_estimation <- simulation_results$results_TE_1_errors[[case]]
  TE_2_estimation <- simulation_results$results_TE_2_errors[[case]]
  TE_4_estimation <- simulation_results$results_TE_4_errors[[case]]
  TE_1_MC_estimation <- simulation_results$results_TE_1_MC[[case]]
  TE_2_MC_estimation <-  simulation_results$results_TE_2_MC[[case]]  
  TE_4_MC_estimation <- simulation_results$results_TE_4_MC[[case]] 
  restrictions <- simulation_results$user_restrictions
  
  # Gathering samples and simulated Monte Carlo Chi2 values for the restricted model(s)
  
  chisq_Re.fitted.MC_TE_1_restricted <- TE_1_MC_estimation$chisq_Re.fitted.MC
  samples_TE_1_under_restricted <- TE_1_MC_estimation$samples
  
  chisq_Re.fitted.MC_TE_2_restricted <- TE_2_MC_estimation$chisq_Re.fitted.MC
  samples_TE_2_under_restricted <- TE_2_MC_estimation$samples
  
  chisq_Re.fitted.MC_TE_4_restricted <- TE_4_MC_estimation$chisq_Re.fitted.MC
  samples_TE_4_under_restricted <- TE_4_MC_estimation$samples
  
  
  # Obtaining Chi2 value(s)  for the unrestricted model(s) 
  
  observed_responses <- NULL
  TE_1_unrestricted_results <- NULL
  TE_2_unrestricted_results <- NULL
  TE_4_unrestricted_results <- NULL
  chisq_TE_1_unrestricted <- NULL
  chisq_TE_2_unrestricted <- NULL
  chisq_TE_4_unrestricted <- NULL
  chisq_TE_1_restricted <- TE_1_estimation$chisq
  chisq_TE_2_restricted <- TE_2_estimation$chisq
  chisq_TE_4_restricted <- TE_4_estimation$chisq 
  
  if (is.null(TE_1_estimation)==0) {
    observed_responses <-TE_1_estimation$contrast.of.frequencies[,1]
    TE_1_unrestricted_results <- TE_estimation(responses=observed_responses, 
                                               restrictions=c("FREE", "FREE", "FREE", "FREE"),
                                               equal.prob.of.error.in.choice="TE-1",
                                               fit_method = fit_method)
    chisq_TE_1_unrestricted <-   TE_1_unrestricted_results$chisq
    
  } 
  
  if (is.null(TE_2_estimation)==0) {
    observed_responses <-TE_2_estimation$contrast.of.frequencies[,1]
    TE_2_unrestricted_results <- TE_estimation(responses=observed_responses, 
                                               restrictions=c("FREE", "FREE", "FREE", "FREE"),
                                               equal.prob.of.error.in.choice="TE-2",
                                               fit_method = fit_method)
    chisq_TE_2_unrestricted <-   TE_2_unrestricted_results$chisq
    
  } 
  
  if (is.null(TE_4_estimation)==0) {
    observed_responses <-TE_4_estimation$contrast.of.frequencies[,1]
    TE_4_unrestricted_results <- TE_estimation(responses=observed_responses, 
                                               restrictions=c("FREE", "FREE", "FREE", "FREE"),
                                               equal.prob.of.error.in.choice="TE-4",
                                               fit_method = fit_method)
    chisq_TE_4_unrestricted <-   TE_4_unrestricted_results$chisq
    
  }
  
  
  # Running Monte Carlo Chi2 difference test
  
  TE_1_Chi2_diff_test <- NULL
  TE_2_Chi2_diff_test <- NULL
  TE_4_Chi2_diff_test <- NULL
  
  if (is.null(TE_1_MC_estimation)==0 ) {
    TE_1_Chi2_diff_test <- MonteCarlo_Chi2_diff(samples_under_restricted=samples_TE_1_under_restricted,
                                                chisq_Re.fitted.MC=chisq_Re.fitted.MC_TE_1_restricted,
                                                chisq_restricted=chisq_TE_1_restricted,
                                                chisq_unrestricted=chisq_TE_1_unrestricted,
                                                restrictions_tested=restrictions ,
                                                equal.prob.of.error.in.choice="TE-1",
                                                fit_method = fit_method) 
    
  }
  
  if (is.null(TE_2_MC_estimation)==0 ) {
    TE_2_Chi2_diff_test <- MonteCarlo_Chi2_diff(samples_under_restricted=samples_TE_2_under_restricted,
                                                chisq_Re.fitted.MC=chisq_Re.fitted.MC_TE_2_restricted,
                                                chisq_restricted=chisq_TE_2_restricted,
                                                chisq_unrestricted=chisq_TE_2_unrestricted,
                                                restrictions_tested=restrictions ,
                                                equal.prob.of.error.in.choice="TE-2",
                                                fit_method = fit_method) 
    
  }
  
  if (is.null(TE_4_MC_estimation)==0 ) {
    TE_4_Chi2_diff_test <- MonteCarlo_Chi2_diff(samples_under_restricted=samples_TE_4_under_restricted,
                                                chisq_Re.fitted.MC=chisq_Re.fitted.MC_TE_4_restricted,
                                                chisq_restricted=chisq_TE_4_restricted,
                                                chisq_unrestricted=chisq_TE_4_unrestricted,
                                                restrictions_tested=restrictions ,
                                                equal.prob.of.error.in.choice="TE-4",
                                                fit_method = fit_method) 
    
  }
  
  print(paste0("Performing Monte Carlo Chi-square (or G2) difference test of selected models for case:  ", case)) 
  
  
  list(  TE_1_Chi2_diff_test=TE_1_Chi2_diff_test,
         TE_2_Chi2_diff_test=TE_2_Chi2_diff_test,
         TE_4_Chi2_diff_test=TE_4_Chi2_diff_test,
         TE_1_unrestricted_results=TE_1_unrestricted_results,
         TE_2_unrestricted_results=TE_2_unrestricted_results,
         TE_4_unrestricted_results=TE_4_unrestricted_results)
  
}  

#///////////////////////////////////////////////////////////////////////////////////////

# +++ Function "MonteCarlo_Chi2_diff" 
#
MonteCarlo_Chi2_diff <- function(samples_under_restricted,
                                 chisq_Re.fitted.MC,
                                 chisq_restricted,
                                 chisq_unrestricted, 
                                 restrictions_tested ,
                                 equal.prob.of.error.in.choice,
                                 fit_method) {
  
  
  # ______________________________________________
  #
  #  Standard Chi2 (or G2) difference test
  # ______________________________________________
  
  Chi2.diff.data <-  chisq_restricted -  chisq_unrestricted
  val=fit_method
  
  # P value from Chi-squared table
  
  df <-   4-sum(restrictions_tested=="FREE")  #how many restrictions are imposed
  Chi2.diff.p.value <-  as.data.frame(1-pchisq(as.numeric(Chi2.diff.data)  ,df=df)) 
  colnames(Chi2.diff.p.value) <-  c(paste("p value, DF=", df))
  
  
  
  
  
  colnames(chisq_unrestricted) <- c(paste(val," - Unrestricted Model"))
  colnames(chisq_restricted) <- c(paste(val," - Restricted Model"))
  colnames(Chi2.diff.data) <-   c(paste("Difference in ", val, " of Restricted and Unrestricted TE"))
  
  
  
  Chi2_diff_results <- cbind(chisq_restricted,chisq_unrestricted,Chi2.diff.data, Chi2.diff.p.value )
  
  
  # If the chi2-diff-value (from Chi-squared table or from Monte Carlo simulations) is significant (p<.05), the 
  # Unrestricted model fits the data better than the Restricted model in which a number of parameters are 
  # fixed to certain value.
  # Otherwise, both models fit equally well and we can accept the Restricted model and fix some parameters.
  
  
  # ______________________________________________
  #
  # Re-Fitted Monte Carlo - Chi2 (or G2) difference test
  # ______________________________________________
  
  # The Re-Fitted Monte Carlo for the restricted model was already computed and set as an input of this function. 
  # But we need the Monte Carlo results for the unrestricted model, using samples created under the restricted model.
  
  # Computing Re-Fitted Monte Carlo for Unrestricted Model
  
  
  responses_montecarlo <- samples_under_restricted
  
  # Determining starting values of the optimization algorithm 
  
  lower_values <- c(0.000001,0.000001,0.000001,0.000001)
  upper_values<-c(1,1,1,1,.5,.5,.5,.5)
  
  initial_values <- (upper_values-lower_values ) / 2
  initial_values <- c(initial_values,.25,.25,.25,.25)
  upper_values <-c(upper_values,.5,.5,.5,.5)
  lower_values <-c(lower_values,0.000001,0.000001,0.000001,0.000001)
  
  
  if (equal.prob.of.error.in.choice=="TE-2") {
    initial_values<-initial_values[c(1:5,6)]
    lower_values<-lower_values[c(1:5,6)]
    upper_values<-upper_values[c(1:5,6)]
  } else if (equal.prob.of.error.in.choice=="TE-1") {
    initial_values<-initial_values[c(1:5)]
    lower_values<-lower_values[c(1:5)]
    upper_values<-upper_values[c(1:5)]
  }
  
  # Performing the optimization for each sample 
  number_samples <- ncol(responses_montecarlo)
  results_list_samples <- lapply(seq_along(1:number_samples), function(i) {
    
    optim_TE_N_times(initial_values=initial_values,
                     lower_values=lower_values, 
                     upper_values=upper_values, 
                     responses=responses_montecarlo[,i], 
                     restrictions=c("FREE","FREE" , "FREE","FREE"),
                     equal.prob.of.error.in.choice=equal.prob.of.error.in.choice,
                     fit_method = fit_method,
                     N=1) 
  })
  
  results_extracted <- lapply(seq_along(1:number_samples), function(i) {
    extract_results(results_list_samples_in=results_list_samples[[i]])
  })
  
  # Extracting results from the optimization
  chi_values_unrestricted <-  matrix(unlist(    lapply(results_extracted , function(l) l[[1]])     ),  
                                     ncol = 1, byrow = TRUE)
  
  
  
  
  
  colnames(chi_values_unrestricted) <- c(paste(val,"of Re-fitted MC - Unrestricted Model"))
  
  chi_values_restricted <- chisq_Re.fitted.MC
  colnames(chi_values_restricted) <- c(paste(val,"of Re-fitted MC - Restricted Modell"))
  
  
  Chi2.diff.MC <- chi_values_restricted - chi_values_unrestricted
  colnames(Chi2.diff.MC) <- c(paste("Difference in ", val, "of Restricted and Unrestricted TE - Re-fitted MC"))
  
  # P value from Monte Carlo simulations
  Chi2.diff.MonteCarlo_p_value <- (length(Chi2.diff.MC[Chi2.diff.MC>=as.numeric(Chi2.diff.data)   ]    ))/length(Chi2.diff.MC)
  
  
  #Storing results in 'Re.fitted' list
  
  MC_Re.fitted_chi2_diff_results <- cbind(chi_values_restricted , chi_values_unrestricted, Chi2.diff.MC)
  Re.fitted <- list(  MC_Re.fitted_chi2_diff_results=MC_Re.fitted_chi2_diff_results, 
                      MC.p_value=Chi2.diff.MonteCarlo_p_value   )
  
  
  
  
  # Output of the function: Standard Chi2 difference test, Re-fitted Monte Carlo test, Conservative Monte Carlo test
  list( Chi2_diff_standard=Chi2_diff_results,Chi2_diff_Re.fitted.MC=Re.fitted )
}





#///////////////////////////////////////////////////////////////////////////////////////

# +++ Function "MonteCarlo_test_chi2_diff_data"

MonteCarlo_test_chi2_diff_data <- function(simulation_results) {
  
  number <- simulation_results$number_cases
  results_list2 <- lapply(seq_along(1:number), function(i) {
    MonteCarlo_test_chi2_diff_case(simulation_results=simulation_results,
                                   case=i)
  })
  
  results_list2
}


#///////////////////////////////////////////////////////////////////////////////////////
# _______________________________________________________________________________________
#
#    FUNCTIONS TO READ THE DATA AND ESTIMATE THE MODELS SELECTED BY THE USER
# _______________________________________________________________________________________
#
#///////////////////////////////////////////////////////////////////////////////////////

# +++ Function "TE_read_inputs_and_MC"

TE_read_inputs_and_MC <- function (data_file) {
  # Specify sheet 
  data_inputs <-read_excel(data_file, sheet = "Inputs")
  response_patterns <-read_excel(data_file, sheet = "Participant responses")
  number_cases <- nrow(response_patterns)
  
  
  
  Model_TE_4<-unname(   unlist(data_inputs[ which(data_inputs$Inputs=='TE-4'),2]))
  Model_TE_2<-unname(   unlist(data_inputs[ which(data_inputs$Inputs=='TE-2'),2]))
  Model_TE_1<-unname(   unlist(data_inputs[ which(data_inputs$Inputs=='TE-1'),2]))
  Model_Independence<-unname(   unlist(data_inputs[ which(data_inputs$Inputs=='Independence Model'),2]))
  
  a_00_restriction <-unname(   unlist(data_inputs[ which(data_inputs$Inputs=='a_00'),2]))
  a_01_restriction <-unname(   unlist(data_inputs[ which(data_inputs$Inputs=='a_01'),2]))
  a_10_restriction <-unname(   unlist(data_inputs[ which(data_inputs$Inputs=='a_10'),2]))
  a_11_restriction <-unname(   unlist(data_inputs[ which(data_inputs$Inputs=='a_11'),2]))
  
  Chi2_diff_test <-unname(   unlist(data_inputs[ which(data_inputs$Inputs=='Chi-square difference test'),2]))
  Chi2_diff_test_MC <-unname(   unlist(data_inputs[ which(data_inputs$Inputs=='Montecarlo Chi2 difference test'),2]))
  
  #There are 4 parameters of true preferences:
  #a_00, a_01, a_10, a_11. Any parameter can be fixed to a value in the range [0:1]
  
  user_restrictions <- c(a_00_restriction,
                         a_01_restriction,
                         a_10_restriction,
                         a_11_restriction)
  
  Montecarlo_P_value <-unname(   unlist(data_inputs[ which(data_inputs$Inputs=='Montecarlo P-value'),2]))
  Bootstrapped_CI  <-unname(   unlist(data_inputs[ which(data_inputs$Inputs=='Bootstrapped CI'),2]))
  Number_of_simulations <-as.numeric(unname(   unlist(data_inputs[ which(data_inputs$Inputs=='Number of simulations'),2])))
  Seed_number <-as.numeric(unname(   unlist(data_inputs[ which(data_inputs$Inputs=='Seed number'),2])))
  Output_file <-unname(   unlist(data_inputs[ which(data_inputs$Inputs=='Output file name'),2]))
  
  
  
  fit_method <-unname(   unlist(data_inputs[ which(data_inputs$Inputs=='CHI2 / G2'),2]))
  
  
  
  
  # Unless a restriction is set, the program by default estimates 4 true probabilitites
  if(is.null(user_restrictions)) { user_restrictions <-c(a_00="FREE", a_01="FREE",a_10="FREE", a_11="FREE") }
  
  suppressWarnings(check_restrictions <- as.numeric(ifelse(user_restrictions=="FREE", 0, user_restrictions)))
  
  if (length(check_restrictions)!=4 |
      any(check_restrictions<0)   | 
      any(user_restrictions=="FREE")==0 |
      any(is.na(check_restrictions))==1 |
      any(check_restrictions>=1) | 
      sum(check_restrictions)>=1 | 
      (Model_TE_4!="YES" & Model_TE_4!="NO") |
      (Model_TE_2!="YES" & Model_TE_2!="NO") |
      (Model_TE_1!="YES" & Model_TE_1!="NO") |
      (Model_Independence!="YES" & Model_Independence!="NO") |
      any(is.na(check_restrictions))) {
    cat("INAPPROPRIATE CONSTRAINED PARAMETERS! 
        
        \nEither a parameter was restricted to be negative or equal or greater than 1. 
        \nOr the sum of some restricted parameters adds to a value greater than 1.
        \nOr there are no free parameters to estimate.
        \nOr restrictions does not containt 4 elements.
        \nOr TE models to estimate have not been selected.
        \nOr there was a typo. \"FREE\" is required to be capitalized")
  } else if ( (sum(user_restrictions=="FREE")==4  & Chi2_diff_test=="YES") |
              (sum(user_restrictions=="FREE")==4 & Chi2_diff_test_MC=="YES") ) {
    cat("INAPPROPRIATE CONSTRAINED PARAMETERS! 
        
        \nChi2 difference test (either the standard test or the Monte Carlo version)
        \nrequires a constrained model. But the parameters set in the input file (Section 2) are all FREE.
        \nFor the test, you need to constrain one or two parameters in input file.
        \nOtherwise, mark the tests: Chi-square difference test (Section 3) and 
        \nthe Montecarlo Chi2 (or G2) difference test (Section 4) with \"NO\"")
  } else {
    
    
    # ESTIMATE PARAMETERS AND FIT THE MODELS
    #
    # The next lines fit the TE model to the data
    results_TE_1_errors <-NULL
    results_TE_2_errors <-NULL
    results_TE_4_errors <-NULL
    results_Independence <- NULL
    
    if (Model_TE_4=="YES") {
      results_TE_4_errors <- TE_estimation_data(data=response_patterns,
                                                restrictions= user_restrictions, 
                                                equal.prob.of.error.in.choice="TE-4", fit_method=fit_method)
    }
    if (Model_TE_2=="YES") {
      results_TE_2_errors <- TE_estimation_data(data=response_patterns,
                                                restrictions= user_restrictions, 
                                                equal.prob.of.error.in.choice="TE-2", fit_method=fit_method)
    }
    if (Model_TE_1=="YES") {
      results_TE_1_errors <- TE_estimation_data(data=response_patterns,
                                                restrictions= user_restrictions, 
                                                equal.prob.of.error.in.choice="TE-1", fit_method=fit_method)
    }
    if (Model_Independence=="YES") {
      results_Independence <- Independence_estimation_data(data=response_patterns, fit_method=fit_method)
    }
    
    
    
    if (Model_TE_4=="YES") {
      all_estimates_TE_4_errors <-aggregate_estimates(results=results_TE_4_errors, model="TE", fit_method = fit_method)
      file_name <- paste0("parameter estimates of TE - 4 errors, fit method ", fit_method, ".csv")
      write.csv(all_estimates_TE_4_errors, file=file_name )
    }
    
    if (Model_TE_2=="YES") {
      all_estimates_TE_2_errors <-aggregate_estimates(results=results_TE_2_errors, model="TE", fit_method = fit_method)
      file_name <- paste0("parameter estimates of TE - 2 errors, fit method ", fit_method, ".csv")
      write.csv(all_estimates_TE_2_errors, file=file_name)
    }
    
    if (Model_TE_1=="YES") {
      all_estimates_TE_1_errors <-aggregate_estimates(results=results_TE_1_errors, model="TE", fit_method = fit_method)
      file_name <- paste0("parameter estimates of TE - 1 error, fit method ", fit_method, ".csv")
      
      write.csv(all_estimates_TE_1_errors, file=file_name)
    }
    
    if (Model_Independence=="YES") { 
      all_estimates_Independence <- aggregate_estimates(results=results_Independence, model="Independence", fit_method = fit_method)
      file_name <- paste0("predictions of Independence, fit method ", fit_method, ".csv")
      
      write.csv(all_estimates_Independence, file=file_name)
    }
    
    
    
    
    
    # MONTE CARLO SECTION
    # ----------------------
    #
    #   We declare the seed for the random generator process prior to any estimation
    
    set.seed(Seed_number)    
    
    results_TE_4_MC <-NULL
    results_TE_2_MC <-NULL
    results_TE_1_MC <-NULL
    
    if (Chi2_diff_test_MC=="YES") {
      Montecarlo_P_value="YES"
      Chi2_diff_test="YES"
    }
    
    
    if (Montecarlo_P_value=="YES") {
      models_monte_carlo<-c(Model_TE_4,Model_TE_2,Model_TE_1)
      model_names <- c("TE-4", "TE-2", "TE-1")
      
      
      models_monte_carlo<-  ifelse(models_monte_carlo=="YES",model_names, "NO")
    }  else {
      models_monte_carlo<-c("NO", "NO", "NO")
      
    }
    
    for (val in models_monte_carlo) {
      if(val !="NO")  {
        results_TE_MC <- TE_MonteCarlo_data(data=response_patterns,
                                            restrictions=user_restrictions,
                                            equal.prob.of.error.in.choice=val,
                                            number_samples=Number_of_simulations,
                                            fit_method=fit_method) 
        
        
        all_p.values_MC_TE <- aggregate_MonteCarlo_p.value(results=results_TE_MC)
        #
        file_MC<- paste0("Monte Carlo p values, Model ", val,", fit method ", fit_method,".csv")
        #       Write out the p-values into a file
        write.csv(all_p.values_MC_TE, file= file_MC)
        
        if (val=="TE-4") {
          results_TE_4_MC<-results_TE_MC
        } else if (val=="TE-2") {
          results_TE_2_MC<-results_TE_MC
        } else {
          results_TE_1_MC<-results_TE_MC
        }
        
      }
    }
    
    # BOOTSTRAPPING SECTION
    # ---------------------
    
    results_TE_4_Boot_CI <-NULL
    results_TE_2_Boot_CI <-NULL
    results_TE_1_Boot_CI <-NULL
    
    if (Bootstrapped_CI=="YES") {
      models_boot<-c(Model_TE_4,Model_TE_2,Model_TE_1)
      model_names <- c("TE-4", "TE-2", "TE-1")
      models_boot<-  ifelse(models_boot=="YES",model_names, "NO")
    }  else {
      models_boot<-c("NO", "NO", "NO")
      
    }
    
    for (val in models_boot) {
      if(val !="NO")  {
        
        
        results.Boot <- TE_Bootstrapping.CI_data (data=response_patterns, 
                                                  restrictions=user_restrictions, 
                                                  equal.prob.of.error.in.choice=val,
                                                  number_samples=Number_of_simulations,
                                                  fit_method=fit_method) 
        all_boot_CI <- aggregate_boot.confidence.intervals(results=results.Boot)
        
        
        
        file_boot<- paste0("Bootstrapped confidence intervals, Model ", val,  ", fit method ", fit_method,".csv")
        write.csv(all_boot_CI, file=file_boot) 
        
        if (val=="TE-4") {
          results_TE_4_Boot_CI<-results.Boot
        } else if (val=="TE-2") {
          results_TE_2_Boot_CI<-results.Boot
        } else {
          results_TE_1_Boot_CI<-results.Boot
        }
        
      }
    }
    
    list(results_TE_1_errors=results_TE_1_errors,
         results_TE_2_errors=results_TE_2_errors,
         results_TE_4_errors=results_TE_4_errors, 
         results_Independence=results_Independence, 
         results_TE_4_MC=results_TE_4_MC, 
         results_TE_2_MC=results_TE_2_MC, 
         results_TE_1_MC=results_TE_1_MC,
         results_TE_4_Boot_CI=results_TE_4_Boot_CI, 
         results_TE_2_Boot_CI=results_TE_2_Boot_CI, 
         results_TE_1_Boot_CI=results_TE_1_Boot_CI,
         Chi2_diff_test_MC=Chi2_diff_test_MC,
         Chi2_diff_test=Chi2_diff_test,
         user_restrictions=user_restrictions,
         number_cases=number_cases,
         Output_file=Output_file,
         fit_method=fit_method)
  } 
}

#///////////////////////////////////////////////////////////////////////////////////////

# +++ Function "TE_READ_DATA"

TE_READ_DATA <- function(data_file) {
  simulation_results <-TE_read_inputs_and_MC(data_file)
  fit_method=simulation_results$fit_method
  
  if ( is.null(simulation_results$results_TE_4_MC) &
       is.null(simulation_results$results_TE_2_MC) &
       is.null(simulation_results$results_TE_1_MC) &
       
       simulation_results$Chi2_diff_test=="NO" & simulation_results$Chi2_diff_test_MC=="NO") {
    
    simulation_results  
  } else {
    
    if (simulation_results$Chi2_diff_test=="YES" & simulation_results$Chi2_diff_test_MC!="YES") {
      
      
      results_Chi2_diff_table <- test_chi2_diff_data(simulation_results=simulation_results) 
      
      n_elements <- length(results_Chi2_diff_table)
      
      Standard.p.value_data_TE_1 <- data.frame()
      Standard.p.value_data_TE_2 <- data.frame()
      Standard.p.value_data_TE_4 <- data.frame()
      unrestricted_all_TE_4 <- list()
      unrestricted_all_TE_2 <- list()
      unrestricted_all_TE_1 <- list()
      
      for (i in 1:n_elements) {
        
        Standard.p.value_results_TE_1 <- as.data.frame(  results_Chi2_diff_table[[i]]$TE_1_Chi2_diff_test$Chi2_diff_standard   )
        Standard.p.value_data_TE_1 <- rbind(Standard.p.value_data_TE_1, Standard.p.value_results_TE_1)
        
        Standard.p.value_results_TE_2 <- as.data.frame(  results_Chi2_diff_table[[i]]$TE_2_Chi2_diff_test$Chi2_diff_standard   )
        Standard.p.value_data_TE_2 <- rbind(Standard.p.value_data_TE_2, Standard.p.value_results_TE_2)
        
        Standard.p.value_results_TE_4 <- as.data.frame(  results_Chi2_diff_table[[i]]$TE_4_Chi2_diff_test$Chi2_diff_standard   )
        Standard.p.value_data_TE_4 <- rbind(Standard.p.value_data_TE_4, Standard.p.value_results_TE_4)
        
        unrestricted_all_TE_4[[i]] <-results_Chi2_diff_table[[i]]$TE_4_unrestricted_results
        unrestricted_all_TE_2[[i]] <-results_Chi2_diff_table[[i]]$TE_2_unrestricted_results
        unrestricted_all_TE_1[[i]] <-results_Chi2_diff_table[[i]]$TE_1_unrestricted_results
        
      }
      
      
      if (length(Standard.p.value_data_TE_1)>0) { 
        
        
        file_name <- paste0("Chi2 difference test, p values, Model TE-1, fit method ",fit_method ,".csv")
        write.csv( Standard.p.value_data_TE_1, file= file_name)
        
        restricted_model_estimates_TE_1 <- aggregate_estimates(results=simulation_results$results_TE_1_errors, model="TE", fit_method = fit_method)
        unrestricted_model_estimates_TE_1 <- aggregate_estimates(results=unrestricted_all_TE_1, model="TE", fit_method = fit_method)
        
        df1 <- cbind(Model="RESTRICTED MODEL ESTIMATES",restricted_model_estimates_TE_1 )
        df2 <- cbind(Model="UNRESTRICTED MODEL ESTIMATES", unrestricted_model_estimates_TE_1)  
        df1[setdiff(names(df2), names(df1)) ] <- NA
        df2[setdiff(names(df1), names(df2)) ] <- NA
        
        file_name <- paste0("Restricted and unrestricted models, estimates, Model TE-1, fit method ", fit_method , ".csv")
        write.csv( rbind(df1,df2), file= file_name)
        
      }
      
      if (length(Standard.p.value_data_TE_2)>0) { 
        
        file_name <- paste0("Chi2 difference test, p values, Model TE-2, fit method ",fit_method ,".csv")
        
        write.csv(Standard.p.value_data_TE_2, file= file_name)
        
        restricted_model_estimates_TE_2 <- aggregate_estimates(results=simulation_results$results_TE_2_errors, model="TE", fit_method = fit_method)
        unrestricted_model_estimates_TE_2 <- aggregate_estimates(results=unrestricted_all_TE_2, model="TE", fit_method = fit_method)
        
        df1 <- cbind(Model="RESTRICTED MODEL ESTIMATES",restricted_model_estimates_TE_2 )
        df2 <- cbind(Model="UNRESTRICTED MODEL ESTIMATES", unrestricted_model_estimates_TE_2)  
        df1[setdiff(names(df2), names(df1)) ] <- NA
        df2[setdiff(names(df1), names(df2)) ] <- NA
        
        file_name <- paste0("Restricted and unrestricted models, estimates, Model TE-2, fit method ", fit_method , ".csv")
        
        write.csv( rbind(df1,df2), file= file_name)
        
      }    
      
      if (length(Standard.p.value_data_TE_4)>0) { 
        
        file_name <- paste0("Chi2 difference test, p values, Model TE-4, fit method ",fit_method ,".csv")
        
        write.csv( Standard.p.value_data_TE_4, file= file_name)
        
        restricted_model_estimates_TE_4 <- aggregate_estimates(results=simulation_results$results_TE_4_errors, model="TE", fit_method = fit_method)
        unrestricted_model_estimates_TE_4 <- aggregate_estimates(results=unrestricted_all_TE_4, model="TE", fit_method = fit_method)
        
        df1 <- cbind(Model="RESTRICTED MODEL ESTIMATES",restricted_model_estimates_TE_4 )
        df2 <- cbind(Model="UNRESTRICTED MODEL ESTIMATES", unrestricted_model_estimates_TE_4)  
        df1[setdiff(names(df2), names(df1)) ] <- NA
        df2[setdiff(names(df1), names(df2)) ] <- NA
        
        file_name <- paste0("Restricted and unrestricted models, estimates, Model TE-4, fit method ", fit_method , ".csv")
        
        
        write.csv( rbind(df1,df2), file= file_name)
        
        
      }    
      
      
      list(results_TE_1_errors=simulation_results$results_TE_1_errors,
           results_TE_2_errors=simulation_results$results_TE_2_errors,
           results_TE_4_errors=simulation_results$results_TE_4_errors, 
           results_Independence=simulation_results$results_Independence, 
           results_TE_4_MC=simulation_results$results_TE_4_MC, 
           results_TE_2_MC=simulation_results$results_TE_2_MC, 
           results_TE_1_MC=simulation_results$results_TE_1_MC,
           results_TE_4_Boot_CI=simulation_results$results_TE_4_Boot_CI, 
           results_TE_2_Boot_CI=simulation_results$results_TE_2_Boot_CI, 
           results_TE_1_Boot_CI=simulation_results$results_TE_1_Boot_CI,
           user_restrictions=simulation_results$user_restrictions,
           fit_method=simulation_results$fit_method,
           results_Chi2_diff_table=results_Chi2_diff_table,
           Output_file = simulation_results$Output_file)
      
    } else if (simulation_results$Chi2_diff_test_MC=="YES") {
      
      results_Chi2_diff_MC<- MonteCarlo_test_chi2_diff_data(simulation_results=simulation_results) 
      
      n_elements <- length(results_Chi2_diff_MC)
      
      MonteCarlo.Re.fitted.p.value_data_TE_1 <- data.frame()
      MonteCarlo.Conservative.p.value_data_TE_1 <- data.frame()
      Standard.p.value_data_TE_1 <- data.frame()
      
      MonteCarlo.Re.fitted.p.value_data_TE_2 <- data.frame()
      MonteCarlo.Conservative.p.value_data_TE_2 <- data.frame()
      Standard.p.value_data_TE_2 <- data.frame()
      
      MonteCarlo.Re.fitted.p.value_data_TE_4 <- data.frame()
      MonteCarlo.Conservative.p.value_data_TE_4 <- data.frame()
      Standard.p.value_data_TE_4 <- data.frame()
      
      unrestricted_all_TE_4 <- list()
      unrestricted_all_TE_2 <- list()
      unrestricted_all_TE_1 <- list()
      
      
      for (i in 1:n_elements) {
        
        Standard.p.value_results_TE_1 <- as.data.frame(  results_Chi2_diff_MC[[i]]$TE_1_Chi2_diff_test$Chi2_diff_standard   )
        Standard.p.value_data_TE_1 <- rbind(Standard.p.value_data_TE_1, Standard.p.value_results_TE_1)
        
        MonteCarlo.Re.fitted.p.value_results_TE_1 <- as.data.frame(  results_Chi2_diff_MC[[i]]$TE_1_Chi2_diff_test$Chi2_diff_Re.fitted.MC$MC.p_value)                
        MonteCarlo.Re.fitted.p.value_data_TE_1 <- rbind(MonteCarlo.Re.fitted.p.value_data_TE_1, MonteCarlo.Re.fitted.p.value_results_TE_1)
        
        Standard.p.value_results_TE_2 <- as.data.frame(  results_Chi2_diff_MC[[i]]$TE_2_Chi2_diff_test$Chi2_diff_standard   )
        Standard.p.value_data_TE_2 <- rbind(Standard.p.value_data_TE_2, Standard.p.value_results_TE_2)
        
        MonteCarlo.Re.fitted.p.value_results_TE_2 <- as.data.frame(  results_Chi2_diff_MC[[i]]$TE_2_Chi2_diff_test$Chi2_diff_Re.fitted.MC$MC.p_value)                
        MonteCarlo.Re.fitted.p.value_data_TE_2 <- rbind(MonteCarlo.Re.fitted.p.value_data_TE_2, MonteCarlo.Re.fitted.p.value_results_TE_2)
        
        Standard.p.value_results_TE_4 <- as.data.frame(  results_Chi2_diff_MC[[i]]$TE_4_Chi2_diff_test$Chi2_diff_standard   )
        Standard.p.value_data_TE_4 <- rbind(Standard.p.value_data_TE_4, Standard.p.value_results_TE_4)
        
        MonteCarlo.Re.fitted.p.value_results_TE_4 <- as.data.frame(  results_Chi2_diff_MC[[i]]$TE_4_Chi2_diff_test$Chi2_diff_Re.fitted.MC$MC.p_value)                
        MonteCarlo.Re.fitted.p.value_data_TE_4 <- rbind(MonteCarlo.Re.fitted.p.value_data_TE_4, MonteCarlo.Re.fitted.p.value_results_TE_4)
        
        unrestricted_all_TE_4[[i]] <-results_Chi2_diff_MC[[i]]$TE_4_unrestricted_results
        unrestricted_all_TE_2[[i]] <-results_Chi2_diff_MC[[i]]$TE_2_unrestricted_results
        unrestricted_all_TE_1[[i]] <-results_Chi2_diff_MC[[i]]$TE_1_unrestricted_results
        
        
      }
      
      
      MonteCarlo.p.value_data_all_TE_1<- cbind(Standard.p.value_data_TE_1 ,  MonteCarlo.Re.fitted.p.value_data_TE_1 )
      MonteCarlo.p.value_data_all_TE_2<- cbind(Standard.p.value_data_TE_2 ,  MonteCarlo.Re.fitted.p.value_data_TE_2 )
      MonteCarlo.p.value_data_all_TE_4<- cbind(Standard.p.value_data_TE_4 ,  MonteCarlo.Re.fitted.p.value_data_TE_4 )
      
      try( {colnames(MonteCarlo.p.value_data_all_TE_1)[5] <- c("Re-fitted MC - p value")}, 
           silent = TRUE)
      try( {colnames(MonteCarlo.p.value_data_all_TE_2)[5] <- c("Re-fitted MC - p value")},
           silent=TRUE)
      try( {colnames(MonteCarlo.p.value_data_all_TE_4)[5] <- c( "Re-fitted MC - p value")},
           silent=TRUE)
      
      
      
      
      if (length(MonteCarlo.p.value_data_all_TE_1)>0) { 
        write.csv( MonteCarlo.p.value_data_all_TE_1, file= "Monte Carlo Chi2 difference test, p values, Model TE-1.csv")
        
        restricted_model_estimates_TE_1 <- aggregate_estimates(results=simulation_results$results_TE_1_errors, model="TE", fit_method = fit_method)
        unrestricted_model_estimates_TE_1 <- aggregate_estimates(results=unrestricted_all_TE_1, model="TE", fit_method = fit_method)
        
        df1 <- cbind(Model="RESTRICTED MODEL ESTIMATES",restricted_model_estimates_TE_1 )
        df2 <- cbind(Model="UNRESTRICTED MODEL ESTIMATES", unrestricted_model_estimates_TE_1)  
        df1[setdiff(names(df2), names(df1)) ] <- NA
        df2[setdiff(names(df1), names(df2)) ] <- NA
        file_name <- paste0("Restricted and unrestricted models, estimates, Model TE-1, fit method ", fit_method , ".csv")
        write.csv( rbind(df1,df2), file= file_name)      
        
      }
      
      if (length(MonteCarlo.p.value_data_all_TE_2)>0) { 
        write.csv( MonteCarlo.p.value_data_all_TE_2, file= "Monte Carlo Chi2 difference test, p values, Model TE-2.csv")
        
        restricted_model_estimates_TE_2 <- aggregate_estimates(results=simulation_results$results_TE_2_errors, model="TE", fit_method = fit_method)
        unrestricted_model_estimates_TE_2 <- aggregate_estimates(results=unrestricted_all_TE_2, model="TE", fit_method = fit_method)
        
        df1 <- cbind(Model="RESTRICTED MODEL ESTIMATES",restricted_model_estimates_TE_2 )
        df2 <- cbind(Model="UNRESTRICTED MODEL ESTIMATES", unrestricted_model_estimates_TE_2)  
        df1[setdiff(names(df2), names(df1)) ] <- NA
        df2[setdiff(names(df1), names(df2)) ] <- NA
        file_name <- paste0("Restricted and unrestricted models, estimates, Model TE-2, fit method ", fit_method , ".csv")
        write.csv( rbind(df1,df2), file= file_name)      
      }    
      
      if (length(MonteCarlo.p.value_data_all_TE_4)>0) { 
        write.csv( MonteCarlo.p.value_data_all_TE_4, file= "Monte Carlo Chi2 difference test, p values, Model TE-4.csv")
        
        restricted_model_estimates_TE_4 <- aggregate_estimates(results=simulation_results$results_TE_4_errors, model="TE", fit_method = fit_method)
        unrestricted_model_estimates_TE_4 <- aggregate_estimates(results=unrestricted_all_TE_4, model="TE", fit_method = fit_method)
        
        df1 <- cbind(Model="RESTRICTED MODEL ESTIMATES",restricted_model_estimates_TE_4 )
        df2 <- cbind(Model="UNRESTRICTED MODEL ESTIMATES", unrestricted_model_estimates_TE_4)  
        df1[setdiff(names(df2), names(df1)) ] <- NA
        df2[setdiff(names(df1), names(df2)) ] <- NA
        file_name <- paste0("Restricted and unrestricted models, estimates, Model TE-4, fit method ", fit_method , ".csv")
        write.csv( rbind(df1,df2), file= file_name)      
        
        
      }    
      
      
      
      
      
      
      
      list(results_TE_1_errors=simulation_results$results_TE_1_errors,
           results_TE_2_errors=simulation_results$results_TE_2_errors,
           results_TE_4_errors=simulation_results$results_TE_4_errors, 
           results_Independence=simulation_results$results_Independence, 
           results_TE_4_MC=simulation_results$results_TE_4_MC, 
           results_TE_2_MC=simulation_results$results_TE_2_MC, 
           results_TE_1_MC=simulation_results$results_TE_1_MC,
           results_TE_4_Boot_CI=simulation_results$results_TE_4_Boot_CI, 
           results_TE_2_Boot_CI=simulation_results$results_TE_2_Boot_CI, 
           results_TE_1_Boot_CI=simulation_results$results_TE_1_Boot_CI,
           user_restrictions=simulation_results$user_restrictions,
           fit_method=simulation_results$fit_method,
           results_Chi2_diff_MC=results_Chi2_diff_MC,
           Output_file = simulation_results$Output_file)
    } else {
      
      
      simulation_results
    }
    
  }
  
}





