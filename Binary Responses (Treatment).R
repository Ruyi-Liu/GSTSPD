library(gsDesign)
library("reshape2") 
library(ggplot2)
library(preference)
library(extraDistr)

# number of iterations to compute average GSD sample size
N_stopping <- 10000
# number of iterations to compute empirical Type I error and empirical power
N_error_rate <- 10000

# Preference Rate
phi = 0.5
# First-stage Randomization Ratio
theta = 0.5
# Treatment Difference
trt_diff = 0.1
# Selection Difference
select_diff = 0.1
# Preference Difference
prefer_diff = 0.1
# Overall Response Proportion
p = 0.5
mu = p
# Significance Level
alpha = 0.05
# Power
power = 0.9
# Maximum Number of Interim Analysis
K = 3
# Information Allocation
timing_info <- c(1)


##################### Treatment Effect ######################### 

# Generate dataset for tesing treatment effect
dataset_treat <- function(n,phi,theta,trt_diff,select_diff,prefer_diff,mu){
  # Generate Patient ID
  cID <- c(1:n)
  # Generate patient preference 
  cPreference <- rbinom(n, 1, phi)
  # Generate Dataset
  df <- data.frame(ID = cID, Preference=cPreference)
  
  # Generate patient treatment (1: treatment A; 0: treatment B)
  # First Stage
  option_index <- sample(nrow(df), theta*n)
  first_stage_randomization <- rep(0, n)
  for (i in 1:n){
    if (is.element(i,option_index)){
      first_stage_randomization[i] <- 1
    }
  }
  df$first_stage <- factor(first_stage_randomization, levels=c(0,1),
                           labels=c("Random","Option"))
  
  # Second Stage  
  df$second_stage[df$first_stage == "Option"] <- df$Preference[df$first_stage == "Option"]
  Random_subset_ID <- df[df$first_stage == "Random",]$ID
  sec_trtA_index <- sample(Random_subset_ID, (1-theta)*n/2)
  df$second_stage[df$first_stage == "Random"] <- 
    ifelse(is.element(df[df$first_stage == "Random",]$ID,sec_trtA_index),1,0)
  df$Treatment <- df$second_stage
  
  # Generate Response Variable 
  tau_1<- trt_diff/2
  tau_2<- -trt_diff/2
  nu_1<- (1-phi)*select_diff
  nu_2<- -phi*select_diff
  pi11<- (1-phi)*prefer_diff
  pi21<- -(1-phi)*prefer_diff
  pi12<- -phi*prefer_diff
  pi22<- phi*prefer_diff
  p = mu 
  prob_11 <- p + tau_1 + nu_1 + pi11 
  prob_12 <- p + tau_1 + nu_2 + pi12
  prob_21 <- p + tau_2 + nu_1 + pi21 
  prob_22 <- p + tau_2 + nu_2 + pi22
  prob1 <- phi*prob_11 + (1-phi)*prob_12
  prob2 <- phi*prob_21 + (1-phi)*prob_22
  
  x_1 <- rbern(nrow(df[df$first_stage == "Option" & df$second_stage == 1,]),prob_11)
  
  x_2 <- rbern(nrow(df[df$first_stage == "Option" & df$second_stage == 0,]),prob_22)
  
  y_1 <- rbern(nrow(df[df$first_stage == "Random" & df$second_stage == 1,]),prob1)
  
  y_2 <- rbern(nrow(df[df$first_stage == "Random" & df$second_stage == 0,]),prob2)
  
  df$Y[df$first_stage == "Option" & df$second_stage == 1] <- x_1
  
  df$Y[df$first_stage == "Option" & df$second_stage == 0] <- x_2
  
  df$Y[df$first_stage == "Random" & df$second_stage == 1] <- y_1
  
  df$Y[df$first_stage == "Random" & df$second_stage == 0] <- y_2 
  
  list_return <- list()
  list_return[[1]] <- df
  list_return[[2]] <- c(prob_11,prob_22,prob1,prob2)
  return(list_return)
}

# Calculate test statistics  
cum_score_treat <- function(dataset_cum, theta,q_cum,K,N,prob_lst){
  tau_lst <- c()
  for (i in 1:K){
    df <- dataset_cum[[i]]
    n_1 <- length(df$Y[df$Treatment == 1 & df$first_stage == "Random"]) 
    n_2 <- length(df$Y[df$Treatment == 0 & df$first_stage == "Random"]) 
    
    s1_sq <- ifelse(n_1 >1, mean(df$Y[df$Treatment == 1 & df$first_stage == "Random"])*
                      (1-mean(df$Y[df$Treatment == 1 & df$first_stage == "Random"])),0)
    
    s2_sq <- ifelse(n_2 >1, mean(df$Y[df$Treatment == 0 & df$first_stage == "Random"])*
                      (1-mean(df$Y[df$Treatment == 0 & df$first_stage == "Random"])),0)
    
    tau_val <- (mean(df$Y[df$Treatment == 1 & df$first_stage == "Random"]) -
                  mean(df$Y[df$Treatment == 0 &
                              df$first_stage == "Random"]))/sqrt((s1_sq/n_1) + (s2_sq/n_2))

    tau_lst[i] <- tau_val
  }
  return(tau_lst)
}

# Calculate empirical Type I error rate and empirical power for GSD
error_rate <- function(num_loops,upper_bound,lower_bound,N,phi,theta,trt_diff,select_diff,
                       prefer_diff,mu,q,K,q_cum){
  
  reject_lst <- c()
  for (n_loop in 1:num_loops){
    dataset_lst <- list()
    # Generate Datasets
    for (i in 1:K){
      dataset_lst_ret <- dataset_treat(ceiling(q[i]*N),phi,theta,trt_diff,select_diff,prefer_diff,mu)
      dataset_lst[[i]] <- dataset_lst_ret[[1]]
      prob_lst<- dataset_lst_ret[[2]]
    }
    dataset_cum <- dataset_lst
    for (i in 2:K){
      dataset_cum[[i]] <- rbind(dataset_cum[[i-1]], dataset_cum[[i]])
    }
    # Calculate test statistics
    test_stats_treat <- cum_score_treat(dataset_cum,theta,q_cum,K,N,prob_lst)
    for (i in 1:K){
      if(test_stats_treat[i] < lower_bound[i]|test_stats_treat[i] > upper_bound[i]){
        reject_lst[n_loop] <- 1
        break
      }
      if (i == K){
        reject_lst[n_loop] <- 0
      }
    }
  }
  return(mean(reject_lst))
}

# Calculate fixed design empirical Type I error rate and empirical power
error_rate_fixed_treat <- function(num_loops,N,phi,theta,trt_diff,select_diff,
                                   prefer_diff,mu,alpha_val){
  
  reject_lst <- c()
  for (n_loop in 1:num_loops){
    # Generate Datasets
    dataset_lst_ret <- dataset_treat(N,phi,theta,trt_diff,select_diff,prefer_diff,mu)
    df <- dataset_lst_ret[[1]]
    prob_lst<- dataset_lst_ret[[2]]
    n_1 <- length(df$Y[df$Treatment == 1 & df$first_stage == "Random"]) 
    n_2 <- length(df$Y[df$Treatment == 0 & df$first_stage == "Random"]) 
    var_y_1 <- prob_lst[3]*(1-prob_lst[3])
    var_y_2 <- prob_lst[4]*(1-prob_lst[4])
    
    test_stats <- (mean(df$Y[df$Treatment == 1 & df$first_stage == "Random"]) -
                     mean(df$Y[df$Treatment == 0 &
                                 df$first_stage == "Random"]))/sqrt((var_y_1/n_1) + (var_y_2/n_2))
    upper_bound <- qnorm(1-alpha_val/2,mean=0,sd=1)
    lower_bound <- qnorm(alpha_val/2,mean=0,sd=1)
    
    if (test_stats < lower_bound | test_stats > upper_bound){
      reject_lst[n_loop] <- 1
    }
    else{
      reject_lst[n_loop] <- 0
    }
  }
  return(mean(reject_lst))
}

# GSD sample size
stopping_sample_size <- function(num_loops,phi,theta,trt_diff,select_diff,
                                 prefer_diff,mu,K,sfu_type,timing_info,N_error_rate){
  delta1 <- trt_diff
  delta0 <- 0
  ratio <- 1
  sd <- 1
  sd2 <- 1
  # Number of analyses
  num_k <- K
  test.type <- 2
  alpha <- 0.025
  beta <- 0.1
  astar <- 0
  timing <- timing_info
  sfu <- sfu_type
  sfupar <- c(0)
  sfl <- sfLDOF
  sflpar <- c(0)
  delta <- 0
  endpoint <- "Binomial"
  alpha_overall <- 0.05
  
  tau_1<- trt_diff/2
  tau_2<- -trt_diff/2
  nu_1<- (1-phi)*select_diff
  nu_2<- -phi*select_diff
  pi11<- (1-phi)*prefer_diff
  pi21<- -(1-phi)*prefer_diff
  pi12<- -phi*prefer_diff
  pi22<- phi*prefer_diff
  p = mu
  prob_11 <- p + tau_1 + nu_1 + pi11
  prob_12 <- p + tau_1 + nu_2 + pi12
  prob_21 <- p + tau_2 + nu_1 + pi21
  prob_22 <- p + tau_2 + nu_2 + pi22
  prob1 <- phi*prob_11 + (1-phi)*prob_12
  prob2 <- phi*prob_21 + (1-phi)*prob_22
  d1 <- prob_11 - prob1
  d2 <- prob_22 - prob2
  
  # Fixed design sample size
  n_tau = 2*ceiling(((qnorm(1-beta,mean=0,sd=1)+qnorm(1-alpha_overall/2,mean=0,sd=1))^2)*
                      (prob1*(1-prob1) + prob2*(1-prob2))/((1-theta)*trt_diff^2))
  n <- ceiling(n_tau)
  
  x <- gsDesign(
    k = num_k,
    test.type = test.type,
    alpha = alpha,
    beta = beta,
    astar = astar,
    timing = timing,
    sfu = sfu,
    sfupar = sfupar,
    sfl = sfl,
    sflpar = sflpar,
    delta = delta,
    delta1 = delta1,
    delta0 = delta0,
    endpoint = endpoint,
    n.fix = n
  )
  n_lst <- ceiling(x$n.I)
  # Maximum GSD sample size
  N <- n_lst[length(x$n.I)]
  
  # information fraction
  q <- c()
  q0 <- 0
  for (i in 1:length(n_lst)){
    q[i] <- (n_lst[i] - q0)/N
    q0<- n_lst[i]
  }
  q<- round(q,5)
  q_cum <- c()
  for (i in 1:K){
    q_cum[i] <- sum(q[1:i])
  }
  
  ret_stopping_size <- c()
  for (n_loop in 1:num_loops){
    dataset_lst <- list()
    # Generate Datasets
    for (i in 1:K){
      dataset_lst_ret <- dataset_treat(ceiling(q[i]*N),phi,theta,trt_diff,select_diff,prefer_diff,mu)
      dataset_lst[[i]] <- dataset_lst_ret[[1]]
      prob_lst<- dataset_lst_ret[[2]]
    }
    dataset_cum <- dataset_lst
    for (i in 2:K){
      dataset_cum[[i]] <- rbind(dataset_cum[[i-1]], dataset_cum[[i]])
    }
    # Calculate test statistics
    test_stats_treat <- cum_score_treat(dataset_cum,theta,q_cum,K,N,prob_lst)
    for (i in 1:K){
      if(test_stats_treat[i] < x$lower$bound[i]|test_stats_treat[i] > x$upper$bound[i]){
        ret_stopping_size[n_loop] <- N*q_cum[i]
        break
      }
      if (i == K){
        ret_stopping_size[n_loop] <- N
      }
    }
  }
  
  # GSD empirical Type I error rate
  type1_error <- error_rate(N_error_rate,x$upper$bound,
                            x$lower$bound,N,phi,theta,0,select_diff,
                            prefer_diff,mu,q,K,q_cum)
  
  # GSD empirical power
  GSD_power <- error_rate(N_error_rate,x$upper$bound,
                          x$lower$bound,N,phi,theta,trt_diff,select_diff,
                          prefer_diff,mu,q,K,q_cum)
  
  if (sfu_type == "OF"){
    # Fixed design empirical Type I error rate
    fixed_type_1 <- error_rate_fixed_treat(N_error_rate,n,phi,theta,0,select_diff,
                                           prefer_diff,mu,alpha_overall)
    
    # Fixed design empirical power
    fixed_power <- error_rate_fixed_treat(N_error_rate,n,phi,theta,trt_diff,select_diff,
                                          prefer_diff,mu,alpha_overall)
  }
  else{
    fixed_type_1 <- 0
    fixed_power <- 0
  }
  return(list(ret_stopping_size,n,type1_error,GSD_power,fixed_type_1,fixed_power))
}

Treatment_Difference <- seq(0.1,0.3,0.05)
Sequential_Sample_Size_Pocock <- c()
Sequential_Sample_Size_OF <- c()
Fixed_Sample_Size <- c()
sd_lst_pocock <- c()
sd_lst_OBF <- c()
type1_pocock <- c()
type1_OBF <- c()
power_pocock <- c()
power_OBF <- c()
fixed_type_1 <- c()
fixed_power <- c()

for (i in 1:length(Treatment_Difference)){
  
  # Pocock result
  ret_pocock <- stopping_sample_size(N_stopping,phi,theta,
                                     Treatment_Difference[i],
                                     select_diff,prefer_diff,mu,K,"Pocock",timing_info,N_error_rate)
  
  Sequential_Sample_Size_Pocock[i] <- ceiling(mean(ret_pocock[[1]]))
  Fixed_Sample_Size[i] <- ceiling(ret_pocock[[2]])
  sd_lst_pocock[i]<- round(sqrt(var(ret_pocock[[1]])),2)
  type1_pocock[i] <- ret_pocock[[3]]
  power_pocock[i] <- ret_pocock[[4]]
  
  # OBF result
  ret_OBF <- stopping_sample_size(N_stopping,phi,theta,
                                  Treatment_Difference[i],
                                  select_diff,prefer_diff,mu,K,"OF",timing_info,N_error_rate)
  
  Sequential_Sample_Size_OF[i] <- ceiling(mean(ret_OBF[[1]]))
  sd_lst_OBF[i]<- round(sqrt(var(ret_OBF[[1]])),2)
  type1_OBF[i] <- ret_OBF[[3]]
  power_OBF[i] <- ret_OBF[[4]]
  fixed_type_1[i] <- ret_OBF[[5]]
  fixed_power[i] <- ret_OBF[[6]]
}

df <- cbind(Treatment_Difference,Sequential_Sample_Size_Pocock,Sequential_Sample_Size_OF,
            Fixed_Sample_Size,sd_lst_pocock,type1_pocock,power_pocock,sd_lst_OBF,type1_OBF,power_OBF,
            fixed_type_1,fixed_power)

df <- data.frame(df)
colnames(df) <- c('Treatment Difference', 'Pocock Sample Size', 'OBF Sample Size',
                  'Fixed Sample Size','Pocock SD', 'Pocock Type 1', 'Pocock Power','OBF SD','OBF Type 1',
                  'OBF Power','Fixed Type 1','Fixed Power')
View(df)
