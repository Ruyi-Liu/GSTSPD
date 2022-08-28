library(gsDesign)
library("reshape2") 
library(ggplot2)
library(preference)
library(ggpubr)

# Preference Rate
phi = 0.5
# First-stage Randomization Ratio
theta = 0.5
# Treatment Difference
trt_diff = 0.4
# Selection Difference
select_diff = 0.2
# Preference Difference
prefer_diff = 0.3
# Overall Mean Response
mu = 0
# Significance Level
alpha = 0.05
# Power
power = 0.9
# Maximum Number of Interim Analysis
K = 3

# Number of iterations
sample_size_num_loop = 5000

# Information Allocation
timing_info <- c(1)


# Function for generate datasets
dataset <- function(n,phi,theta,trt_diff,select_diff,prefer_diff,mu){
  # Generate Patient ID
  n <- ceiling(n)
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
  
  df$Y[df$Treatment == 1 & df$Preference == 1] <- mu + tau_1 + nu_1 + pi11 +
    rnorm(nrow(df[df$Treatment == 1 & df$Preference == 1,]),mean=0,sd=1)
  
  df$Y[df$Treatment == 1 & df$Preference == 0] <- mu + tau_1 + nu_2 + pi12 +
    rnorm(nrow(df[df$Treatment == 1 & df$Preference == 0,]),mean=0,sd=1)
  
  df$Y[df$Treatment == 0 & df$Preference == 1] <- mu + tau_2 + nu_1 + pi21 +
    rnorm(nrow(df[df$Treatment == 0 & df$Preference == 1,]),mean=0,sd=1)
  
  df$Y[df$Treatment == 0 & df$Preference == 0] <- mu + tau_2 + nu_2 + pi22 +
    rnorm(nrow(df[df$Treatment == 0 & df$Preference == 0,]),mean=0,sd=1)
  
  return(df)
}




##################### Treatment effect ######################### 

## Calculate test statistics
cum_score_treat <- function(dataset_cum, theta,q_cum,K,N){
  tau_lst <- c()
  for (i in 1:K){
    df <- dataset_cum[[i]]
    tau_val <- (mean(df$Y[df$Treatment == 1 & df$first_stage == "Random"]) -
                  mean(df$Y[df$Treatment == 0 &
                              df$first_stage == "Random"]))/sqrt(2*1/((1-theta)*q_cum[i]*N/2))
    #print((1-theta)*q_cum[i]*N/2)
    #tau_lst[i] <- tau_val*sqrt(i)
    tau_lst[i] <- tau_val
  }
  
  return(tau_lst)
}

# Function for computing empirical Type I error and empirical power
error_rate <- function(num_loops,upper_bound,lower_bound,N,phi,theta,trt_diff,select_diff,
                       prefer_diff,mu,q,K,q_cum){
  
  reject_lst <- c()
  for (n_loop in 1:num_loops){
    dataset_lst <- list()
    # Generate Datasets
    for (i in 1:K){
      dataset_lst[[i]] <- dataset(ceiling(q[i]*N),phi,theta,trt_diff,select_diff,prefer_diff,mu)
    }
    dataset_cum <- dataset_lst
    for (i in 2:K){
      dataset_cum[[i]] <- rbind(dataset_cum[[i-1]], dataset_cum[[i]])
    }
    # Calculate test statistics
    test_stats_treat <- cum_score_treat(dataset_cum,theta,q_cum,K,N)
    for (i in 1:K){
      # Compare with boundary statistics
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



# Compute average group sequential sample size
stopping_sample_size <- function(num_loops,phi,theta,trt_diff,select_diff,
                                 prefer_diff,mu,K,sfu_type,timing_info){
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
  endpoint <- 'normal' 
  
  n <- nNormal(
    delta1 = delta1,
    delta0 = delta0,
    ratio = ratio,
    sd = sd,
    sd2 = sd2,
    alpha = alpha,
    beta = beta
  )
  
  # Fixed design sample size
  n <- ceiling(n/(1-theta))
  
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
  # Maximum group sequential sample size
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
      dataset_lst[[i]] <- dataset(ceiling(q[i]*N),phi,theta,trt_diff,select_diff,prefer_diff,mu)
    }
    dataset_cum <- dataset_lst
    for (i in 2:K){
      dataset_cum[[i]] <- rbind(dataset_cum[[i-1]], dataset_cum[[i]])
    }
    # Calculate test statistics
    test_stats_treat <- cum_score_treat(dataset_cum,theta,q_cum,K,N)
    for (i in 1:K){
      # Compare with boundary statistics
      if(test_stats_treat[i] < x$lower$bound[i]|test_stats_treat[i] > x$upper$bound[i]){
        ret_stopping_size[n_loop] <- N*q_cum[i]
        break
      }
      if (i == K){
        ret_stopping_size[n_loop] <- N
      }
    }
    
  }
  
  # Compute empirical Type I error
  type1_error <- error_rate(2,x$upper$bound,
                            x$lower$bound,N,phi,theta,0,select_diff,
                            prefer_diff,mu,q,K,q_cum)
  
  # Compute empirical power
  emp_power <- error_rate(2,x$upper$bound,
                            x$lower$bound,N,phi,theta,trt_diff,select_diff,
                            prefer_diff,mu,q,K,q_cum)
  
  return(list(ret_stopping_size,n,type1_error,emp_power))
}

Treatment_Difference <- seq(0.2,1,0.1)
Sequential_Sample_Size_Pocock <- c()
Sequential_Sample_Size_OF <- c()
Fixed_Sample_Size <- c()
sd_lst_pocock <- c()
sd_lst_OBF <- c()
type1_pocock <- c()
type1_OBF <- c()
emp_power_pocock <- c()
emp_power_OBF <- c()


for (i in 1:length(Treatment_Difference)){
  # results given by Pocock boundary
  ret_pocock <- stopping_sample_size(sample_size_num_loop,phi,theta,
                                     Treatment_Difference[i],
                                     select_diff,prefer_diff,mu,K,"Pocock",timing_info)
  
  Sequential_Sample_Size_Pocock[i] <- ceiling(mean(ret_pocock[[1]]))
  Fixed_Sample_Size[i] <- ceiling(ret_pocock[[2]])
  sd_lst_pocock[i]<- round(sqrt(var(ret_pocock[[1]])),2)
  type1_pocock[i] <- ret_pocock[[3]]
  emp_power_pocock[i] <- ret_pocock[[4]]
  
  # results given by OBF boundary
  ret_OBF <- stopping_sample_size(sample_size_num_loop,phi,theta,
                                  Treatment_Difference[i],
                                  select_diff,prefer_diff,mu,K,"OF",timing_info)
  
  Sequential_Sample_Size_OF[i] <- ceiling(mean(ret_OBF[[1]]))
  sd_lst_OBF[i]<- round(sqrt(var(ret_OBF[[1]])),2)
  type1_OBF[i] <- ret_OBF[[3]]
  emp_power_OBF[i] <- ret_OBF[[4]]
}



df <- cbind(Treatment_Difference,Sequential_Sample_Size_Pocock,Sequential_Sample_Size_OF,
            Fixed_Sample_Size,sd_lst_pocock,type1_pocock,emp_power_pocock,sd_lst_OBF,
            type1_OBF,emp_power_OBF)

df <- data.frame(df)
colnames(df) <- c('Treatment Difference', 'Pocock Sample Size', 'OBF Sample Size',
                  'Fixed Sample Size','Pocock SD', 'Pocock Type 1','Pocock Power','OBF SD',
                  'OBF Type 1','OBF Power')
View(df)

values <- c(Sequential_Sample_Size_Pocock/Fixed_Sample_Size,Sequential_Sample_Size_OF/Fixed_Sample_Size)
group <- c(rep("Pocock/Fixed",length(Treatment_Difference)),rep("OBF/Fixed",length(Treatment_Difference)))
results <- data.frame(Treatment_Difference,group,values)
print(results)

# plot
g_treatment <- ggplot(data = results, aes(x=Treatment_Difference, y=values, group = group)) +geom_line(aes(color=group))+
  geom_point(aes(color=group))  + labs(title = "\n", x = "Treatment Effect",
                                       y = "Ratio")


################ Selection Difference ####################

## Calculate test statistics
test_select_diff <- function(dataset_cum,K){
  
  nu_lst <- c()
  for (i in 1:K){
    df <- dataset_cum[[i]]
    m_1 <- length(df$Y[df$Treatment == 1 & df$first_stage == "Option"]) 
    m_2 <- length(df$Y[df$Treatment == 0 & df$first_stage == "Option"]) 
    m = m_1 + m_2  
    n_1 <- length(df$Y[df$Treatment == 1 & df$first_stage == "Random"]) 
    n_2 <- length(df$Y[df$Treatment == 0 & df$first_stage == "Random"]) 
    
    # Empirical Variance
    s11_sq <- ifelse(m_1 > 1,var(df$Y[df$Treatment == 1 & df$first_stage == "Option"]),0 )
    s22_sq <- ifelse(m_2 > 1, var(df$Y[df$Treatment == 0 & df$first_stage == "Option"]),0)
    s1_sq <- ifelse(n_1 >1, var(df$Y[df$Treatment == 1 & df$first_stage == "Random"]),0)
    s2_sq <- ifelse(n_2 >1, var(df$Y[df$Treatment == 0 & df$first_stage == "Random"]),0)
    
    # Define z_1 and z_2
    z_1 <- sum(df$Y[df$Treatment == 1 & df$first_stage == "Option"]) -
      m_1 *mean(df$Y[df$Treatment == 1 & df$first_stage == "Random"])
    z_2 <- sum(df$Y[df$Treatment == 0 & df$first_stage == "Option"]) -
      m_2 *mean(df$Y[df$Treatment == 0 & df$first_stage == "Random"])
    
    # Calculate variance of z_1-z_2
    if (m_1 != 0 & m_2 != 0){
      var_z1 <- m_1*s11_sq+(m_1*s1_sq/n_1)*(1+m_1*(m-1)/m)+ (m_1*m_2/m)*
        (mean(df$Y[df$Treatment == 1 & df$first_stage == "Option"])-
           mean(df$Y[df$Treatment == 1 & df$first_stage == "Random"]))^2
      
      var_z2 <- m_2*s22_sq+(m_2*s2_sq/n_2)*(1+m_2*(m-1)/m)+ (m_1*m_2/m)*
        (mean(df$Y[df$Treatment == 0 & df$first_stage == "Option"])-
           mean(df$Y[df$Treatment == 0 & df$first_stage == "Random"]))^2
      
      cov_z1_z2 <- -(m_1*m_2/m)*(mean(df$Y[df$Treatment == 1 & df$first_stage == "Option"])-
                                   mean(df$Y[df$Treatment == 1 & df$first_stage == "Random"]))*
        (mean(df$Y[df$Treatment == 0 & df$first_stage == "Option"])-
           mean(df$Y[df$Treatment == 0 & df$first_stage == "Random"]))
      
    }
    
    if (m_1 == 0){
      var_z1 <- 0
      var_z2 <- m_2*s22_sq+(m_2*s2_sq/n_2)*(1+m_2*(m-1)/m)
      cov_z1_z2 <- 0
    }
    if (m_2 == 0){
      var_z1 <- m_1*s11_sq+(m_1*s1_sq/n_1)*(1+m_1*(m-1)/m)
      var_z2 <- 0
      cov_z1_z2 <- 0
    }
    var_z1z2 <- var_z1 + var_z2 - cov_z1_z2
    # Test statistics
    nu_lst[i]<- (z_1 - z_2)/sqrt(var_z1z2)
    
  }
  return(nu_lst)
}

# Function for computing empirical Type I error rate and empirical power
error_rate_select <- function(num_loops,upper_bound,lower_bound,N,phi,theta,trt_diff,
                              select_diff,prefer_diff,mu,q,K,q_cum){
  
  reject_lst <- c()
  for (n_loop in 1:num_loops){
    dataset_lst <- list()
    # Generate Datasets
    for (i in 1:K){
      dataset_lst[[i]] <- dataset(q[i]*N,phi,theta,trt_diff,select_diff,prefer_diff,mu)
    }
    dataset_cum <- dataset_lst
    for (i in 2:K){
      dataset_cum[[i]] <- rbind(dataset_cum[[i-1]], dataset_cum[[i]])
    }
    # Calculate test statistics
    test_stats <- test_select_diff(dataset_cum,K)
    for (i in 1:K){
      if(test_stats[i] < lower_bound[i]|test_stats[i] > upper_bound[i]){
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

# Average GSD sample size
stopping_sample_size_select <- function(num_loops,phi,theta,trt_diff,select_diff,
                                        prefer_diff,mu,K,sfu_type,timing_info){
  # Number of analyses
  k <- K
  # See gsDesign() help for description of boundary types
  test.type <- 2
  # 1-sided Type I error
  alpha <- 0.025
  # Type 2 error (1 - targeted power)
  beta <- 0.1
  astar <- 0
  # Timing (information fraction) at interim analyses
  timing <- timing_info
  # Efficacy bound spending function
  sfu <- sfu_type
  # Upper bound spending function parameters, if any
  sfupar <- c(0)
  # Lower bound spending function, if used (test.type > 2)
  sfl <- sfLDOF
  # Lower bound spending function parameters, if any
  sflpar <- c(0)
  endpoint <- 'normal' 
  # Variance
  sigma_2 <- 1
  alpha_overall <- 0.05
  delta <- 0
  
  n <- nNormal(
    delta1 = trt_diff,
    delta0 = 0,
    ratio = 1,
    sd = 1,
    sd2 = 1,
    alpha = alpha,
    beta = beta
  )
  n <- ceiling(n/(1-theta))
  
  x <- gsDesign(
    k = k,
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
    delta1 = trt_diff,
    delta0 = 0,
    endpoint = endpoint,
    n.fix = n
  )
  
  # Fixed design sample size
  n <- ((qnorm(1-beta,mean=0,sd=1)+qnorm(1-alpha_overall/2,mean=0,sd=1))^2/
          (4*theta*phi^2*(1-phi)^2*select_diff^2))*
    (sigma_2 + phi*(1-phi)*(select_diff*(2*phi-1) + prefer_diff)^2 +
       2*(theta/(1-theta))*(phi^2 + (1-phi)^2)*sigma_2)
  
  x_inflate <- gsDesign(
    k = k,
    test.type = test.type,
    alpha = alpha,
    beta = beta,
    astar = astar,
    timing = timing,
    sfu = sfu,
    sfupar = sfupar,
    sfl = sfl,
    sflpar = sflpar,
    delta = 0,
    endpoint = endpoint,
    n.fix=1
  )
  
  inflation <- x_inflate$n.I
  n_lst <- ceiling(inflation*n)
  # Maximum group sequential sample size
  N <- n_lst[length(n_lst)]
  
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
      dataset_lst[[i]] <- dataset(q[i]*N,phi,theta,trt_diff,select_diff,prefer_diff,mu)
    }
    dataset_cum <- dataset_lst
    for (i in 2:K){
      dataset_cum[[i]] <- rbind(dataset_cum[[i-1]], dataset_cum[[i]])
    }
    
    # Calculate test statistics
    test_stats <- test_select_diff(dataset_cum,K)
    for (i in 1:K){
      # Compare with boundary statistics
      if(test_stats[i] < x$lower$bound[i]|test_stats[i] > x$upper$bound[i]){
        ret_stopping_size[n_loop] <- N*q_cum[i]
        break
      }
      if (i == K){
        ret_stopping_size[n_loop] <- N
      }
    }
  }
  
  # Empirical Type I error
  type1_error <- error_rate_select(2,x$upper$bound,
                            x$lower$bound,N,phi,theta,trt_diff,0,
                            prefer_diff,mu,q,K,q_cum)
  
  # Empirical power
  emp_power_select <- error_rate_select(2,x$upper$bound,
                                        x$lower$bound,N,phi,theta,trt_diff,select_diff,
                                        prefer_diff,mu,q,K,q_cum)
  
  return(list(ret_stopping_size,ceiling(n),type1_error,emp_power_select))
}

Selection_Difference <- seq(0.2,1,0.1)
Seq_Size_Pocock_select <- c()
Seq_Size_OF_select <- c()
Fixed_Size_select <- c()
sd_lst_pocock <- c()
sd_lst_OBF <- c()
type1_pocock <- c()
type1_OBF <- c()
emp_power_pocock <- c()
emp_power_OBF <- c()



for (i in 1:length(Selection_Difference)){
  # Pocock Stopping Sample Size 
  ret_pocock <- stopping_sample_size_select(sample_size_num_loop,phi,theta,trt_diff,
                                     Selection_Difference[i],
                                     prefer_diff,mu,K,"Pocock",timing_info)
  Seq_Size_Pocock_select[i] <- ceiling(mean(ret_pocock[[1]]))
  Fixed_Size_select[i] <- ret_pocock[[2]]
  sd_lst_pocock[i]<- round(sqrt(var(ret_pocock[[1]])),2)
  type1_pocock[i] <- ret_pocock[[3]]
  emp_power_pocock[i] <- ret_pocock[[4]]
  
  # OBF Stopping Sample Size 
  ret_OBF <- stopping_sample_size_select(sample_size_num_loop,phi,theta,trt_diff,
                                         Selection_Difference[i],
                                         prefer_diff,mu,K,"OF",timing_info)
  Seq_Size_OF_select[i] <- ceiling(mean(ret_OBF[[1]]))
  sd_lst_OBF[i]<- round(sqrt(var(ret_OBF[[1]])),2)
  type1_OBF[i] <- ret_OBF[[3]]
  emp_power_OBF[i] <- ret_OBF[[4]]
}


df_select <- cbind(Selection_Difference,Seq_Size_Pocock_select,Seq_Size_OF_select,Fixed_Size_select,
                   sd_lst_pocock,type1_pocock,emp_power_pocock,sd_lst_OBF,type1_OBF,emp_power_OBF)
df_select <- data.frame(df_select)
colnames(df_select) <- c('Selection Difference', 'Pocock Sample Size', 'OBF Sample Size',
                  'Fixed Sample Size','Pocock SD', 'Pocock Type 1','Pocock Power','OBF SD',
                  'OBF Type 1','OBF Power')
View(df_select)

values <- c(Seq_Size_Pocock_select/Fixed_Size_select,Seq_Size_OF_select/Fixed_Size_select)
group <- c(rep("Pocock/Fixed",length(Selection_Difference)),rep("OBF/Fixed",length(Selection_Difference)))
results <- data.frame(Selection_Difference,group,values)
print(results)

# plot
g_selection <- ggplot(data = results, aes(x=Selection_Difference, y=values, group = group)) +geom_line(aes(color=group))+
  geom_point(aes(color=group))  + labs(title = "\n", x = "Selection Effect",
                                       y = "Ratio")

################ Preference Difference ####################  

# Function for calculating empirical power and empirical Type I error rate
error_rate_prefer <- function(num_loops,upper_bound,lower_bound,N,phi,theta,trt_diff,
                              select_diff,prefer_diff,mu,q,K,q_cum){
  reject_lst <- c()
  for (n_loop in 1:num_loops){
    dataset_lst <- list()
    # Generate Datasets
    for (i in 1:K){
      dataset_lst[[i]] <- dataset(q[i]*N,phi,theta,trt_diff,select_diff,prefer_diff,mu)
    }
    dataset_cum <- dataset_lst
    for (i in 2:K){
      dataset_cum[[i]] <- rbind(dataset_cum[[i-1]], dataset_cum[[i]])
    }
    # Calculate test statistics
    test_stats <- test_prefer_diff(dataset_cum,K)
    for (i in 1:K){
      if(test_stats[i] < lower_bound[i]|test_stats[i] > upper_bound[i]){
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

# Test statistics
test_prefer_diff <- function(dataset_cum,K){
  
  pi_lst <- c()
  for (i in 1:K){
    df <- dataset_cum[[i]]
    # Define m_1, m_2, n_1, and n_2
    m_1 <- length(df$Y[df$Treatment == 1 & df$first_stage == "Option"]) 
    m_2 <- length(df$Y[df$Treatment == 0 & df$first_stage == "Option"]) 
    m = m_1 + m_2  
    n_1 <- length(df$Y[df$Treatment == 1 & df$first_stage == "Random"]) 
    n_2 <- length(df$Y[df$Treatment == 0 & df$first_stage == "Random"]) 
    
    # Empirical Variance
    s11_sq <- ifelse(m_1 > 1,var(df$Y[df$Treatment == 1 & df$first_stage == "Option"]),0 )
    s22_sq <- ifelse(m_2 > 1, var(df$Y[df$Treatment == 0 & df$first_stage == "Option"]),0)
    s1_sq <- ifelse(n_1 >1, var(df$Y[df$Treatment == 1 & df$first_stage == "Random"]),0)
    s2_sq <- ifelse(n_2 >1, var(df$Y[df$Treatment == 0 & df$first_stage == "Random"]),0)
    
    # Define z_1 and z_2
    z_1 <- sum(df$Y[df$Treatment == 1 & df$first_stage == "Option"]) -
      m_1 *mean(df$Y[df$Treatment == 1 & df$first_stage == "Random"])
    z_2 <- sum(df$Y[df$Treatment == 0 & df$first_stage == "Option"]) -
      m_2 *mean(df$Y[df$Treatment == 0 & df$first_stage == "Random"])
    
    # Calculate variance of z_1-z_2
    if (m_1 != 0 & m_2 != 0){
      var_z1 <- m_1*s11_sq+(m_1*s1_sq/n_1)*(1+m_1*(m-1)/m)+ (m_1*m_2/m)*
        (mean(df$Y[df$Treatment == 1 & df$first_stage == "Option"])-
           mean(df$Y[df$Treatment == 1 & df$first_stage == "Random"]))^2
      
      var_z2 <- m_2*s22_sq+(m_2*s2_sq/n_2)*(1+m_2*(m-1)/m)+ (m_1*m_2/m)*
        (mean(df$Y[df$Treatment == 0 & df$first_stage == "Option"])-
           mean(df$Y[df$Treatment == 0 & df$first_stage == "Random"]))^2
      
      cov_z1_z2 <- -(m_1*m_2/m)*(mean(df$Y[df$Treatment == 1 & df$first_stage == "Option"])-
                                   mean(df$Y[df$Treatment == 1 & df$first_stage == "Random"]))*
        (mean(df$Y[df$Treatment == 0 & df$first_stage == "Option"])-
           mean(df$Y[df$Treatment == 0 & df$first_stage == "Random"]))
    }
    
    if (m_1 == 0){
      var_z1 <- 0
      var_z2 <- m_2*s22_sq+(m_2*s2_sq/n_2)*(1+m_2*(m-1)/m)
      cov_z1_z2 <- 0
    }
    if (m_2 == 0){
      var_z1 <- m_1*s11_sq+(m_1*s1_sq/n_1)*(1+m_1*(m-1)/m)
      var_z2 <- 0
      cov_z1_z2 <- 0
    }
    var_z1z2 <- var_z1 + var_z2 + cov_z1_z2
    pi_lst[i]<- (z_1 + z_2)/sqrt(var_z1z2)
  }
  return(pi_lst)
}

# Average GSD sample size
stopping_sample_size_prefer <- function(num_loops,phi,theta,trt_diff,select_diff,
                                        prefer_diff,mu,K,sfu_type,timing_info){
  
  # Number of analyses
  k <- K
  # See gsDesign() help for description of boundary types
  test.type <- 2
  # 1-sided Type I error
  alpha <- 0.025
  # Type 2 error (1 - targeted power)
  beta <- 0.1
  astar <- 0
  # Timing (information fraction) at interim analyses
  timing <- timing_info
  # Efficacy bound spending function
  sfu <- sfu_type
  # Upper bound spending function parameters, if any
  sfupar <- c(0)
  # Lower bound spending function, if used (test.type > 2)
  sfl <- sfLDOF
  # Lower bound spending function parameters, if any
  sflpar <- c(0)
  delta <- 0
  endpoint <- 'normal'  
  
  n <- nNormal(
    delta1 = trt_diff,
    delta0 = 0,
    ratio = 1,
    sd = 1,
    sd2 = 1,
    alpha = alpha,
    beta = beta
  )
  n <- ceiling(n/(1-theta))
  
  x <- gsDesign(
    k = k,
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
    delta1 = trt_diff,
    delta0 = 0,
    endpoint = endpoint,
    n.fix = n
  )
  sigma_2 <- 1
  alpha_overall <- 0.05
  
  # Fixed design sample size
  n <- ((qnorm(1-beta,mean=0,sd=1)+qnorm(1-alpha_overall/2,mean=0,sd=1))^2/
          (4*theta*phi^2*(1-phi)^2*prefer_diff^2))*
    (sigma_2 + phi*(1-phi)*(prefer_diff*(2*phi-1) + select_diff)^2 +
       2*(theta/(1-theta))*(phi^2 + (1-phi)^2)*sigma_2)
  
  x_inflate <- gsDesign(
    k = k,
    test.type = test.type,
    alpha = alpha,
    beta = beta,
    astar = astar,
    timing = timing,
    sfu = sfu,
    sfupar = sfupar,
    sfl = sfl,
    sflpar = sflpar,
    delta = 0,
    endpoint = endpoint,
    n.fix=1
  )
  inflation <- x_inflate$n.I
  n_lst <- ceiling(inflation*n)
  # Maximum group sequential sample size
  N <- n_lst[length(n_lst)]
  
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
      dataset_lst[[i]] <- dataset(q[i]*N,phi,theta,trt_diff,select_diff,prefer_diff,mu)
    }
    dataset_cum <- dataset_lst
    for (i in 2:K){
      dataset_cum[[i]] <- rbind(dataset_cum[[i-1]], dataset_cum[[i]])
    }
    
    # Calculate test statistics
    test_stats <- test_prefer_diff(dataset_cum,K)
    for (i in 1:K){
      if(test_stats[i] < x$lower$bound[i]|test_stats[i] > x$upper$bound[i]){
        ret_stopping_size[n_loop] <- N*q_cum[i]
        break
      }
      if (i == K){
        ret_stopping_size[n_loop] <- N
      }
    }
  }
  
  # Empirical Type I error
  type1_error <- error_rate_prefer(2,x$upper$bound,
                            x$lower$bound,N,phi,theta,trt_diff,select_diff,
                            0,mu,q,K,q_cum)

  # Empirical power
  emp_power_prefer <- error_rate_prefer(2,x$upper$bound,
                                   x$lower$bound,N,phi,theta,trt_diff,select_diff,
                                   prefer_diff,mu,q,K,q_cum)
  
  return(list(ret_stopping_size,ceiling(n),type1_error,emp_power_prefer))
}

Preference_Difference <- seq(0.2,1,0.1)
Seq_Size_Pocock_prefer <- c()
Seq_Size_OF_prefer <- c()
Fixed_Size_prefer <- c()
sd_lst_pocock <- c()
sd_lst_OBF <- c()
type1_pocock <- c()
type1_OBF <- c()
emp_power_pocock <- c()
emp_power_OBF <- c()

for (i in 1:length(Preference_Difference)){
  # Pocock Stopping Sample Size 
  ret_pocock <- stopping_sample_size_prefer(sample_size_num_loop,phi,theta,trt_diff,
                                     select_diff,
                                     Preference_Difference[i],mu,K,"Pocock",timing_info)
  
  Seq_Size_Pocock_prefer[i] <- ceiling(mean(ret_pocock[[1]]))
  Fixed_Size_prefer[i] <- ret_pocock[[2]]
  sd_lst_pocock[i]<- round(sqrt(var(ret_pocock[[1]])),2)
  type1_pocock[i] <- ret_pocock[[3]]
  emp_power_pocock[i] <- ret_pocock[[4]]
  
  # OBF Stopping Sample Size
  ret_OBF <- stopping_sample_size_prefer(sample_size_num_loop,phi,theta,trt_diff,
                                         select_diff,
                                         Preference_Difference[i],mu,K,"OF",timing_info)
  
  Seq_Size_OF_prefer[i] <- ceiling(mean(ret_OBF[[1]]))
  sd_lst_OBF[i]<- round(sqrt(var(ret_OBF[[1]])),2)
  type1_OBF[i] <- ret_OBF[[3]]
  emp_power_OBF[i] <- ret_OBF[[4]]
}

df_prefer <- cbind(Preference_Difference,Seq_Size_Pocock_prefer,Seq_Size_OF_prefer,Fixed_Size_prefer,
                   sd_lst_pocock,type1_pocock,emp_power_pocock, sd_lst_OBF,type1_OBF,emp_power_OBF)
df_prefer <- data.frame(df_prefer)
colnames(df_prefer) <- c('Preference Difference', 'Pocock Sample Size', 'OBF Sample Size',
                  'Fixed Sample Size','Pocock SD', 'Pocock Type 1','Pocock Power','OBF SD','OBF Type 1','OBF Power')
View(df_prefer)

values <- c(Seq_Size_Pocock_prefer/Fixed_Size_prefer,Seq_Size_OF_prefer/Fixed_Size_prefer)
group <- c(rep("Pocock/Fixed",length(Preference_Difference)),rep("OBF/Fixed",length(Preference_Difference)))
results <- data.frame(Preference_Difference,group,values)
print(results)

# plot
g_preference <- ggplot(data = results, aes(x=Preference_Difference, y=values, group = group)) +geom_line(aes(color=group))+
  geom_point(aes(color=group))  + labs(title = "\n", x = "Preference Effect",
                                     y = "Ratio")


# Combine plots
g_combine <- ggarrange(g_treatment,g_selection,g_preference,ncol=1,heights=2)
annotate_figure(g_combine, top = text_grob("Change of Average Sample Size Ratio by Test Effect Size", 
                                           color = "black", face = "bold", size = 12))


