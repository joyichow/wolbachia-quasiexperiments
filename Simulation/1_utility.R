# assumptions:
# 1. homogeneous mixing mosquitoes population
# 2. no dispersal of mosquitoes
# 3. Wolbachia does not influence the immature stage - eggs, larvae and pupae are grouped into one aquatic stage
# 4. Enviroment conditions remain constant - motality rates constant
# 4. W-infected male mosquitoes are equally successful in finding and mating with female
################################################################################
# simulate evenly release of mosquitoes
# match to mature wolbachia/non-wobachia mosquito death rate
generate_intervention <- function(c_intervention, start, end, population, mortality) { 
  d_intervention <- numeric (iteration)
  # Calculate the mean increment per iteration
  mean_increment <- population / (end - start) 
  # Initialize the first value
  c_intervention[start] <- c_intervention[start] + mean_increment
  # generate cumulative sequence
  for (i in (start + 1):iteration) {
    if (i <= end) {
      c_intervention[i] <- c_intervention[i-1] * (1 - μw) + mean_increment
    } 
    else {
      c_intervention[i] <- c_intervention[i-1] * (1 - μw)
    }
  }
  # generate derivative sequence
  for (i in 1:(iteration+1)) {
    if (i == (iteration+1)) {
      d_intervention[i] <- c_intervention[i]
    }
    else{
      d_intervention[i] <-  c_intervention[i+1] - c_intervention[i]
    }
  }
  # combine result
  intervention <- data.frame(cumulative=c_intervention, derivative=d_intervention)
  return(intervention)
}



# add poison noise return integer values
add_poisson_noise <- function(sequence) {
  noisy_sequence <- numeric(length(sequence))
  
  for (i in 1:length(sequence)) {
    noisy_sequence[i] <- rpois(1, sequence[i]) # generate random number given value as mean per interval
  }
  
  return(noisy_sequence)
}


run_iit <- function(iteration, mosquitoes){
  for (j in seq_len(iteration)) {
    ########################################################
    # calculating constant for each iteration
    Nh <- Sh[j] + Eh[j] + Ih[j] + Rh[j] #👨👩: entire population
    K <- 100000                         #🦟: carrying capacity
    Fn <- Sn[j] + En[j] + In[j]         #🦟:㊛ non-wolbachia 
    Fw <- 0                             #🦟:㊛ wolbachia
    Mn <- Sn[j] + En[j] + In[j]         #🦟:㊚ non-wolbachia 
    Mw <- intervention$cumulative[j]    #🦟:㊚ wolbachia
    # Sh are infected when bitten by In and Iw at rate of λn and λw
    λnh <- (Bn*Tn*In[j])/Nh
    λhn <- (Bn*Tn*Ih[j])/Nh
    ########################################################
    # Human total Nh population in SEIR compartments
    
    dSh <- B*Nh - λnh*Sh[j] - μh*Sh[j] + Ω*Rh[j]
    # susceptible <- born - infected(n)- die + return
    
    dEh <- λnh*Sh[j]- γh*Eh[j] - μh*Eh[j]
    # exposed <- infected(n)- progress - die
    
    dIh <- γh*Eh[j] - σ*Ih[j] - μh*Ih[j]
    # infectious <- progress - recover - die
    
    dRh <- σ*Ih[j] - μh*Rh[j] - Ω*Rh[j]
    # recovered <- recover - die - return
    
    ########################################################
    # mosquito K population, F denoted for female mosquitoes, 
    # Non-wolbachia mosquitoes An population in ASEI compartments
    
    dAn <- ρn*((K-An[j])/K)*(Fn*Mn/(Fn+Mn+Mw)) - τn*An[j] - μna*An[j]
    # Aquatic(n) <- produced - mature - die(n_aquatic)
    #               produced assumption: F:M = 1:1 so, 
    #               produced = rep_ratio * capacity * (Fn*Mn/(Fn+Fw+Mn+Mw))
    
    dSn <- (1/2)*τn*An[j] - λhn*Sn[j] - μn*Sn[j]
    # Suspected(n) <- female_mature  - bite - die(n_mature)
    
    dEn <-  λhn*Sn[j] - γn*En[j] - μn*En[j]
    # Exposed(n) <- bite - progress - die(n_mature)
    
    dIn <- γn*En[j] - μn*In[j]
    # Infected(n) <- progress - die(n_mature)
    
    #########################################################
    dfSh[j] <- dSh
    dfEh[j] <- dEh
    dfIh[j] <- dIh
    dfRh[j] <- dRh

    dfAn[j] <- dAn
    dfSn[j] <- dSn
    dfEn[j] <- dEn
    dfIn[j] <- dIn

    
    Sh[j+1] <- Sh[j] + dSh
    Eh[j+1] <- Eh[j] + dEh
    Ih[j+1] <- Ih[j] + dIh
    Rh[j+1] <- Rh[j] + dRh

    An[j+1] <- An[j] + dAn
    Sn[j+1] <- Sn[j] + dSn
    En[j+1] <- En[j] + dEn
    In[j+1] <- In[j] + dIn

    
  }

  result <- data.frame(Sh, Eh, Ih, Rh, An, Sn, En, In, dfSh, dfEh, dfIh, dfRh, dfAn, dfSn, dfEn, dfIn)
  return(result)
}



run_introgression <- function(iteration, mosquitoes){
  for (j in seq_len(iteration)) {
    ########################################################
    # calculating constant for each iteration
    Nh <- Sh[j] + Eh[j] + Ih[j] + Rh[j]
    K  <- 100000 #🦟: carrying capacity
    Fn <- Sn[j] + En[j] + In[j]
    Fw <- Sw[j] + Ew[j] + Iw[j] 
    Mn <- Sn[j] + En[j] + In[j]
    Mw <- Sw[j] + Ew[j] + Iw[j]
    # Sh are infected when bitten by In and Iw at rate of λn and λw
    λnh <- (Bn*Tn*In[j])/Nh
    λwh <- (Bw*Thw*Iw[j])/Nh
    λhn <- (Bn*Tn*Ih[j])/Nh
    λhw <- (Bw*Thw*Ih[j])/Nh
    ########################################################
    # Human total Nh population in SEIR compartments
    
    dSh <- B*Nh - λnh*Sh[j] - λwh*Sh[j] - μh*Sh[j] + Ω*Rh[j]
  # susceptible <- born - infected(n) - infected(w) - die + return
  
    dEh <- λnh*Sh[j] + λwh*Sh[j] - γh*Eh[j] - μh*Eh[j]
  # exposed <- infected(n) - infected(w) - progress - die
    
    dIh <- γh*Eh[j] - σ*Ih[j] - μh*Ih[j]
  # infectious <- progress - recover - die
    
    dRh <- σ*Ih[j] - μh*Rh[j] - Ω*Rh[j]
  # recovered <- recover - die - return
    
    ########################################################
    # mosquito K population, F denoted for female mosquitoes, 
    # Non-wolbachia mosquitoes An population in ASEI compartments
    
    dAn <- ρn*((K-An[j]-Aw[j])/K)*((Fn*Mn)/(Fn+Mn+Fw+Mw)) - τn*An[j] - μna*An[j]
  # Aquatic(n) <- produced - mature - die(n_aquatic)
  #               assumption: F:M = 1:1 so, 
  #               produced = reproductive_ratio * capacity * (Fn*Mn/(Fn+Fw+Mn+Mw))
    
    dSn <- (1/2)*τn*An[j] + (1-ɑ)*τw*(1/2)*Aw[j] - λhn*Sn[j] - μn*Sn[j]
  # Suspected(n) <- female_mature + maternal_transmission_leakage - bite - die(n_mature)
    
    dEn <-  λhn*Sn[j] - γn*En[j] - μn*En[j]
  # Exposed(n) <- bite - progress - die(n_mature)
    
    dIn <- γn*En[j] - μn*In[j]
  # Infected(n) <- progress - die(n_mature)
    
    #########################################################
    # Wolbachia mosquitoes Aw population in ASEI compartments
    dAw <- ρw*((K-An[j]-Aw[j])/K)*((Fw*Mw+Fw*Mn)/(Fw+Mw+Fn+Mn))- τw*Aw[j] - μwa*Aw[j]

    # Aquatic(w) <- produced - mature - die(w_aquatic)
    #             produced = reproductive_ratio * capacity * (Fw*Mw/Fw+Mw)
    
    dSw <- (1/2)*ɑ*τw*Aw[j] + (1/2)*intervention$derivative[j] - λhw*Sw[j] - μw*Sw[j]
  # Suspected(w) <- female_mature(maternal_transmission) + released - bite - die(w_mature)
    
    dEw <-  λhw*Sw[j] - γw*Ew[j] - μw*Ew[j]
  # Exposed(w) <- bite - progress - die(w_mature)
    
    dIw <- γw*Ew[j] - μw*Iw[j]
  # Infected(w) <- progress - die(w_mature)
    
    #########################################################
    Sh[j+1] <- Sh[j] + dSh
    Eh[j+1] <- Eh[j] + dEh
    Ih[j+1] <- Ih[j] + dIh
    Rh[j+1] <- Rh[j] + dRh

    An[j+1] <- An[j] + dAn
    Sn[j+1] <- Sn[j] + dSn
    En[j+1] <- En[j] + dEn
    In[j+1] <- In[j] + dIn

    Aw[j+1] <- Aw[j] + dAw
    Sw[j+1] <- Sw[j] + dSw
    Ew[j+1] <- Ew[j] + dEw
    Iw[j+1] <- Iw[j] + dIw
    

  }
  
  result <- data.frame(Sh, Eh, Ih, Rh, An, Sn, En, In, Aw, Sw, Ew, Iw)
  return(result)
}


