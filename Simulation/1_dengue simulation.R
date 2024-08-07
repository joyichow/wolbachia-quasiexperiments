library(ggplot2)
library(openxlsx)
source("/your/local/filepath/dengue simulator/utility.R")
###############################################################################
# initialization - create record frame with iteration times
{
  iteration <- 1000
  c_intervention <- numeric(iteration + 1)

  Sh <- numeric(iteration + 1)
  Eh<- numeric(iteration + 1)
  Ih<- numeric(iteration + 1)
  Rh<- numeric(iteration + 1)

  An<- numeric(iteration + 1)
  Sn<- numeric(iteration + 1)
  En<- numeric(iteration + 1)
  In<- numeric(iteration + 1)

  Aw<- numeric(iteration + 1)
  Sw<- numeric(iteration + 1)
  Ew<- numeric(iteration + 1)
  Iw<- numeric(iteration + 1)
  
  dfSh <- numeric(iteration + 1)
  dfEh<- numeric(iteration + 1)
  dfIh<- numeric(iteration + 1)
  dfRh<- numeric(iteration + 1)
  
  dfAn<- numeric(iteration + 1)
  dfSn<- numeric(iteration + 1)
  dfEn<- numeric(iteration + 1)
  dfIn<- numeric(iteration + 1)
  
  dfAw<- numeric(iteration + 1)
  dfSw<- numeric(iteration + 1)
  dfEw<- numeric(iteration + 1)
  dfIw<- numeric(iteration + 1)
}

# set parameter value
{
  # biological rate
  Ω  <- 1/(365*2)      # waning immunity
  σ  <- 1/10           # human recover rate
  τn <- 0.12           # non-wolbachia juvenile mature rate
  τw  <- 0.12          # wolbachia juvenile mature rate
  ɑ  <- 0.9            # maternal transmission rate
  Bn  <- 0.7           # biting rate of non-wolbachia mosquito
  Bw  <- 0.5 * Bn     # biting rate of wolbachia mosquito
  # dengue rate
  γh  <- 1/6           # human dengue progression rate
  γn  <- 1/10          # non-wolbachia dengue progression rate
  γw  <- 1/10          # wolbachia dengue progression rate
  Tn  <- 0.2614        # transmission probability 
  Thw  <- 0.5 * Tn     # transmission probability from infectious wolbachia mosquito to human
  # mortality rate
  μh  <- 1/(83.5*365)  # human
  μna  <- 0.06         # non-wolbachia-aquatic
  μn  <- 0.02709       # non-wolbachia-mature
  μwa  <- 0.06         # wolbachia-aquatic
  μw <- 1.1 * μn       # wolbachia-mature
  # birth rate
  B <- 1/(83.5*365)    # human birth rate
  #ρn  <- 5.7           # non-wolbachia
  #ρw  <- 3.24          # wolbachia
  ρn  <- 1.25           # non-wolbachia
  ρw  <- 0.95*ρn          # wolbachia
}
###############################################################################
# simulaten male-wolbachia mosquitoes release
intervention<- generate_intervention(c_intervention, 50, 150, 100000, μw)
plot.ts(intervention$cumulative)
plot.ts(intervention$derivative)

###############################################################################
# Incompatible Insect Technology Strategy (IIT) release male-wolbachia
# male-wolbachia mosquitoes mate with wildtype female mosquitoes but produce
# no offspring due to Cytoplasmic Incompatibility (CI)

# set up initials for population
{
  Sh[1]<- 9990
  Eh[1]<- 0
  Ih[1]<- 10
  Rh[1]<- 0
  
  An[1]<- 180
  Sn[1]<- 0
  En[1]<- 0
  In[1]<- 20
  
  Aw[1]<- 0
  Sw[1]<- 0
  Ew[1]<- 0
  Iw[1]<- 0
}

# run compartments model with population initials and released mosquito track
result <- run_iit(iteration, mosquitoes)
result$cumulative_intervention <- intervention$cumulative
write.xlsx(result, file = "result.xlsx")

# get poissoned report
poissioned_Ih<-add_poisson_noise(result$Ih)
plot.ts(poissioned_Ih)

# plot human compartments
ggplot(data = result, aes(x = seq_len(iteration + 1))) +
  geom_line(aes(y = Sh, color = "Susceptible")) +
  geom_line(aes(y = Eh, color = "Exposed")) +
  geom_line(aes(y = Ih, color = "Infected")) +
  geom_line(aes(y = Rh, color = "Recovered")) +
  labs(title = "SEIR Model Simulation",
       x = "Time",
       y = "Population",
       color = "Legend") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = c("Susceptible" = "blue", 
                                "Exposed" = "orange",
                                "Infected" = "red", 
                                "Recovered" = "green")) +
  theme_minimal()

# plot mosquito compartments
ggplot(data = result, aes(x = seq_len(iteration + 1))) +
  geom_line(aes(y = An, color = "Aquatic")) +
  geom_line(aes(y = Sn, color = "Suspectible")) +
  geom_line(aes(y = En, color = "Exposed")) +
  geom_line(aes(y = In, color = "Infectious")) +
  geom_line(aes(y = cumulative_intervention, color = "Intervention"))+
  labs(title = "IIT mosquito simulation",
       x = "Time",
       y = "Population",
       color = "Legend") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = c("Aquatic" = "green", 
                                "Suspectible" = "blue",
                                "Exposed" = "orange", 
                                "Infectious" = "red",
                                "Intervention" = "black")) +
  theme_minimal()




################################################################################
# Introgression approach release female-wolbachia mosquitoes 
# female-wolbachia mosquitoes pass wolbachia to offspring at maternal transmission rate

# set up initials for population
{
  Sh[1]<- 9900
  Eh[1]<- 0
  Ih[1]<- 100
  Rh[1]<- 0
  
  An[1]<- 1800
  Sn[1]<- 0
  En[1]<- 0
  In[1]<- 200
  
  Aw[1]<- 0
  Sw[1]<- 0
  Ew[1]<- 0
  Iw[1]<- 0
}


# run simulation
result <- run_introgression(iteration, mosquitoes)
result$cumulative_intervention <- intervention$cumulative

# plot poissoned report
poissioned_Ih<-add_poisson_noise(result$Ih)
plot.ts(poissioned_Ih)

# plot human compartments
ggplot(data = result, aes(x = seq_len(iteration + 1))) +
  geom_line(aes(y = Sh, color = "Susceptible")) +
  geom_line(aes(y = Eh, color = "Exposed")) +
  geom_line(aes(y = Ih, color = "Infected")) +
  geom_line(aes(y = Rh, color = "Recovered")) +
  labs(title = "Introgression: Human",
       x = "Time/day",
       y = "Accumulative Population",
       color = "Legend") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = c("Susceptible" = "blue", 
                                "Exposed" = "orange",
                                "Infected" = "red", 
                                "Recovered" = "green")) +
  theme_minimal()

# plot non-Wolbachia mosquito compartments
ggplot(data = result, aes(x = seq_len(iteration + 1))) +
  geom_line(aes(y = An, color = "Aquatic")) +
  geom_line(aes(y = Sn, color = "Suspectible")) +
  geom_line(aes(y = En, color = "Exposed")) +
  geom_line(aes(y = In, color = "Infectious")) +
  geom_line(aes(y = cumulative_intervention, color = "Intervention"))+
  labs(title = "Introgression: Non-Wolbachia mosquito",
       x = "Time/day",
       y = "Accumulative Population",
       color = "Legend") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = c("Aquatic" = "green", 
                                "Suspectible" = "blue",
                                "Exposed" = "orange", 
                                "Infectious" = "red",
                                "Intervention" = "black")) +
  theme_minimal()

# plot Wolbachia mosquito compartments
ggplot(data = result, aes(x = seq_len(iteration + 1))) +
  geom_line(aes(y = Aw, color = "Aquatic")) +
  geom_line(aes(y = Sw, color = "Suspectible")) +
  geom_line(aes(y = Ew, color = "Exposed")) +
  geom_line(aes(y = Iw, color = "Infectious")) +
  geom_line(aes(y = cumulative_intervention, color = "Intervention"))+
  labs(title = "Introgression: Wolbachia mosquito",
       x = "Time/day",
       y = "Accumulative Population",
       color = "Legend") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = c("Aquatic" = "green", 
                                "Suspectible" = "blue",
                                "Exposed" = "orange", 
                                "Infectious" = "red",
                                "Intervention" = "black")) +
  theme_minimal()


