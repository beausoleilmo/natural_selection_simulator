# Description  ------------------------------------------------------------
#### ### ### ## #### ### ### ## #### ### ### ## 
# Genetic drift simulation
# Created by Marc-Olivier Beausoleil
# 2022-01-12
# Why: 
# Requires:
# NOTES: 
# Drift is (from Futuyma)
# - unbiased
# - random fluctuations in allele frequency are larger in smaller populations
# - drift causes genetic variation to be lost
# - drift causes populations that are initially identical to become different
# - an allele can become fixed without the benefit of natural selection
# Reference : 
# Futuyma p. 167, figure 7.2
#### ### ### ## #### ### ### ## #### ### ### ## 

# graphing parameters -----------------------------------------------------
par(mfrow = c(2,2))

# Random seed -------------------------------------------------------------
# set.seed(1245)

# Simulation parameters ---------------------------------------------------
# Number of gametes to chose from 
n.sperm = 2
n.eggs = 2
# Number of generations (x axis)
gen = 500
# Number of replicate populations 
popu = 5

# variance = p*(1-p)/(2*N)
# Variation is smaller when the population size is bigger 

# Loops -------------------------------------------------------------------
# Loops for all population replicates, tracking allele frequency change over the generations 
# number of individual per population 
n.id.pop = 5*10^seq(0,3, by=1)
# Loop that will change the maximum number of individual per population 
for (l in n.id.pop) {
  # Initial allele frequency 
  p.init = .5 
  # Maximum population size 
  max.pop = l
  # Total number of gametes in the population 
  n.gametes = c(max.pop*(n.sperm+n.eggs))
  # Make an empty object to record the population information 
  all.pops = NULL
  # Loop to track the all population allele frequency change 
  for (j in 1:popu) {
    all.fq.change = .5
    # Loop to track the within population allele frequency change 
    for (i in 1:gen) {
      # If the first iteration, make the probability equal the initial allele frequency 
      if (i == 1) {prob.p = p.init} else {prob.p = prop.all[2]}
      # binomial function to generate the new allele frequency (0 = q, 1 = p)
      allele.fq = rbinom(n = n.gametes, size = 1, prob = prob.p)
      # Randomly sample the population (this is the drift, a random sample of the population)
      all.drift = sample(x = allele.fq, size = max.pop, replace = F)
      # Get the proportion of the alleles in the new population 
      prop.all = prop.table(table(all.drift))
      # Record the p allele only 
      all.fq.change = c(all.fq.change, prop.all[2])
      # If there is an allele that goes to fixation, it'll print NA. In this case, break the for loop and go to the next iteration
      if(is.na(prop.all[2])) {break}
    } # End i
    # Record all population information
    one.pop = data.frame(p.fq = as.numeric(all.fq.change), pop = j)
    all.pops = rbind(all.pops,one.pop)
  } # End j
  
  
  # Remove all NAs ----------------------------------------------------------
  all.pops = na.omit(all.pops)
  
  
  # Plot --------------------------------------------------------------------
  # Make the empty plot 
  plot(all.pops$p.fq~c(1:nrow(all.pops)), 
       col = as.factor(all.pops$pop),
       main = paste("Population N =",max.pop, "\nStarting p =",p.init),
       ylab = "Allele frquency p",
       xlab = "Generations",
       ylim = c(0,1),
       xlim = c(1,gen),
       type = "n")
  
  # Add the lines per population and colour them 
  for (k in 1:popu) {
    pttmp = all.pops[all.pops$pop==k,]
    points(c(1:nrow(pttmp)), pttmp$p.fq, 
           type = "l",
           col = k)
  } # End k
} # End l