# Description  ------------------------------------------------------------
#### ### ### ## #### ### ### ## #### ### ### ## 
# Genetic drift simulation
# Created by Marc-Olivier Beausoleil
# 2022-01-07
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






# Generate coalescence ----------------------------------------------------
# Maximum population size 
pop = 5 # The maximum at the moment is 13 
n.genes = c(pop*2)
# Total number of gametes in the population 
n.sperm = 20
n.eggs = 20
n.gametes = c(pop*(n.sperm+n.eggs))

# Get the fake genes names
x <- LETTERS[1:n.genes]
# Make data frame with genes in it 
coalescence <- data.frame(gene = matrix(c(x),ncol = 1, nrow = n.genes))
# First generation of genes 
coalescence$gen1 = 1
# Number of generations 
gen = 100

# Loop to track the within population allele frequency change 
for (i in 1:gen) {
  # Get the last column index 
  last.col.index = ncol(coalescence)
  # Find the genes present from the last generation
  genes.present = coalescence$gene[coalescence[,last.col.index]>0]
  # Get the probability of having a gene reproduce based on its frequency 
  probs.genes = coalescence[,last.col.index][coalescence[,last.col.index]>0]/n.genes
  
  # binomial function to generate the new allele frequency (0 = q, 1 = p)
  # allele.fq = rbinom(n = n.gametes, size = 1, prob = prob.p)
  
  # Number of genes present 
  max.nb.gene = length(genes.present)
  
  # Get the genes reproduce
  drifted.genes = sample(x = genes.present, size = n.gametes, replace = T, prob = probs.genes)
  # gametes.to.combine = drifted.genes[as.logical(allele.fq)]
  
  # Randomly sample the population (this is the drift, a random sample of the population)
  all.drift = sample(x = drifted.genes, size = n.genes, replace = F)
  
  # Combine the data 
  temp.df = as.data.frame(table(all.drift))
  coalescence = merge(x = coalescence,y = temp.df,
                      by.x = "gene",by.y = "all.drift", 
                      all.x = TRUE)
  # Replace NAs with 0
  coalescence[is.na(coalescence)] <- 0
  # Make last column name
  lastname = paste0("gen",i+1)
  names(coalescence)[names(coalescence) == 'Freq'] <- lastname
  
  # Check if the genes are at fixation 
  test = length(which(unlist(coalescence[names(coalescence) == lastname]) %in% 10))
  if(test > 0) {break}
}
coalescence


# Gene tree ---------------------------------------------------------------
# Futuyma p. 171 figure 7.5
# Make a plot with the genes in all generations 
gene.seq = 1:n.genes
y = 1:(ncol(coalescence)-1)
# Make empty plot 
plot(NA,type = "n",
     ylab = "Generations",
     xlab = "Genes",asp = 1,
     ylim = range(y),xlim = range(gene.seq))

# Make the data frame with all genes in each generation
newdat = NULL
for (j in 2:ncol(coalescence)) {
  tmp = rep(coalescence$gene,coalescence[,j])
  df.tmp = data.frame(tmp,id = 1:n.genes)
  newdat = rbind(newdat,tmp)
  text(gene.seq,y = j-1,tmp,col = as.numeric(factor(tmp,levels = LETTERS[1:n.genes]))+1)
  newdat[nrow(newdat)-1,]
  newdat[nrow(newdat),]
}
# Add the segments based on the relationship between the genes 
for (z in 1:n.genes) {
  # Gene to be tracked 
  check.gene = LETTERS[z]
  for (k in 1:(ncol(coalescence)-2)) {
    # check the position of each gene in 2 consecutive generations 
    gen1 = (1:n.genes)[newdat[k,] ==check.gene]
    gen2 = (1:n.genes)[newdat[k+1,] ==check.gene]
    # If nothing in the next generation, skip the iteration 
    if(length(gen2)==0) {break}
    # If there is less elements in gen2, subset gen1 to the length of gene2
    if (length(gen2)<length(gen1)) {gen1 = gen1[1:length(gen2)]}
    # make a temporary data to get the position of gene in the 2 generations 
    ddd = rbind(gen1,gen2)
    # Reorder the genes 
    ddd[1,] = sort(ddd[1,])
    
    text(x = unique(ddd[1,][duplicated(ddd[1,])]),y = k,labels = check.gene,font =2)
    # Draw the segment from generation 1 to generation 2 
    segments(x0 = ddd[1,],
             y0 = k,
             x1 = ddd[2,],
             y1 = k+1,
             col = z+1) # colour based on the iteration number (the gene)
  } # End of k
} # End of z
