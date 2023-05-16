# Description -------------------------------------------------------------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created March 18, 2023

# Why:
  # selection changes allele frequencies (See how selection changes allele frq)

# Requires 

# NOTES: 
# Source: 
  # See p. 194 in Evolution: making sense of life 3rd edition
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

select.change.allele.fq <- function(generations = 10, 
                                    p = .5,
                                    w.geno = c(1.6, 1.2, 1)
                                    ) {
  pop = NULL
  delta_p.all = NULL
  avg.excess.fit.all = NULL
  i_plus = 0
  for (i in 1:generations) {
    # cat("Doing pop",i,fill = TRUE)
    gen1 = p^2
    gen2 = 2*p*(1-p)
    gen3 = (1-p)^2
    gen = c(gen1,gen2,gen3)
    w = c(w.geno[1],
          w.geno[2],
          w.geno[3])
    
    # Average fitness in population 
    w_bar = sum(gen*w)
    # Get new genotype frequencies 
    freq_t1 = gen*w/w_bar
    # Table of genotype frequencies before selection, fitness, genotype frequencies after selection
    # rbind(gen,w,freq_t1)
    
    # Get the allele fq after selection 
    p_t1 = freq_t1[1] + .5 * freq_t1[2]
    q_t1 = freq_t1[3] + .5 * freq_t1[2]
    
    # If 1st generation, keep the initial parameters
    if (i == 1) {
      tmp = data.frame(p = c(p, p_t1),q = c((1-p), q_t1), 
                       row.names = c(paste("t",i-1,sep = "_"), 
                                     paste("t",i,sep = "_")))  
    } else { # Add only the new pop 
      tmp = data.frame(p = p_t1, q = q_t1, 
                       row.names = c(paste("t",i,sep = "_")))
    }
    
    # Average excess of fitness for allele 
    avg.excess.fit = ((p*(w.geno[1]-w_bar))+((1-p)*(w.geno[2]-w_bar)))
    delta_p = p/w_bar*avg.excess.fit
    avg.excess.fit.all = c(avg.excess.fit.all, avg.excess.fit)
    delta_p.all = c(delta_p.all, delta_p)
    # Combine the data 
    pop = rbind(pop, tmp)
    # Reset the allele frequency 
    p = p_t1
  }
  # Make generations as column 
  pop$time = 0:c(nrow(pop)-1)
  # Return the data 
  return(list(pop = pop, delta_p = delta_p.all, avg.excess.fit = avg.excess.fit.all))
}

par(mfrow = c(1,2))
max.gen = 100
w.geno = c(1.9,1,1)
plot(NA, xlim = c(0,max.gen), ylim = c(0,1), main = "Allele fq change \n through generations", ylab = "Allele frequency p", type = "l")
for (k in seq(0.01, .9, by = .01)) {
  sel.dat = select.change.allele.fq(generations = max.gen, p = k, w.geno = w.geno)
  with(sel.dat$pop, points(time, p, type = "l", col = grey(level = k)))
}
legend("topright", legend = paste("Fitness of genotype",1:3,"—",w.geno), cex = .6)

sel.dat = select.change.allele.fq(generations = 10, p = .001, w.geno = w.geno)
sel.dat$avg.excess.fit
plot(sel.dat$avg.excess.fit, main = 'Average excess fitness', type = "l", ylim = c(-1,1))
