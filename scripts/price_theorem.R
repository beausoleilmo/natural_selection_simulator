# Price equations

# ∆E_{i element of I} = Cov_{i element of I}(w_iz_i)+E_{i element of I}(w_i∆z_i)

price.theorem <- function(iter = 100, 
                          n = 100, 
                          n.off = 100,
                          beta1 = 0,
                          mean.p = 10, 
                          mean.o = 0, # 
                          max.off = 10, 
                          sd.p = 1, sd.o = 1,
                          plothist = TRUE, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Assign relative abundance to each parent (must be between 0 and 1 )
  
  rec.exp.price = NULL
  rec.select = NULL
  rec.transm = NULL
  for(i in 1:iter){
    pop.df = data.frame(parent.id = 1:n,
                        rel.abuns = rep(1/n, n) ,# round(sum(rel.abuns)) == 1
                        z.par = rnorm(n = n, mean = mean.p, sd = sd.p))
    
    # Nb offspring given by binomial 
    pop.df$z.par.scale = scale(pop.df$z.par)
    z = beta1*pop.df$z.par.scale
    pr = 1/(1+exp(-z)) # Transform to get the LOG(odds) # inverse-logit function; 
    pop.df$off = rbinom(n = n, size = max.off, prob = pr)
    
    pop.df$off.abuns = pop.df$off/sum(pop.df$off) # sum(pop.df$off.abuns)
    pop.df$w = pop.df$off.abuns/pop.df$rel.abuns
    pop.df$z.off = pop.df$z.par + rnorm(n = n, mean = mean.o, sd = sd.o)
    
    # Difference offspring parents 
    pop.df$delta.z.par.off = pop.df$z.off-pop.df$z.par
    
    # Selection part of the price equation which correspond to the covariance of fitness and the parental phenotype
    selection = cov(pop.df$w, pop.df$z.par)
    
    # Transmission part of the price equation 
    transmission =  mean(pop.df$w*(pop.df$delta.z.par.off))
    
    # Final result 
    delta_z = selection + transmission
    
    # Combine the results 
    rec.exp.price = c(rec.exp.price, delta_z)
    rec.select = c(rec.select, selection)
    rec.transm = c(rec.transm, transmission)
  }
  
  if (plothist) {
    old.par = par(no.readonly = TRUE)
    par(mfrow =c(2,3), mar = c(4,4,2,2))
    hist(rec.exp.price, main = "Expected difference"); abline(v = delta_z, col = "red")
    hist(rec.select, main = "Selection (covariance w and z)") ; abline(v = selection, col = "red")  
    hist(rec.transm, main = "Transmission parent-offspring")  ; abline(v = mean.o, col = "red")
    plot(pop.df[,"w"]~pop.df[,"z.par"], 
         pch = 19, 
         xlab = "Parent phenotype (z)",
         ylab = "Parent fitness (w)",
         main = "Simulated population",
         col = scales::alpha(colour = "black", .5)
         )#, 
         # cex = pop.df$off.abuns*100)
    abline(v = mean(pop.df$z.par), h = mean(pop.df$w), lty = 3)
    
    plot(pop.df$off~pop.df$z.par,
         col = scales::alpha(colour = "black", .5), pch = 19, 
         xlab = "Parent phenotype (z)",
         ylab = "Number of offspring",
         main = "Simulated population",
    ); abline(v = mean(pop.df$z.par), h = mean(pop.df$off), lty = 3)
    
    # Get parent offspring regression 
    lm.out = lm(pop.df$z.off~pop.df$z.par)
    plot(pop.df$z.off~pop.df$z.par,
         col = scales::alpha(colour = "black", .5), pch = 19, 
         xlab = "Parent phenotype (z)",
         ylab = "Offspring phenotype (z)",
         main = "Simulated population",
    ); abline(lm.out)
    abline(v = mean(pop.df$z.par), h = mean(pop.df$z.off), lty = 3)
    
    par(old.par)
  }
 return(list(delta.exp.price = rec.exp.price, 
             selection = rec.select, 
             transmission = rec.transm,
             pop.df = pop.df)) 
}
out.price = price.theorem(iter = 500, 
                          n = 100, beta1 = .50, 
                          mean.p = 10, mean.o = 0, 
                          max.off = 2,
                          sd.p = 1, sd.o = 1, seed = NULL)
round(out.price$selection[1:10],2)
round(out.price$transmission[1:10],2)
out.price$pop.df
