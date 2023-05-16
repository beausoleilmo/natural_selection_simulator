# Description  ------------------------------------------------------------
#### ### ### ## #### ### ### ## #### ### ### ## 
# Created by Marc-Olivier Beausoleil
# Why: 
# Requires:
# NOTES: 
# Reference : 
#### ### ### ## #### ### ### ## #### ### ### ## 
set.seed(12345)

iter = 10

for (i in 1:iter) {
  # Population size 
  n = 1e2
  
  # Number of chromosomes 
  nchr = 2*n
  
  # Allele frequency
  alldfq = .50
  
  # Get alleles based on allele frequency for that population 
  alleles = rbinom(nchr, 1, alldfq)
  
  # Function to get expected number of genotypes 
  genohwe = function(q, n = NULL) {
    p = 1 - q
    hop = p^2
    he = 2*p*q
    hoq = q^2
    geno = cbind(hop,he,hoq)
    return((geno))
  }
  
  # Constructing genotypes
  # of all the alleles take a random sample from the chromosomes
  pos.c1 = sample(1:nchr,size = n)
  pos.c2 = c(1:nchr)[!(c(1:nchr) %in% pos.c1)]
  # sort(c(pos.c1,pos.c2))
  c1 = alleles[pos.c1] # store each allele for a chromosome 
  c2 = alleles[pos.c2]
  
  # Construct individual genotypes
  genotypes = cbind(c1,c2) |> as.data.frame()
  
  # Genotypes 0, 1 or 2 meaning homozygous '0,0', heterozygous '0,1' or '1,0', and homozygous '1,1' 
  genotypes$g = apply(genotypes,1,sum)
  
  # count the number of each genotype 
  (genotypes.table = table(genotypes$g))
  
  # From the allele frequency, find the expected number of genotypes 
  (expected.geno = genohwe(alldfq) * n)
  
  # This will let you know how 'far' we are from the expected: 
  # 100% means they are the same. 
  # <100% means less genotypes than expected, 
  # >100% means more than expected 
  check.missing.genotype <- function(data) {
    if (length(data) <3) {
      if (length(which(names(data) %in% "0")) == 0) {
        data = rbind(`0` = 0, as.matrix(data))
      }
    }
    if (length(data) <3) {
      if (length(which(names(data) %in% "1")) == 0) {
        data = rbind(`1` = 0, as.matrix(data))
      }
    }
    if (length(data) <3) {
      if (length(which(names(data) %in% "2")) == 0) {
        data = rbind(`2` = 0, as.matrix(data))
      }
    }  
    return(data)
  }
  genotypes.table = check.missing.genotype(genotypes.table)
  t(as.matrix(genotypes.table))/expected.geno *100
  
  table(alleles) / length(alleles)
  
  # Make the individuals reproduce together (selfing is allowed here...)
  repro = cbind(genotypes[,1:2], # 1 parent 
                genotypes[sample(x = 1:nrow(genotypes)), 1:2]) # second parent 
  
  # Select columns (e.g., equivalent of Mendelian segregation )
  col.sel1 = c(rbinom(nrow(repro), 1, 0.5)+1)
  col.sel2 = c(rbinom(nrow(repro), 1, 0.5)+3)
  
  # Make offspring table 
  offspring = cbind(
    c3 = repro[cbind(seq_along(col.sel1), col.sel1)], 
    c4 = repro[cbind(seq_along(col.sel2), col.sel2)]
  )
  
  offspring = as.data.frame(offspring)
  offspring$g = apply(offspring,1,sum)
  (offspring.table = table(offspring$g))
  offspring.table = check.missing.genotype(offspring.table)
  
  t(as.matrix(offspring.table))/expected.geno *100
  
  table(c(offspring$c3, offspring$c4))/ c(nrow(offspring)*2)
  
  offspring.table
  expected.geno
  
  
  repro$g1 = apply(repro[,1:2],1,sum)
  repro$g2 = apply(repro[,3:4],1,sum)
  repro$g3 = offspring$g
  all.fq1 = table(unlist(repro[,1:2]))/length(unlist(repro[,1:2]))
  all.fq2 = table(unlist(repro[,3:4]))/length(unlist(repro[,3:4]))
  
  gfq = c(all.fq1[1]^2,
          2*all.fq1[1]*all.fq1[2],
          all.fq1[2]^2)
  
  if (i==1) {
    matplot(x = seq(0,1, by = .01), 
            y = genohwe(q = seq(0,1, by = .01)), 
            main = paste("Population =",n, "Allele fq =", alldfq), 
            lty = 1, pch = 19, type = "l")
    abline(v = alldfq, lty = 3)
  }
  points(x = rep(all.fq1[2], 3), y = gfq, pch = 19, 
         cex= 4, col = scales::alpha(c("black","red","green"), .5))
  
}

