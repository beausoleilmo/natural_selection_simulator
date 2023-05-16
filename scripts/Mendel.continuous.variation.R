# Description -------------------------------------------------------------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created March 18, 2023
# Why:
  # Mendel matrix
# Requires 
# NOTES: 
# Source: 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Figure 2.25  from Genetics from genes to genome p. 51 
# (Michael Goldberg and Janice Fischer and Leroy Hood and Leland Hartwell - Genetics From Genes to Genomes-McGraw Hill (2021))

library(tidyverse)

mendel.continuous <- function(nb.gene = 1, 
                              nb.alleles = 2) {
  # Get genes ID 
  gene = LETTERS[1:nb.gene]
  # Get allele ID 
  allele = c(nb.alleles-1):0
  # Make a grid of all gene-allele combinations 
  ga = expand.grid(gene, allele)
  # Make a new column that will paste all gene-allele combinations 
  ga = ga %>% mutate(ga = paste(Var1, Var2, sep = "")) %>%  arrange(Var1)
  # Separate the data for all genes so that I can loop them 
  out <- split( ga , f = ga$Var1 )
  # Make initial condition (first gene-allele)
  var.init = out[[1]]$ga
  # If nb.gene == 1, no need to loop! Obviously! 
  if (nb.gene > 1) {
    # From the second gene, an on... 
    for (i in 2:length(out)) {
      # Get the init values as first gene-allele
      tmp.g1 = var.init
      # Take the next gene-allele to combine
      tmp.g2 = out[[i]]$ga
      # Make a grid of the 2 gene-allele selected 
      gall = expand.grid(tmp.g1,tmp.g2)
      # Combine them again
      ga = gall %>% mutate(ga = paste(Var1, Var2, sep = ""))
      # Order to make it look good 
      ga = ga[order(ga$Var1),]
      # NEW var.init so that we ACCUMULATE the gene-allele 
      var.init = ga$ga
    }
  }
  # Make a punnet square type 
  punintended = outer(ga$ga, ga$ga, paste)
  # Rename matrix row/columns
  rownames(punintended) <- ga$ga
  colnames(punintended) <- ga$ga
  # Get only the numbers 
  nb.only = gsub('[[:alpha:]]+', '', punintended)
  # Get dimensions of data to reconstruct the matrix later 
  mdim = dim(nb.only)
  # Split the strings to get all the numbers 
  nb.list = strsplit(nb.only, "")
  
  # Reconstructing punnett square
  matpun = do.call(rbind, nb.list)[,order(c(gene," ", gene))][,2:(length(gene)*2+1)]
  # punnett in order 
  mypun = t(apply(matpun, 1, function(x) paste(sort(c(gene, gene)), x)))
  mypun2 = t(apply(mypun, 1, function(x) gsub(pattern = " ", replacement = "", x)))
  mypun3 = apply(mypun2, 1, function(x) paste(x, collapse = ""))
  pun.pretty.neat = matrix(unlist(mypun3), nrow = mdim[1], ncol = mdim[2])
  # Rename matrix row/columns
  rownames(pun.pretty.neat) <- ga$ga
  colnames(pun.pretty.neat) <- ga$ga
  
  # Sum of additive effect of alleles 
  nb.sum = lapply(nb.list,  function(x) sum(as.numeric(x), na.rm = TRUE))
  # Remake the matrix
  mmat = matrix(unlist(nb.sum), nrow = mdim[1], ncol = mdim[2])
  # Order the matrix to get symetrical data 
  mmat = mmat[,order(unlist(mmat[1,]), decreasing = T)]
  mmat = mmat[order(unlist(mmat[,1]), decreasing = T),]
  # Rename matrix row/columns
  rownames(mmat) <- ga$ga
  colnames(mmat) <- ga$ga
  
  # Plotting 
  par(mfrow=c(1,2), mar = c(2,2,1,1))
  # Make an image of the matrix to show the additivity of the genes and alleles 
  image(x = 1:mdim[1], y = 1:mdim[2], z = (mmat), 
        ylim = rev(range(1:mdim[2])), 
        asp = 1, xlab = "", ylab = "")
  txt.lab = expand.grid(1:mdim[1],1:mdim[2])
  text(x = txt.lab[,1], y = txt.lab[,2], labels = mmat)
  # mmat
  # Get the barplot height data 
  height.bp = table(unlist(nb.sum))
  # Plot the barplot of the number of times a certain additive effect is found (keep the bp data for adding labels)
  bp.dat = barplot(height.bp, 
                   ylim = c(0, max(height.bp)+(10^log10(max(height.bp))*.1)))
  # Add baseline to barplot 
  abline(h=0)
  
  # This is a way to add plotting information that SCALES nicely when adding more data 
  # 10^log10(max(height.bp))
  # max(height.bp)/.5
  text(bp.dat, height.bp+(10^log10(max(height.bp))*.05), 
       labels = round(height.bp/sum(height.bp)*100, digits = 1))
  
  # Return the data generated from the function 
  return(list(punnett.square.labels = pun.pretty.neat,
              additive.effects.alleles = mmat))
} 

nb.gene = 2
nb.alleles = 2
# DIMENSIONS 
nb.alleles^nb.gene
punintended = mendel.continuous(nb.gene = nb.gene, 
                  nb.alleles = nb.alleles)
dim(punintended$punnett.square.labels)
punintended$additive.effects.alleles
