library(tidyverse)
n.sub = 4
n = 400 * n.sub
# set.seed(12)
pop.geno = data.frame(geno = rbinom(n = n,size = 2, prob = .1), 
           sub.pop = factor(rep(x = LETTERS[1:n.sub], each = n/n.sub)))
probs = c(.1,.2,.8,.75)
# pop.geno$geno = rep(x = c(1,2,0,1), each = n/n.sub)
pop.geno$geno = ifelse(pop.geno$sub.pop == "A",rbinom(n = n,size = 2, prob = probs[1]),pop.geno$geno)
pop.geno$geno = ifelse(pop.geno$sub.pop == "B",rbinom(n = n,size = 2, prob = probs[2]),pop.geno$geno)
pop.geno$geno = ifelse(pop.geno$sub.pop == "C",rbinom(n = n,size = 2, prob = probs[3]),pop.geno$geno)
pop.geno$geno = ifelse(pop.geno$sub.pop == "D",rbinom(n = n,size = 2, prob = probs[4]),pop.geno$geno)

geno.fq = table(pop.geno$geno)
tot.alleles = sum(geno.fq)*2
p.count = geno.fq[1]*2 + geno.fq[2]
q.count = geno.fq[3]*2 + geno.fq[2]
p = p.count/tot.alleles
q = q.count/tot.alleles

ht = 2*p*q

geno.table = pop.geno %>% 
  group_by(sub.pop) %>% 
  summarise(h0 = length(which(geno==0)),
            h1 = length(which(geno==1)),
            h2 = length(which(geno==2))) %>% 
  mutate(p.count = 2*h0+h1,
         q.count = 2*h2+h1,
         tot.all = p.count+q.count, 
         p = p.count/tot.all,
         q = q.count/tot.all,
         hs = 2*p*q) %>% 
  as.data.frame()

fst.cal <- function(hs, ht) {
  fst = (ht-(hs))/ht
  return(fst)
}
# fst.cal(hs = 0.1424, ht = 0.2371)

fst.cal(hs = mean(geno.table$hs), ht = ht)

###
hert= c(0.4893,
        0.4058,
        0.5000,
        0.4507,
        0.4216,
        0.4482,
        0.0000,
        0.0620,
        0.0139,
        0.0159,
        0.0100,
        0.0178,
        0.0100,
        0.0198,
        0.1268,
        0.0040,
        0.0080,
        0.2202,
        0.1895,
        0.3476,
        0.4842,
        0.0276)

sum(hert)/30
