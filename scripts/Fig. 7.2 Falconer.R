# Description  ------------------------------------------------------------
##########################################################################################
# 
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created 2022
# Why: 
# Requires 
# NOTES: 
##########################################################################################

q = 1/4 
p = 1-q
a = 1
d = 3/4*a
c = -a
aa = 0
Aa = 1
AA = 2

pop.mean = p^2*a + 2*p*q*d + q^2*c
# mean.geno = p*a+q*d

alpha=a+d*(q-p)

A.avg = 2*q*alpha
D.avg = (q-p)*alpha
C.avg = -2*p*alpha

q2 = q^2
p2 = p^2
pq2 = 2*p*q

gen.eff= data.frame(geno = c(aa,Aa,AA),geno.val = c(c,d,a), avg.eff = c(C.avg,D.avg,A.avg))

plot(gen.eff$geno.val~gen.eff$geno, ylim = range(c(a,d,c,C.avg+d, D.avg+d, A.avg+d)), pch = 19)
abline(h = 0, lty = 2)
abline(h = d, lty = 3)
par(new = TRUE)
plot(gen.eff$geno,gen.eff$avg.eff+d, ylim = range(c(a,d,c,C.avg+d, D.avg+d, A.avg+d)), pch = 21, axes = FALSE, ann = FALSE)
axis(4, at = gen.eff$avg.eff+d, labels = gen.eff$avg.eff)
abline(lm(c(gen.eff$avg.eff+d)~gen.eff$geno))
points(1.5,pop.mean, pch =3, cex = 3)
