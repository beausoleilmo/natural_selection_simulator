# Description ----------- -------------------------------------------------------------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Created by Marc-Olivier Beausoleil                                                  #
# McGill University                                                                   #
# Created Thursday, February 1, 2023 
# Why:                                                                                #
  # Check a 3 allele HWE 
# Requires                                                                            #
# NOTES:                                                                              #
  # Computing HWE for 3 alleles 
  # The proportions still need to equate 1 
  # (p+q+r) = 1
  # To get the genotypes 
  # (p+q+r)^2 = 1^2
  # (p+q+r)(p+q+r) = 1
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
hwe.3alleles = function(p,q,r) {
  p2 = p^2
  q2 = q^2
  r2 = r^2  
  pq2 = 2*p*q
  pr2 = 2*p*r
  qr2 = 2*q*r
  geno.fq = c(p2,q2,r2,pq2,pr2,qr2)
  # return(geno.fq)
  return(list(p2=p2,q2=q2,r2=r2,pq2=pq2,pr2=pr2,qr2=qr2))
}

mysq = seq(0,1,by = .3)
alll2=(1/3)/seq(1,5,by = .5)
ref0 = c(rep(mysq[1], length(alll2)),
         rep(mysq[2], length(alll2)),
         rep(mysq[3], length(alll2)),
         rep(mysq[4], length(alll2)))
q.all = c((1)/seq(1,5,by = .5),
          (1/3)/seq(1,5,by = .5),
          (1/6)/seq(1,5,by = .5),
          (1/9)/seq(1,5,by = .5))
df.pq = data.frame(ref0,q.all)
df.pq$r = 1-(df.pq$ref0 + df.pq$q.all)
# exp.grsq = expand.grid(mysq,mysq)
hwe.out=hwe.3alleles(df.pq[,1],df.pq[,2],df.pq[,3])
mat.hwe=do.call(cbind,hwe.out)
all.fq =rep(mysq,nrow(mat.hwe)/length(mysq))

# matplot(x = ref0, y = mat.hwe, type = "l")

by = 0.01
x = seq(0,1, by = by)
i=1
rec.var = NULL
for (i in 1:length(x)) {
  remain = 1 - x[i] # What to separate between 2 other 
  y = seq(0, remain, by = by)
  remain2 = 1-(x[i] + y)
  rec.var = rbind(rec.var, cbind(rep(x[i], length(y)), y, remain2))
}
colnames(rec.var) <- c("p","q","r")

hwe.out=hwe.3alleles(rec.var[,1],rec.var[,2],rec.var[,3])
mat.hwe=do.call(cbind,hwe.out)

# matplot(x = rec.var[,1], y = mat.hwe[,4:6], type = "l", lwd = 3)
# legend("topleft", legend = c("p2","q2","r2","pq2","pr2","qr2"), col = 1:6, lty = 1)

# for (j in 1:ncol(mat.hwe)) {
#   if (j==1) {
#     scatterplot3d::scatterplot3d(x = rec.var[,1], y = rec.var[,2], mat.hwe[,j], pch = 19, zlim = c(0,1))
#   }
#   if (j>1) {
#     par(new = TRUE)
#     scatterplot3d::scatterplot3d(x = rec.var[,1], y = rec.var[,2], mat.hwe[,j], pch = 19,color = j, zlim = c(0,1))
#   }
# }
mat.all.geno = cbind(rec.var, mat.hwe) |> as.data.frame()

tab.rec.var = table(rec.var[,1])
tab.rec.var
mat.geno = matrix(data = NA, 
                  nrow = length(tab.rec.var), 
                  ncol = length(tab.rec.var))

mat.geno.list = NULL
w = 1
for (m in 1:ncol(mat.hwe)) {
  for (k in 1:length(tab.rec.var)) {
    mat.geno[k, 1:tab.rec.var[k]] <- mat.hwe[cumsum(tab.rec.var)[k],m]
    w = w + cumsum(tab.rec.var)[k]
  }
  mat.geno.list = c(mat.geno.list, list(mat.geno))
}

# Easy plotting for plotly 
library(plotly)
mat.all.geno.long = mat.all.geno %>% 
  pivot_longer(cols = c(p2,q2,r2,pq2,pr2,qr2))
oppp = 1

# Make up all matrices for plotly 
zmatp2 = mat.all.geno %>% 
  dplyr::select(-c(r,q2,r2,pq2,pr2,qr2)) %>% 
  pivot_wider(names_from = c(p), values_from = p2) %>% 
  dplyr::select(-q) %>% 
  as.data.frame()
zmatq2 = mat.all.geno %>% 
  dplyr::select(-c(p2,r,r2,pq2,pr2,qr2)) %>% 
  pivot_wider(names_from = c(p), values_from = q2) %>% 
  dplyr::select(-q) %>% 
  as.data.frame()
zmatr2 = mat.all.geno %>% 
  dplyr::select(-c(r,q2,p2,pq2,pr2,qr2)) %>% 
  pivot_wider(names_from = c(p), values_from = r2) %>% 
  dplyr::select(-q) %>% 
  as.data.frame()
zmatpq2 = mat.all.geno %>% 
  dplyr::select(-c(p2,r,r2,q2,pr2,qr2)) %>% 
  pivot_wider(names_from = c(p), values_from = pq2) %>% 
  dplyr::select(-q) %>% 
  as.data.frame()
zmatpr2 = mat.all.geno %>% 
  dplyr::select(-c(p2,r,r2,q2,pq2,qr2)) %>% 
  pivot_wider(names_from = c(p), values_from = pr2) %>% 
  dplyr::select(-q) %>% 
  as.data.frame()
zmatqr2 = mat.all.geno %>% 
  dplyr::select(-c(p2,r,r2,q2,pr2,pq2)) %>% 
  pivot_wider(names_from = c(p), values_from = qr2) %>% 
  dplyr::select(-q) %>% 
  as.data.frame()

zoom = 1.7
plotly::plot_ly(x = as.numeric(names(tab.rec.var)), 
                y = as.numeric(names(tab.rec.var)),
                z = as.matrix(zmatp2), 
              colorbar=list(title='p2'),
                type="surface", colors = colorRamp(c("yellow", "red"))) %>% 
  layout(title = 'Hardy-Weinberg 3 alleles',
         scene=list(
           xaxis = list(title = 'p',
                        autorange = TRUE),
           yaxis = list(title = 'q'),
           # autorange = "reversed"),
           zaxis=list(title = 'Proportion genotypes',nticks = 10,range = c(0,1)),# Number of years observed in population
           camera = list(eye = list(x = cos(0.8*pi)*zoom,
                                     y = sin(pi*1.3)*zoom#,
                                    # z= 2.00
                                    )
                         )
           )
         ) %>% 
  add_surface(x = as.numeric(names(tab.rec.var)), 
              y = as.numeric(names(tab.rec.var)),
              z = as.matrix(zmatq2), inherit = FALSE,
              colorbar=list(title='q2'),
              colorscale = list(c(0, 1), c("green", "blue")),
              type="surface") %>% 
  add_surface(x = as.numeric(names(tab.rec.var)), 
              y = as.numeric(names(tab.rec.var)),
              z = as.matrix(zmatr2), 
              colorbar=list(title='r2'),
              colorscale = list(c(0, 1), c("black", "white")),
              type="surface") %>% 
  add_surface(x = as.numeric(names(tab.rec.var)), 
              y = as.numeric(names(tab.rec.var)),
              z = as.matrix(zmatpq2), 
              colorbar=list(title='pq2'),
              colorscale = list(c(0, 1), c(viridis(5)[1], viridis(5)[5])),
              type="surface") %>% 
  add_surface(x = as.numeric(names(tab.rec.var)), 
              y = as.numeric(names(tab.rec.var)),
              z = as.matrix(zmatpr2), 
              colorbar=list(title='pr2'),
              colorscale = list(c(0, 1), c(rainbow(5)[1], rainbow(5)[5])),
              type="surface") %>% 
  add_surface(x = as.numeric(names(tab.rec.var)), 
              y = as.numeric(names(tab.rec.var)),
              z = as.matrix(zmatqr2), 
              colorbar=list(title='qr2'),
              colorscale = list(c(0, 1), c(terrain.colors(5)[1], terrain.colors(5)[5])),
              type="surface") 

table(mat.all.geno.long[mat.all.geno.long$name == "pq2",c("p","q")])

# Same as above but printing points instead of surface 
plotly::plot_ly() %>% 
  add_trace(data = mat.all.geno.long,
            x = ~p, y = ~q,
            z = ~value,
            # size = 1,
            type="scatter3d", 
            # type="scatter3d", 
            mode = "markers",
            showlegend = TRUE, inherit = TRUE,
            opacity = oppp,
            # marker = list(color = "red",opacity = 1, size = 12),
            # text = ~paste('ID: ', BANDFINAL),
            colors = viridis(length(unique(mat.all.geno.long$name))),
            color = ~name
            )  