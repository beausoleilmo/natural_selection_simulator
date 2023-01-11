
# Source https://stackoverflow.com/questions/71547573/draw-a-vector-field-from-matrix-multiplication-r 

library(fields) # for image.plot 
library(plotly)
library(raster)

# Genotype fitness matrix -------------------------------------------------
geno.fit = matrix(c(0.791,1.000,0.834,
                    0.670,1.006,0.901,
                    0.657,0.657,1.067), 
                  nrow = 3, 
                  ncol = 3,
                  byrow = T)
# Resolution 
res = 0.01
# Sequence of X 
seq.x = seq(0,1,by = res)
# Make a matrix 
space = outer(seq.x,seq.x,"*") 

# Function to calculate the AVERAGE fitness for a given frequency of an allele to get the expected frequency of genotypes in a population
all.p <- function(p) { # Takes frequency of an allele in the population 
  
  if (length(p)>1) { # Has to be only 1 number 
    q = 1-p
    p2 = p^2
    pq2 = p*q*2
    q2 = q^2
    geno.df = data.frame(p2, pq2,q2)
    return((geno.df))
  } else {
    p2 = p^2 # Gets the AA 
    pq2 = 2*p*(1-p) # gets the Aa
    q2 = (1-p)^2 # Gets the aa 
    return(list=c(p2 = p2, pq2 = pq2, q2 = q2 # Return the values 
    )) 
  }
}

# Examples 
all.p(0)
all.p(1)
all.p(.5)
all.p(seq(0,1, by =.1))

# Plot the matrix of all combinations of genotype frequencies
image.plot(space,
           ylim=c(1.05,-0.05), 
           ylab= "Percentage of Chromosome EF of TD form",
           xlab= "Percentage of Chromosome CD of BL form")
# Backup the data 
space2 = space

# calculate the average fitness for EVERY combination of frequency of 2 genotypes 
for (i in 1:length(seq.x)) {
  for (j in 1:length(seq.x)) {
    # Calculate mean fitness 
    space[i,j] = all.p(1-seq.x[i]) %*% geno.fit %*% all.p(1-seq.x[j])
  }
}
# Show the result 
round(t(space),3)

# Transform the space
new.space = t(space)
image.plot(new.space, 
           ylim = rev(c(0,1)),
           # ylim=c( 1.01,-0.01), 
           ylab= "Percentage of Chromosome EF of TD (Tidbinbilla) form",
           xlab= "Percentage of Chromosome CD of BL (Blundell) form")
# Add the numbers to get a better sense of the average fitness values at each point 
by.text = 8
for (i in seq(1,length(seq.x),by = by.text)) {
  for (j in seq(1,length(seq.x),by = by.text)) {
    text(seq.x[i],seq.x[j],
         labels = round(new.space[i,j],4),
         cex = new.space[i,j]/2, 
         col = "black") # col = "gray70"
  }
}
# Add contour lines 
contour(new.space,ylim=c(1,0),add = T, nlevels = 50)

dim.mat = dim(new.space)
# set.seed(1236)
# set.seed(1239)
# set.seed(1240)
# set.seed(1241)
set.seed(1246)
x1 = sample(x = 1:dim.mat[1], size = 1) # 90
x2 = sample(x = 1:dim.mat[2], size = 1) # 90
x1 = which(round(seq.x,2)==.05)
x2 = which(round(seq.x,2)==.05)
# x1 = which(round(seq.x,2)==.90)
# x2 = which(round(seq.x,2)==.02)
# x1 = which(round(seq.x,2)==.90)
# x2 = which(round(seq.x,2)==.2)
# x1 = which(round(seq.x,2)==.3)
# x2 = which(round(seq.x,2)==.90)
new.space[x1,x2]

# Circle method 

# library(plotrix) # The function circle comes from this package 
circle = function (x, y, radius, nv = 100, border = NULL, col = NA, lty = 1, 
                   density = NULL, angle = 45, lwd = 1) 
{
  xylim <- par("usr")
  plotdim <- par("pin")
  ymult <- getYmult()
  angle.inc <- 2 * pi/nv
  angles <- seq(0, 2 * pi - angle.inc, by = angle.inc)
  if (length(col) < length(radius)) 
    col <- rep(col, length.out = length(radius))
  for (circle in 1:length(radius)) {
    xv <- cos(angles) * radius[circle] + x
    yv <- sin(angles) * radius[circle] * ymult + y
    polygon(xv, yv, border = border, col = col[circle], 
            lty = lty, density = density, angle = angle, lwd = lwd)
  }
  return(list(x = xv, y = yv))
}

geno.hwe = function(geno.q) {
  q = sqrt(geno.q)
  p=1-q
  p2 = p^2
  pq2 = 2*p*q
  q2 = q^2
  return(list=c(p2 = p2, pq2 = pq2, q2 = q2)) 
}

i=1
for (i in 1:100) {
  print(i)
  
  if (i==1) {
    points(seq.x[x1],seq.x[x2], pch =19, col = "red")
    circl.dat = circle(seq.x[x1],seq.x[x2],.01)
  } else {
    circl.dat = circle(x1,x2,.01)
  }
  circl.dat$fit = rep(NA,length(circl.dat$x))
  # calculate the average fitness for EVERY combination of frequency of 2 genotypes 
  for (k in 1:length(circl.dat$x)) {
    # Calculate mean fitness 
    circl.dat$fit[k] = all.p(1-circl.dat$y[k]) %*% geno.fit %*% all.p(1-circl.dat$x[k])
  }
  x1 = circl.dat$x[which.max(circl.dat$fit)]
  x2 = circl.dat$y[which.max(circl.dat$fit)]
  points(circl.dat$x,circl.dat$y,cex = -log(circl.dat$fit/max(circl.dat$fit))*12)
  points(x1, x2, col = "green", pch = 19)
}

# Square method (NOT OPTIMALLL)
# for (i in 1:1000) {
#   points(seq.x[x1],seq.x[x2], pch =19, col = "red")
#   # Fitness at that location 
#   new.space[x1,x2]
#   seq.x[c(x1+1,x2)]
#   right.pt = c(x1+1,x2)
#   left.pt = c(x1-1, x2)
#   bottom.pt = c(x1, x2+1)
#   top.pt = c(x1, x2-1)
#   dright.pt = c(x1+1, x2+1)
#   dleft.pt = c(x1-1, x2-1)
#   dbottom.pt = c(x1-1, x2+1)
#   dtop.pt = c(x1+1, x2-1)
#   
#   new.space[x1,x2]
#   # all.pos=rbind(top.pt,bottom.pt,right.pt,left.pt,dright.pt,dleft.pt,dbottom.pt,dtop.pt)
#   
#   pos.move = matrix(c(seq.x[top.pt[1]], seq.x[top.pt[2]],
#                       seq.x[bottom.pt[1]], seq.x[bottom.pt[2]],
#                       seq.x[right.pt[1]], seq.x[right.pt[2]],
#                       seq.x[left.pt[1]], seq.x[left.pt[2]],
#                       seq.x[dright.pt[1]], seq.x[dright.pt[2]],
#                       seq.x[dleft.pt[1]], seq.x[dleft.pt[2]],
#                       seq.x[dbottom.pt[1]], seq.x[dbottom.pt[2]],
#                       seq.x[dtop.pt[1]], seq.x[dtop.pt[2]]), nrow = 8, byrow = TRUE)
#   pos.move.coord = matrix(c(top.pt[1], top.pt[2],
#                             bottom.pt[1], bottom.pt[2],
#                             right.pt[1], right.pt[2],
#                             left.pt[1], left.pt[2],
#                             dright.pt[1], dright.pt[2],
#                             dleft.pt[1], dleft.pt[2],
#                             dbottom.pt[1], dbottom.pt[2],
#                             dtop.pt[1], dtop.pt[2]), nrow = 8, byrow = TRUE)
#   
#   dim.test = apply(pos.move.coord, 2, max)
#   if (dim.mat[1] <dim.test[1] | dim.mat[2] <dim.test[2]) {
#     stop("Bang! I think you just hit the wall!")  
#   }
#   fit.val.near = c(new.space[right.pt[1], right.pt[2]],
#                    new.space[left.pt[1], left.pt[2]],
#                    new.space[bottom.pt[1], bottom.pt[2]],
#                    new.space[top.pt[1], top.pt[2]],
#                    new.space[dright.pt[1], dright.pt[2]],
#                    new.space[dleft.pt[1], dleft.pt[2]],
#                    new.space[dbottom.pt[1], dbottom.pt[2]],
#                    new.space[dtop.pt[1], dtop.pt[2]])
#   
#   
#   # points(pos.move)
#   pos.move[which.max(fit.val.near),]
#   
#   x1 = pos.move.coord[which.max(fit.val.near),1]
#   x2 = pos.move.coord[which.max(fit.val.near),2]
# }

# Plotly 3D graph  --------------------------------------------------------
# To get the 3D plane in an INTERACTIVE graph 
xyz=cbind(expand.grid(seq.x,
                      seq.x),
          as.vector(new.space))

plot_ly(x = rev(xyz[,1]),y = xyz[,2],z = xyz[,3],
        color = xyz[,3], type = "scatter3d", mode = "markers") %>% 
  layout(yaxis  = list(range = c(1, 0), autorange = F, autorange="reversed"),showlegend = F)



# Vector field on the Adaptive landscape ----------------------------------
library(tidyverse)
library(ggquiver)
raster2quiver <- function(rast, aggregate = 50, colours = terrain.colors(6), contour.breaks = 200)
{
  names(rast) <- "z"
  quiv <- aggregate(rast, aggregate)
  terr <- terrain(quiv, opt = c('slope', 'aspect'))
  quiv$u <- -terr$slope[] * sin(terr$aspect[])
  quiv$v <- -terr$slope[] * cos(terr$aspect[])
  quiv_df <- as.data.frame(quiv, xy = TRUE)
  rast_df <- as.data.frame(rast, xy = TRUE)
  
  print(ggplot(mapping = aes(x = x, y = y, fill = z)) + 
          geom_raster(data = rast_df, na.rm = TRUE) + 
          geom_contour(data = rast_df, 
                       aes(z=z, color=..level..),
                       breaks = seq(0,3, length.out = contour.breaks), 
                       size = 1.4)+
          scale_color_gradient(low="blue", high="red")+
          geom_quiver(data = quiv_df, aes(u = u, v = v), vecsize = 1.5) +
          scale_fill_gradientn(colours = colours, na.value = "transparent") +
          theme_bw())
  
  return(quiv_df)
}

r <-raster(
  space,
  xmn=range(seq.x)[1], xmx=range(seq.x)[2],
  ymn=range(seq.x)[1], ymx=range(seq.x)[2],
  crs=CRS("+proj=utm +zone=11 +datum=NAD83")
)

# Draw the adaptive landscape
raster2quiver(rast = r, aggregate = 2, colours = tim.colors(100)) 
x1 = which(round(seq.x,2)==.05)
x2 = which(round(seq.x,2)==.05)
# x1 = which(round(seq.x,2)==.90)
# x2 = which(round(seq.x,2)==.02)
# x1 = which(round(seq.x,2)==.90)
# x2 = which(round(seq.x,2)==.2)
# x1 = which(round(seq.x,2)==.3)
# x2 = which(round(seq.x,2)==.90)
new.space[x1,x2]

# Circle method 

for (i in 1:100) {
  print(i)
  
  if (i==1) {
    points(seq.x[x1],seq.x[x2], pch =19, col = "red")
    circl.dat = circle(seq.x[x1],seq.x[x2],.01)
  } else {
    circl.dat = circle(x1,x2,.01)
  }
  circl.dat$fit = rep(NA,length(circl.dat$x))
  # calculate the average fitness for EVERY combination of frequency of 2 genotypes 
  for (k in 1:length(circl.dat$x)) {
    # Calculate mean fitness 
    circl.dat$fit[k] = all.p(1-circl.dat$y[k]) %*% geno.fit %*% all.p(1-circl.dat$x[k])
  }
  x1 = circl.dat$x[which.max(circl.dat$fit)]
  x2 = circl.dat$y[which.max(circl.dat$fit)]
  points(circl.dat$x,circl.dat$y,cex = -log(circl.dat$fit/max(circl.dat$fit))*12)
  points(x1, x2, col = "green", pch = 19)
}