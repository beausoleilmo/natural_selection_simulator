
# Source https://stackoverflow.com/questions/71547573/draw-a-vector-field-from-matrix-multiplication-r 


# loading libraries -------------------------------------------------------
library(fields) # for image.plot 
library(plotly)
library(plotrix)
library(raster)

library(tidyverse)
library(ggquiver)


# Functions ---------------------------------------------------------------
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


adaptive.land.geno = function(geno.mat, res = 0.01, 
                              x = NULL, # POSITION from 1 to 100 (placing a ball on the landscape)
                              vector.field.plot = FALSE, 
                              chrx=NULL, chry=NULL) {
  # if (vector.field.plot) {
  #   par(mfrow = c(2,1))
  # }
  # Sequence of X 
  seq.x = seq(from = 0, to = 1, by = res)
  # Make a matrix 
  space = outer(seq.x,seq.x,"*") 
  # Plot the matrix of all combinations of genotype frequencies
  # image.plot(space, ylim=c(1.05,-0.05), ylab= "Percentage of Chromosome EF of TD form", xlab= "Percentage of Chromosome CD of BL form")
  # Backup the data 
  space2 = space
  
  if (dim(geno.mat)[1]==1) {
    # calculate the average fitness for EVERY combination of frequency of 2 genotypes 
    # From Lamicchaney 
    fill.seq = seq.x
    for (i in 1:length(seq.x)) {
      # Calculate mean fitness 
      fill.seq[i] =   (geno.fit2/max(geno.fit2)) %*% all.p(1-seq.x[i])
    }
    plot(fill.seq~seq.x, pch = 19);abline(lm(fill.seq~seq.x), lwd = 4, col = "red")
    x = NULL
  } else {
    # calculate the average fitness for EVERY combination of frequency of 2 genotypes 
    for (i in 1:length(seq.x)) {
      for (j in 1:length(seq.x)) {
        # Calculate mean fitness 
        space[i,j] = all.p(1-seq.x[i]) %*% geno.mat %*% all.p(1-seq.x[j])
      }
    }
    # Show the result 
    round(t(space),3)
    
    if (is.null(chrx)) {  chrx = "A"  }
    if (is.null(chry)) {  chry = "B"  }
    
    # Transform the space
    new.space = t(space)
    image.plot(new.space, 
               ylim = rev(c(0,1)),
               # ylim=c( 1.01,-0.01), 
               ylab= paste("Percentage of Chromosome",chry, sep = " "),
               xlab= paste("Percentage of Chromosome",chrx, sep = " "))
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
    
  }
  
  dim.mat = dim(new.space)
  
  if(!is.null(x)){
    x1 = x[1]
    x2 = x[2]
    set.seed(1246)
    # Circle method 
    for (i in 1:100) {
      print(i)
      
      # add the starting point 
      if (i==1) {
        points(seq.x[x1],seq.x[x2], pch =19, col = "red")
        circl.dat = circle(seq.x[x1],seq.x[x2],.01)
      } else {
        circl.dat = circle(x1,x2,.01)
      }
      # recalculate the fitness 
      circl.dat$fit = rep(NA,length(circl.dat$x))
      # calculate the average fitness for EVERY combination of frequency of 2 genotypes 
      for (k in 1:length(circl.dat$x)) {
        # Calculate mean fitness 
        circl.dat$fit[k] = all.p(1-circl.dat$y[k]) %*% geno.mat %*% all.p(1-circl.dat$x[k])
      }
      x1 = circl.dat$x[which.max(circl.dat$fit)]
      x2 = circl.dat$y[which.max(circl.dat$fit)]
      points(circl.dat$x,circl.dat$y,cex = -log(circl.dat$fit/max(circl.dat$fit))*12)
      points(x1, x2, col = "green", pch = 19)
    }
  }
  # To get the 3D plane in an INTERACTIVE graph 
  xyz=cbind(expand.grid(x = seq.x,
                        y = seq.x),
            fit = as.vector(new.space))
  
  # Vector field on the Adaptive landscape ----------------------------------
  r <- raster::raster(
    space,
    xmn=range(seq.x)[1], xmx=range(seq.x)[2],
    ymn=range(seq.x)[1], ymx=range(seq.x)[2],
    crs=CRS("+proj=utm +zone=11 +datum=NAD83")
  )
  
  # Draw the adaptive landscape
  if (vector.field.plot) {
    raster2quiver(rast = r, aggregate = 2, colours = tim.colors(100)) 
  }
  
  return(list(new.space= new.space, xyz= xyz, raster = r))
}


# Genotype fitness matrix -------------------------------------------------
geno.fit = matrix(c(0.791,1.000,0.834,
                    0.670,1.006,0.901,
                    0.657,0.657,1.067), 
                  nrow = 3, 
                  ncol = 3,
                  byrow = T)
geno.fit.space = adaptive.land.geno(geno.mat = geno.fit, x = c(90,10))
geno.fit2 = matrix(c(30.0,53.1,73.7),
                   nrow = 1, 
                   ncol = 3,
                   byrow = T)
geno.fit2.space = adaptive.land.geno(geno.mat = geno.fit2)

genot.mat.fit = structure(c(0, 0.4, 0.206896551724138, 0.444444444444444, 0.55421686746988, 
                            0.323943661971831, 0.508771929824561, 0.522556390977444, 0.777777777777778), dim = c(3L, 3L))
geno.mat.fit.space = adaptive.land.geno(geno.mat = genot.mat.fit, x = c(50,50), 
                                        chrx = "LR761574.1_33199177",chry = "LR761574.1_33200027")

set.seed(1235)
fake.gen = matrix(data = runif(n = 9), nrow = 3, ncol = 3)
(fake.gen = matrix(data = c(0,3,0,
                            1,0,3.1,
                            2,2,0), nrow = 3, ncol = 3, byrow = T))
# fake.gen = matrix(data = 1:9, nrow = 3, ncol = 3, byrow = T)
tmp = adaptive.land.geno(geno.mat = fake.gen/max(fake.gen), x = c(90,12))
tmp = adaptive.land.geno(geno.mat = fake.gen/max(fake.gen), x = c(12,33))
tmp = adaptive.land.geno(geno.mat = fake.gen/max(fake.gen), x = c(12,41))

# Plotly 3D graph  --------------------------------------------------------
plotly.data = tmp
fit.mat = matrix(plotly.data$xyz$fit,
                 ncol = nrow(plotly.data$new.space), 
                 nrow = nrow(plotly.data$new.space), byrow = T)
plotly.data$fit.mat = fit.mat

# adds points 
# plot_ly(x = rev(plotly.data$xyz[,1]),
#         y = plotly.data$xyz[,2],
#         z = plotly.data$xyz[,3],
#         color = plotly.data$xyz[,3], type = "scatter3d", mode = "markers") %>%
# Add surface 
plot_ly() %>% 
  add_surface(data = plotly.data,  
              x=rev(plotly.data$xyz$x), 
              y=unique(plotly.data$xyz$y), 
              z=plotly.data$fit.mat) %>% 
  layout(yaxis  = list(range = c(1, 0), autorange = F, autorange="reversed"),showlegend = F)
