# Description  ------------------------------------------------------------
#### ### ### ## #### ### ### ## #### ### ### ## 
# Theoretical fitness landscape 
# Drawn from multiple perspectives (linear, quadratic, normal, etc.)
# Created by Marc-Olivier Beausoleil
# 20??
# Why: 
  # This script generates 3D graphs with hypothetical fitness values for 2 traits 
  # This was originally made to make a point of what we expect to see from natural data 
# Requires:
# NOTES: 
# Reference : 
#### ### ### ## #### ### ### ## #### ### ### ## 

dir.create("output/images", recursive = TRUE)

pdf("output/images/test_fitland.pdf", 
    height = 5,width = 5)
library(rgl)
color = "green"
border = NA
par(mfrow=c(2,2))
# Kimura's Neutral theory "flat landscape" --------------------------------
x = seq(1,100,10)
y = x
# z = rep(0.5,100)
# z = matrix(0.5, nrow = 100, ncol = 100)
# flat = c(x,y,z)
# persp3d(x,y,z,zlim = c(0,1),aspect = c(1, 1,1), col = color,
#         xlab = "Trait 1", ylab = "Trait 2", zlab = "Fitness")
# surface3d(x,y,z, back = "lines")
# surface3d(x,y,z, front = "lines")

z <- outer(x, y, function(x,y){x/x*1/2})
persp(x, y, z, col = color, 
      xlab = "Trait 2", 
      ylab = "Trait 1", zlab = "Fitness", 
      theta = 60, d = 10,
      phi = 15,expand = 1,shade = 0.6,
      border = border,
      zlim=c(0,1))

### Linear 
x = seq(1,100,10)
y = x
# z = rep(0.5,100)
# z = matrix(0.5, nrow = 100, ncol = 100)
# flat = c(x,y,z)
# persp3d(x,y,z,zlim = c(0,1),aspect = c(1, 1,1), col = color,
#         xlab = "Trait 1", ylab = "Trait 2", zlab = "Fitness")
# surface3d(x,y,z, back = "lines")
# surface3d(x,y,z, front = "lines")

z <- outer(x, y, function(x,y){2*y})
persp(x, y, z, col = color, 
      xlab = "Trait 2", 
      ylab = "Trait 1", zlab = "Fitness", 
      theta = 60, d = 10,
      phi = 15,expand = 1,shade = 0.6,
      border = border,
      zlim=c(0,200))

### Linear 
x = seq(1,100,10)
y = x
# z = rep(0.5,100)
# z = matrix(0.5, nrow = 100, ncol = 100)
# flat = c(x,y,z)
# persp3d(x,y,z,zlim = c(0,1),aspect = c(1, 1,1), col = color,
#         xlab = "Trait 1", ylab = "Trait 2", zlab = "Fitness")
# surface3d(x,y,z, back = "lines")
# surface3d(x,y,z, front = "lines")

# z <- outer(x, y, function(x,y){-x+y+95})
# persp(x, y, z, col = color, 
#       xlab = "Trait 2", 
#       ylab = "Trait 1", zlab = "Fitness", 
#       theta = 60, d = 10,
#       phi = 15,expand = 1,shade = 0.6,
#       border = border,
#       zlim=c(0,200))

plot3d(x, y, z, col = color, 
      xlab = "Trait 2", 
      ylab = "Trait 1", zlab = "Fitness", 
      border = border,
      zlim=c(0,200))

x = seq(-10,10,0.1)
y = x

# z <- outer(x, y, function(x,y){x^2+y^2})
# persp(x, y, z, col = color, 
#       xlab = "Trait 2", 
#       ylab = "Trait 1", zlab = "Fitness", 
#       theta = 60, d = 10,
#       phi = 15,expand = 1,shade = 0.6,
#       border = border,
#       zlim=c(0,200))

z <- outer(x, y, function(x,y){200-(x^2+y^2)})
persp(x, y, z, col = color, 
      xlab = "Trait 2", 
      ylab = "Trait 1", zlab = "Fitness", 
      theta = 60, d = 10,
      phi = 15,expand = 1,shade = 0.6,
      border = border,
      zlim=c(0,200))


# points3d(x = x,y = y,z = z,col="green")



# One peak only (Fisher) --------------------------------------------------
x <- seq(-10, 10, length = 100)
y <- x
norm <- function(x,y, mean1 = 0, var1 = 2, mean2 = 0, var2 = 2) {
  1/sqrt(2*var1^2*1)*exp(-((x-mean1)^2/(2*var1^2))) *
    1/sqrt(2*var2^2*1)*exp(-((y-mean2)^2/(2*var2^2)))}
z <- outer(x, y, norm)
z[is.na(z)] <- 1
persp(x, y, z, col = color, 
      xlab = "Trait 2", 
      ylab = "Trait 1", zlab = "Fitness", 
      theta = 60,
      phi = 15,d = 10,
      border = border,
      shade = 0.6,
      zlim=c(0,.2))

# Library
library(plotly)
# Plot
p <- plot_ly(z = z, type = "surface")
p 
# save the widget
# library(htmlwidgets)
# saveWidget(p, file=paste0( getwd(), "/HtmlWidget/3dSurface.html"))

#### With environment
x <- seq(-10, 10, length = 100)
y <- x
norm <- function(x,y, mean1 = 0, var1 = 2, mean2 = 0, var2 = 2) {
  1/sqrt(2*var1^2*1)*exp(-((x-mean1)^2/(2*var1^2))) *
    1/sqrt(2*var2^2*1)*exp(-((y-mean2)^2/(2*var2^2)))}
z <- outer(x, y, norm)
z[is.na(z)] <- 1
persp(x, y, z, col = color, 
      xlab = "Trait 2", 
      ylab = "Trait 1", zlab = "Fitness", 
      theta = 60,
      phi = 15,d = 10,
      border = border,
      shade = 0.6,
      zlim=c(0,.2))
write.csv(z,"~/Desktop/z.csv",row.names = FALSE)
z2 = as.matrix(read.csv("~/Desktop/z.csv"))
z2[1,] = 0.13 # wall1 trait1 "min" trait 2
# z2[dim(z2)[1],] = 0.13 # wall2 trait1 "max" trait 2
z2[,dim(z2)[2]] = 0.13 # wall2 trait2 "max" trait 1
# z2[,1] = 0.13 # wall2 trait2 "min" trait 1
par(new=TRUE)
red.a = adjustcolor( "red", alpha.f = 0.999)
persp(x, y, z2 + 0.7e-1, col = red.a, 
      xlab = "", 
      ylab = "", zlab = "", 
      theta = 60,
      phi = 15,d = 10,
      border = border,
      shade = 0.6,
      zlim=c(0,.2))

x = seq(1,100,10)
y = x
z <- outer(x, y, function(x,y){x/x*1})
par(new=TRUE)
red.a = adjustcolor( "red", alpha.f = 0.1)
persp(x, y, z, col = red.a, 
      xlab = "", 
      ylab = "", zlab = "", 
      theta = 60, d = 10,
      phi = 15,expand = 1,shade = 0.6,
      border = border,
      zlim=c(0,1))



#### Only environment 
# x <- seq(-10, 10, length = 100)
# y <- x
# red.a = adjustcolor( "red", alpha.f = 1)
# persp(x, y, z2 + 0.7e-1, col = red.a, 
#       xlab = "Trait 2", 
#       ylab = "Trait 1", zlab = "Fitness", 
#       theta = 60,
#       phi = 15,d = 10,
#       border = border,
#       shade = 0.6,
#       zlim=c(0,.2))

x = seq(1,100,10)
y = x
z <- outer(x, y, function(x,y){x/x*1})
par(new=TRUE)
red.a = adjustcolor( "red", alpha.f = 0.1)
persp(x, y, z, col = red.a, 
      xlab = "", 
      ylab = "", zlab = "", 
      theta = 60, d = 10,
      phi = 15,expand = 1,shade = 0.6,
      border = border,
      zlim=c(0,1))




# Double peaks ------------------------------------------------------------
x <- seq(-10, 10, length = 100)
y <- x
# var = 10
# mean = 0 
norm <- function(x,y, mean1 = 3, var1 =2,mean2 = -3, var2 =2) {
  1/sqrt(2*var1^2*1)*exp(-((y-mean2)^2/(2*var2^2))) * 
    1/sqrt(2*var2^2*1)*exp(-((x-mean2)^2/(2*var1^2)))+
  1/sqrt(2*var2^2*1)*exp(-((y-mean1)^2/(2*var2^2))) * 
    1/sqrt(2*var1^2*1)*exp(-((x-mean1)^2/(2*var1^2)))      }

z <- outer(x, y, norm)
z[is.na(z)] <- 1
# open3d()
# bg3d("white")
# material3d(col = "black")
# persp3d(x, y, z, aspect = c(1, 1, 0.5), col = color,
#         xlab = "Trait 2", ylab = "Trait 1", zlab = "Fitness")

persp(x, y, z, col = color, 
      xlab = "Trait 2", 
      ylab = "Trait 1", 
      zlab = "Fitness", 
      theta = 60,
      phi = 15,d = 10,
      border = border,
      shade = 0.6,
      zlim=c(0,.2))




# Independent fitness functions  ------------------------------------------
x <- seq(-10, 10, length = 100)
y <- x

norm <- function(x,y, mean1 = 3, var1 =2,mean2 = -3, var2 =2) {
  1/sqrt(2*var1^2*1)*exp(-((y-mean2)^2/(2*var2^2))) * 
    1/sqrt(2*var2^2*1)*exp(-((x-mean2)^2/(2*var1^2)))+
    1/sqrt(2*var2^2*1)*exp(-((y-mean1)^2/(2*var2^2))) * 
    1/sqrt(2*var1^2*1)*exp(-((x-mean1)^2/(2*var1^2)))      }

z <- outer(x, y, norm)
z[is.na(z)] <- 1

rbPal <- colorRampPalette(c('red','blue'))
Col <- rbPal(10)[as.numeric(cut(x = x,breaks = 2))]

persp(x, y, z, 
      col = "blue", 
      xlab = "Trait 2", 
      ylab = "Trait 1", 
      zlab = "Fitness", 
      theta = 60,
      phi = 15,d = 10,
       border = border,
      shade = 0.6,
      zlim=c(0,.2))

z2 <- z
for(i in 1:100){
  for(j in 1:(100-i)){
    z2[i,j] <- NA
  }
}



par(new=T)
graphics::persp(x, y, z2, col = "yellow", 
                xlab = "Trait 2", 
                ylab = "Trait 1", 
                zlab = "Fitness", 
                theta = 60,
                phi = 15,d = 10,
                 border = border,
                shade = 0.6,
                zlim=c(0,.2))



# Rugged ------------------------------------------------------------------
x <- seq(-10, 10, length = 100)
y <- x
var = 10
mean = 0 
# norm <- function(x,y, mean1 = -1, var1 =1.5,
#                  mean2 = rnorm(1,i,0.5), var2 =1,
#                  mean3 = 4, var3 =1.8,
#                  mean4 = 0, var4 =1.6,
#                  mean5 = -4, var5 =1.1,
#                  mean6 = 3, var6 =0.9,
#                  mean7 = -2, var7 =1.5,
#                  mean8 = 3, var8 =1.2,
#                  mean9 = 1, var9 =1.1,
#                  mean10 = -4, var10 =1.2) {
#   1/sqrt(2*var1^2*1)*exp(-((x-mean1)^2/(2*var1^2))) *
#     1/sqrt(2*var2^2*1)*exp(-((y-mean2)^2/(2*var2^2))) + 
#     1/sqrt(2*var3^2*1)*exp(-((x-mean3)^2/(2*var1^2))) *
#     1/sqrt(2*var4^2*1)*exp(-((y-mean4)^2/(2*var2^2))) +
#   1/sqrt(2*var5^2*1)*exp(-((x-mean5)^2/(2*var1^2))) *
#   1/sqrt(2*var6^2*1)*exp(-((y-mean6)^2/(2*var2^2))) +
#   1/sqrt(2*var7^2*1)*exp(-((x-mean7)^2/(2*var1^2))) *
#     1/sqrt(2*var8^2*1)*exp(-((y-mean8)^2/(2*var2^2))) +
#   1/sqrt(2*var9^2*1)*exp(-((x-mean9)^2/(2*var1^2))) *
#     1/sqrt(2*var10^2*1)*exp(-((y-mean10)^2/(2*var2^2))) }
vec = rep(1,10)
png(filename = "output/images/fit.land.example.png",
    width = 10,height = 6,units = "in", res = 300, pointsize = 12, bg = "white")
set.seed(1214)
for(i in 2){
norm <- function(x,y, 
                 mean1 = rnorm(1,vec[i],2), var1 =1.5,
                 mean2 = rnorm(1,vec[i],2), var2 =1.5,
                 mean3 = rnorm(1,vec[i],5), var3 =1.8,
                 mean4 = rnorm(1,vec[i],2), var4 =1.8,
                 mean5 = rnorm(1,vec[i],1), var5 =1.1,
                 mean6 = rnorm(1,vec[i],3), var6 =1.1,
                 mean7 = rnorm(1,vec[i],4), var7 =1.5,
                 mean8 = rnorm(1,vec[i],4), var8 =1.5,
                 mean9 = rnorm(1,vec[i],1), var9 =1.1,
                 mean10 = rnorm(1,vec[i],3), var10 =1.1) {
  1/sqrt(2*var1^2*1)*exp(-((x-mean1)^2/(2*var1^2))) *
    1/sqrt(2*var2^2*1)*exp(-((y-mean2)^2/(2*var2^2))) + 
    1/sqrt(2*var3^2*1)*exp(-((x-mean3)^2/(2*var1^2))) *
    1/sqrt(2*var4^2*1)*exp(-((y-mean4)^2/(2*var2^2))) +
    1/sqrt(2*var5^2*1)*exp(-((x-mean5)^2/(2*var1^2))) *
    1/sqrt(2*var6^2*1)*exp(-((y-mean6)^2/(2*var2^2))) +
    1/sqrt(2*var7^2*1)*exp(-((x-mean7)^2/(2*var1^2))) *
    1/sqrt(2*var8^2*1)*exp(-((y-mean8)^2/(2*var2^2))) +
    1/sqrt(2*var9^2*1)*exp(-((x-mean9)^2/(2*var1^2))) *
    1/sqrt(2*var10^2*1)*exp(-((y-mean10)^2/(2*var2^2)))}

z <- outer(x, y, norm)
z[is.na(z)] <- 1
# open3d()
# bg3d("white")
# material3d(col = "black")
# persp3d(x, y, z, aspect = c(1, 1, 0.5), col = color, contour = TRUE, #scale = FALSE, # front="line",
#         xlab = "Trait 1", ylab = "Trait 2", zlab = "Fitness",  ticktype="detailed",box=TRUE, axes=TRUE)

# Create a function interpolating colors in the range of specified colors
flame.colors <- colorRampPalette( c("yellow", "red") )
# Generate the desired number of colors from this palette
nbcol <- 100
color <- flame.colors(nbcol)
nrz <- nrow(z)
ncz <- ncol(z)
# Compute the z-value at the facet centres
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
# Recode facet z-values into color indices
facetcol <- cut(zfacet, nbcol)

par(mfrow = c(1,2), mar = c(3,3,1,1))
cex.labs = 2
p1 = persp(x, y, z, 
           col = color[facetcol],
           # col = color, 
      xlab = "",#"Trait 2", 
      ylab = "",#"Trait 1", 
      zlab = "",#"Fitness", 
      theta = 60,
      phi = 15,d = 10,
      border = border,
      shade = 0.6,
      zlim=c(0,1))
# persp3D(x, y, z, aspect = c(1, 1, 0.5), col = color, contour = FALSE, scale = FALSE, # front="line",
#         xlab = "Trait 1", ylab = "Trait 2", zlab = "Fitness",  ticktype="detailed",box=FALSE, axes=TRUE)
labels <- c('Trait 2')
x.axis = mean(x) 
min.x <- min(x) 
max.x <- max(x) 
y.axis <- mean(z)
min.y <- min(y) 
max.y <- max(y) 
z.axis <- seq(-100, 100, by=25)
min.z <- min(z) 
max.z <- max(z) 

label.pos <- trans3d(x = x.axis-4, y = (min.y - 4.0), z = min.z, pmat = p1)
text(label.pos$x, label.pos$y, labels=c('Trait 2'), adj=c(0, NA), srt=330, cex=cex.labs)
label.pos <- trans3d(x = (min.x+ 26.0), y = (y.axis -8.5), z = min.z, pmat = p1)
text(label.pos$x, label.pos$y, labels=c('Trait 1'), adj=c(0, NA), srt=14, cex=cex.labs)

label.pos <- trans3d(x = x.axis-55, y = y.axis+15, z = min.z+.1, pmat = p1)
text(label.pos$x, label.pos$y, labels=c('Fitness'), adj=c(0, NA), srt=270, cex=cex.labs)

}

# par(mfrow = c(1,1), mar = c(4,4,1,1))
image(x,y,z, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
title(xlab = "Trait 2", ylab = "Trait 1", line = 1, cex.lab = cex.labs)
dev.off()



vec = 1:100
for(i in 1){
  norm <- function(x,y, 
                   mean1 = -5, var1 =2.,
                   mean2 = -2, var2 =1.9,
                   mean3 = 3, var3 =1.4,
                   mean4 = -3, var4 =1.8,
                   mean5 = -3, var5 =1.3,
                   mean6 = -2, var6 =1.3,
                   mean7 = 2, var7 =1.8,
                   mean8 = 2, var8 =1.5,
                   mean9 = -2, var9 =1.3,
                   mean10 = 3, var10 =1.8) {
    1/sqrt(2*var1^2*1)*exp(-((x-mean1)^2/(2*var1^2))) *log(abs(rnorm(x,10)/2)) *
      1/sqrt(2*var2^2*1)*exp(-((y-mean2)^2/(2*var2^2))) + 
      1/sqrt(2*var3^2*1)*exp(-((x-mean3)^2/(2*var1^2))) *log(abs(rnorm(x,10))) *
      1/sqrt(2*var4^2*1)*exp(-((y-mean4)^2/(2*var2^2))) +
      1/sqrt(2*var5^2*1)*exp(-((x-mean5)^2/(2*var1^2))) *log(abs(rnorm(x,5))) *
      1/sqrt(2*var6^2*1)*exp(-((y-mean6)^2/(2*var2^2))) +
      1/sqrt(2*var7^2*1)*exp(-((x-mean7)^2/(2*var1^2))) *log(abs(rnorm(x,10))) *
      1/sqrt(2*var8^2*1)*exp(-((y-mean8)^2/(2*var2^2))) +
      1/sqrt(2*var9^2*1)*exp(-((x-mean9)^2/(2*var1^2))) *log(abs(rnorm(x,10))) *
      1/sqrt(2*var10^2*1)*exp(-((y-mean10)^2/(2*var2^2))) + 
      1/sqrt(2*var3^2*1)*exp(-((x-mean1)^2/(2*var1^2))) *log(abs(rnorm(x,10))) *
      1/sqrt(2*var4^2*1)*exp(-((y-mean2)^2/(2*var2^2))) +
      1/sqrt(2*var3^2*1)*exp(-((x-mean3)^2/(2*var1^2))) *log(abs(rnorm(x,10))) *
      1/sqrt(2*var4^2*1)*exp(-((y-mean4)^2/(2*var2^2))) +
    1/sqrt(2*var3^2*1)*exp(-((x-mean8)^2/(2*var1^2))) *log(abs(rnorm(x,10))) *
      1/sqrt(2*var4^2*1)*exp(-((y-mean8)^2/(2*var2^2))) +
      1/sqrt(2*var3^2*1)*exp(-((x-mean9)^2/(2*var1^2))) *log(abs(rnorm(x,10))) *
      1/sqrt(2*var4^2*1)*exp(-((y-mean9)^2/(2*var2^2))) }
  z <- outer(x, y, norm)
  # z * log(abs(x+1)) * sin(x)
  z[is.na(z)] <- 1
  # open3d()
  # bg3d("white")
  # material3d(col = "black")
  # persp3d(x, y, z, aspect = c(1, 1, 0.5), col = color, contour = TRUE, #scale = FALSE, # front="line",
  #         xlab = "Trait 1", ylab = "Trait 2", zlab = "Fitness",  ticktype="detailed",box=TRUE, axes=TRUE)

  persp(x, y, z, col = color,
        xlab = "Trait 2",
        ylab = "Trait 1", zlab = "Fitness",
        theta = 60,
        phi = 15,d = 10,
        border = border,
        shade = 0.6#,
        #zlim=c(0,10)
        )

  
    # persp3D(x, y, z, aspect = c(1, 1, 0.5), col = color, contour = FALSE, scale = FALSE, # front="line",
  #         xlab = "Trait 1", ylab = "Trait 2", zlab = "Fitness",  ticktype="detailed",box=FALSE, axes=TRUE)
}

# write.csv(z, "~/Desktop/treedimfitland.csv")
# z = read.csv("~/Desktop/treedimfitland.csv")

par(mfrow = c(1,1))

# Double peaks variation through years ------------------------------------------------------------
x <- seq(-10, 10, length = 100)
y <- x
# var = 10
# mean = 0 
norm <- function(x,y, mean1 = 3, var1 =1.7,mean2 = -3, var2 =2.3) {
  1/sqrt(2*var2^2*1)*exp(-((y-mean2)^2/(2*var2^2))) * 
    1/sqrt(2*var2^2*1)*exp(-((x-mean2)^2/(2*var2^2)))+
    1/sqrt(2*var1^2*1)*exp(-((y-mean1)^2/(2*var1^2))) * 
    1/sqrt(2*var1^2*1)*exp(-((x-mean1)^2/(2*var1^2)))      }

z <- outer(x, y, norm)
z[is.na(z)] <- 1
# open3d()
# bg3d("white")
# material3d(col = "black")
# persp3d(x, y, z, aspect = c(1, 1, 0.5), col = color,
#         xlab = "Trait 2", ylab = "Trait 1", zlab = "Fitness")

persp(x, y, z, col = color, 
      xlab = "Trait 2", 
      ylab = "Trait 1", 
      zlab = "Fitness", 
      theta = 60,
      phi = 15,d = 10,
      border = border,
      shade = 0.6,
      zlim=c(0,.2))


x <- seq(-10, 10, length = 100)
y <- x
# var = 10
# mean = 0 
norm <- function(x,y, mean1 = 3, var1 =2.3,mean2 = -3, var2 =1.7) {
  1/sqrt(2*var2^2*1)*exp(-((y-mean2)^2/(2*var2^2))) * 
    1/sqrt(2*var2^2*1)*exp(-((x-mean2)^2/(2*var2^2)))+
    1/sqrt(2*var1^2*1)*exp(-((y-mean1)^2/(2*var1^2))) * 
    1/sqrt(2*var1^2*1)*exp(-((x-mean1)^2/(2*var1^2)))      }

z <- outer(x, y, norm)
z[is.na(z)] <- 1
# open3d()
# bg3d("white")
# material3d(col = "black")
# persp3d(x, y, z, aspect = c(1, 1, 0.5), col = color,
#         xlab = "Trait 2", ylab = "Trait 1", zlab = "Fitness")

persp(x, y, z, col = color, 
      xlab = "Trait 2", 
      ylab = "Trait 1", 
      zlab = "Fitness", 
      theta = 60,
      phi = 15,d = 10,
      border = border,
      shade = 0.6,
      zlim=c(0,.2))


dev.off()


# 5 peaks  ----------------------------------------------------------------
par(mar = c(0,0,0,0))
x <- seq(-15, 15, length = 150)
y <- x
norm <- function(x,y, 
                 mean1 =  9,   var1 = 2,
                 mean2 = -9,   var2 = 2,
                 mean3 = -2.5, var3 = 2.3,
                 mean4 =  2.5, var4 = 2, 
                 mean5.1 = -6, var5.1 = 2.5,
                 mean5.2 =  6, var5.2 = 2.5) {
    1/sqrt(2*var2^2*1)*exp(-((y-mean2)^2/(2*var2^2))) * 
    1/sqrt(2*var2^2*1)*exp(-((x-mean2)^2/(2*var2^2))) +
    1/sqrt(2*var1^2*1)*exp(-((y-mean1)^2/(2*var1^2))) * 
    1/sqrt(2*var1^2*1)*exp(-((x-mean1)^2/(2*var1^2))) + 
    
    1/sqrt(2*var3^2*1)*exp(-((y-mean3)^2/(2*var3^2))) * 
    1/sqrt(2*var3^2*1)*exp(-((x-mean3)^2/(2*var3^2))) +
    1/sqrt(2*var4^2*1)*exp(-((y-mean4)^2/(2*var4^2))) * 
    1/sqrt(2*var4^2*1)*exp(-((x-mean4)^2/(2*var4^2))) +
    
    1/sqrt(2*var5.1^2*1)*exp(-((y-mean5.1)^2/(2*var5.2^2))) * 
    1/sqrt(2*var5.2^2*1)*exp(-((x-mean5.2)^2/(2*var5.2^2)))}

z <- outer(x, y, norm)
z[is.na(z)] <- 1

persp(x, y, z, col = color, 
      xlab = "Trait 2", 
      ylab = "Trait 1", 
      zlab = "Fitness", 
      theta = 60,
      phi = 15,d = 10,
      border = border,
      shade = 0.6,
      zlim=c(0,.2))



x <- seq(-15, 15, length = 150)
y <- x
norm <- function(x,y, 
                 mean1 =  9,   var1 = 2.2,
                 mean2 = -9,   var2 = 1.8,
                 mean3 = -2.5, var3 = 2,
                 mean4 =  2.5, var4 = 2.3, 
                 mean5.1 = -6, var5.1 = 2.2,
                 mean5.2 =  6, var5.2 = 2.2) {
  1/sqrt(2*var2^2*1)*exp(-((y-mean2)^2/(2*var2^2))) * 
    1/sqrt(2*var2^2*1)*exp(-((x-mean2)^2/(2*var2^2))) +
    1/sqrt(2*var1^2*1)*exp(-((y-mean1)^2/(2*var1^2))) * 
    1/sqrt(2*var1^2*1)*exp(-((x-mean1)^2/(2*var1^2))) + 
    
    1/sqrt(2*var3^2*1)*exp(-((y-mean3)^2/(2*var3^2))) * 
    1/sqrt(2*var3^2*1)*exp(-((x-mean3)^2/(2*var3^2))) +
    1/sqrt(2*var4^2*1)*exp(-((y-mean4)^2/(2*var4^2))) * 
    1/sqrt(2*var4^2*1)*exp(-((x-mean4)^2/(2*var4^2))) +
    
    1/sqrt(2*var5.1^2*1)*exp(-((y-mean5.1)^2/(2*var5.2^2))) * 
    1/sqrt(2*var5.2^2*1)*exp(-((x-mean5.2)^2/(2*var5.2^2)))}

z <- outer(x, y, norm)
z[is.na(z)] <- 1

persp(x, y, z, col = color, 
      xlab = "Trait 2", 
      ylab = "Trait 1", 
      zlab = "Fitness", 
      theta = 60,
      phi = 15,d = 10,
      border = border,
      shade = 0.6,
      zlim=c(0,.2))


x <- seq(-15, 15, length = 150)
y <- x
norm <- function(x,y, 
                 mean1 =  9,   var1 = 1.8,
                 mean2 = -9,   var2 = 2,
                 mean3 = -2.5, var3 = 1.8,
                 mean4 =  2.5, var4 = 2.1, 
                 mean5.1 = -6, var5.1 = 2.5,
                 mean5.2 =  6, var5.2 = 2.5) {
  1/sqrt(2*var2^2*1)*exp(-((y-mean2)^2/(2*var2^2))) * 
    1/sqrt(2*var2^2*1)*exp(-((x-mean2)^2/(2*var2^2))) +
    1/sqrt(2*var1^2*1)*exp(-((y-mean1)^2/(2*var1^2))) * 
    1/sqrt(2*var1^2*1)*exp(-((x-mean1)^2/(2*var1^2))) + 
    
    1/sqrt(2*var3^2*1)*exp(-((y-mean3)^2/(2*var3^2))) * 
    1/sqrt(2*var3^2*1)*exp(-((x-mean3)^2/(2*var3^2))) +
    1/sqrt(2*var4^2*1)*exp(-((y-mean4)^2/(2*var4^2))) * 
    1/sqrt(2*var4^2*1)*exp(-((x-mean4)^2/(2*var4^2))) +
    
    1/sqrt(2*var5.1^2*1)*exp(-((y-mean5.1)^2/(2*var5.2^2))) * 
    1/sqrt(2*var5.2^2*1)*exp(-((x-mean5.2)^2/(2*var5.2^2)))}

z <- outer(x, y, norm)
z[is.na(z)] <- 1

persp(x, y, z, col = color, 
      xlab = "Trait 2", 
      ylab = "Trait 1", 
      zlab = "Fitness", 
      theta = 60,
      phi = 15,d = 10,
      border = border,
      shade = 0.6,
      zlim=c(0,.2))
