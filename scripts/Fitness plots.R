# Description  ------------------------------------------------------------
#### ### ### ## #### ### ### ## #### ### ### ## 
# Fake fitness landscapes 
# Fitness test theoretical fitness landscapes 
# Created by Marc-Olivier Beausoleil
# 20??
# Why: 
  # Simple schematic representations of fitness landscapes in 2D 
  # Gets the basics of directional, stabilizing and disruptive selection
# Requires:
# NOTES: 
# Reference : 
#### ### ### ## #### ### ### ## #### ### ### ## 

par(mfrow=c(1,3))
q.fun = function(x, xlab = xlab,
                 ylab = "Fitness",
                 lwd = 3,cex.text = 1,
                 col.main = "black", 
                 col.line = "red", 
                 col.box = "black", 
                 col.text = "black") {
  y = -x^2
  plot(y~x, type = "l", axes=FALSE, xaxs="i",yaxs="i",
       frame.plot=FALSE, 
       xlab="", ylab = "",lwd = lwd, col = col.main,
       ylim = c(min(y),max(y)+5000))
  mtext(side=2,text=ylab,line = 1.1, col = col.text, cex = cex.text)
  mtext(side=1,text=xlab,line = 1.1, col = col.text, cex = cex.text)
  axis(side = 1, labels = FALSE, tck = 0.000000001, lwd = 4, col = col.box)
  axis(side = 2, labels = FALSE, tck = 0.000000001, lwd = 4, col = col.box)
  
  y2 = -3*x^2+3500
  points(y2~x, type = "l", lty = 2, lwd = lwd,
         col = col.line)
}
# col.main = "white"
col.main = "black"
xlab = ""
q.fun(-100:100,xlab = xlab,lwd = 5,cex.text = 2,
      ylab = "Fitness", col.line = NA,
      col.main = col.main,col.box = col.main, col.text = col.main)


l.fun = function(x, xlab = "Phenotype",ylab = "Fitness",lwd = 3,cex.text = 1,
                 col.main = "black", col.line = "red", col.box = "black", col.text = "black") {
  y = -x^2
  plot(y ~ x, type = "l", 
       axes=FALSE, xaxs="i",yaxs="i",
       frame.plot=FALSE, 
       xlab="", ylab = "",
       lwd = lwd, col = col.main,
       ylim = c(min(y),max(y)+5000))
  mtext(side=2,text=ylab,line = 1.1, col = col.text, cex = cex.text)
  mtext(side=1,text=xlab,line = 1.1, col = col.text, cex = cex.text)
  axis(side = 1, labels = FALSE, tck = 0.000000001, lwd = 4, col = col.box)
  axis(side = 2, labels = FALSE, tck = 0.000000001, lwd = 4, col = col.box)
  
  y2 = -(x-25)^2
  points(y2~x, type = "l", lty = 2, lwd = lwd,
         col = col.line)
}
par(bg = NA)
l.fun(-100:100,
      xlab = xlab,lwd = 5,cex.text = 2,
      ylab = "Fitness", 
      col.line = "red",
      col.main = "white",
      col.box = "white", 
      col.text = "white")



par(mfrow=c(1,3),bg = NA)
d.fun = function(x, xlab = xlab,
                 ylab = "Fitness",
                 lwd = 3,
                 cex.text = 1,
                 col.main = "black", 
                 col.line = "red", 
                 col.box = "black", 
                 col.text = "black") {
  y = -x^2
  plot(y~x, type = "l", 
       axes=FALSE,
       frame.plot=FALSE, 
       xlab="", 
       ylab = "",
       xaxs="i",
       yaxs="i",
       ylim = c(min(y), max(y)+5000),
       lwd = lwd, col = col.main)
  mtext(side=2,text = ylab,line = 1.1, 
        col = col.text, cex = cex.text)
  mtext(side=1,text = xlab,line = 1.1, 
        col = col.text, cex = cex.text)
  # box(lty = 1, col = col.box, lwd = lwd)
  axis(side = 1, labels = FALSE, tck = 0.000000001, lwd = 4, col = col.box)
  axis(side = 2, labels = FALSE, tck = 0.000000001, lwd = 4, col = col.box)
  
  y2 =-x*-x*(-x/35-2)*(x/35-2)-4250
  points(y2~x, type = "l", 
         lty = 2, lwd = lwd,
         col = col.line)
  return(y2)
}

pt = d.fun(x = -100:100,
           xlab = xlab,
           lwd = 5, 
           cex.text = 2,
           ylab = "Fitness", 
           col.line = "red",
           col.main = col.main,
           col.box = col.main, 
           col.text = col.main)

# x = -100:100
# points(pt~x, type = "l", lty = 2,col = "red")


# dev.off()
par(bg = NA)
d.fun = function(x, xlab = "Phenotype",ylab = "Fitness",lwd = 3,cex.text = 1,
                 col.main = "black", col.line = "red", col.box = "black", col.text = "black") {
  y = -3*x*-3*x*(-x/40-2)*(x/40-2)-50000
  plot(y~x, type = "l", 
       axes=FALSE,
       frame.plot=FALSE, xlab="", ylab = "",xaxs="i",yaxs="i",
       ylim = c(min(y),max(y)+99000),lwd = lwd, col = col.main)
  mtext(side=2,text = ylab,line = 1.1, col = col.text, cex = cex.text)
  mtext(side=1,text = xlab,line = 1.1, col = col.text, cex = cex.text)
  # box(lty = 1, col = col.box, lwd = lwd)
  axis(side = 1, labels = FALSE, tck = 0.000000001, lwd = 4, col = col.box)
  axis(side = 2, labels = FALSE, tck = 0.000000001, lwd = 4, col = col.box)
  
  y2 =-5*x*-5*x*(-x/36-2)*(x/36-2)-99000
  points(y2~x, type = "l", lty = 2, lwd = lwd,
         col = col.line)
  return(y2)
}

pt = d.fun(x = -100:100,
           xlab = xlab,lwd = 5,cex.text = 2,
           ylab = "", 
           col.line = "red",
           col.main = col.main,
           col.box = col.main, col.text = col.main)

# x = -100:100
# points(pt~x, type = "l", lty = 2,col = "red")

# dev.off()
par(bg = NA)
d.fun = function(x, xlab = "Phenotype",ylab = "Fitness",lwd = 3,cex.text = 1,
                 col.main = "black", col.line = "red", col.box = "black", col.text = "black") {
  y = -10*x*-10*x*(-x/50-2)*(x/50-2)-50000
  plot(y~x, type = "l", 
       axes=FALSE,
       frame.plot=FALSE, xlab="", ylab = "",xaxs="i",yaxs="i",
       ylim = c(min(y),max(y)+99000),lwd = lwd, col = col.main)
  mtext(side=2,text = ylab,line = 1.1, col = col.text, cex = cex.text)
  mtext(side=1,text = xlab,line = 1.1, col = col.text, cex = cex.text)
  # box(lty = 1, col = col.box, lwd = lwd)
  axis(side = 1, labels = FALSE, tck = 0.000000001, lwd = 4, col = col.box)
  axis(side = 2, labels = FALSE, tck = 0.000000001, lwd = 4, col = col.box)
  
  y2 =-8*x*-8*x*(-x/45-2)*(x/45-2)+300000
  points(y2~x, type = "l", lty = 2, lwd = lwd,
         col = col.line)
  return(y2)
}

pt = d.fun(x = -100:100,
           xlab = xlab,lwd = 5,cex.text = 2,
           ylab = "", col.line = "red",col.main = col.main,col.box = col.main, col.text = col.main)

# x = -100:100
# points(pt~x, type = "l", lty = 2,col = "red")

