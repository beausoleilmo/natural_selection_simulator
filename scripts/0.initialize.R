# Load libraries and functions  -------------------------------------------
# Remove all plots 
tryCatch(expr = dev.off(dev.list()["RStudioGD"]),error =function(e) print("Nothing to be closed"))

# Load libraries 
library(tidyverse)
# Load functions 
rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}
# Make plot 
resetplot <- function(title = NULL,
                      xlim=range.x, ylim=range.y, col = "black") {
  plot(NA, 
       xlim = xlim, ylab="y",xlab="x",
       ylim = ylim,
       asp =1, main = title)
  # Add background 
  rect(xleft = range.x[1], xright = range.x[2],
       ybottom = range.y[1], ytop = range.y[2],
       col = col)
}

add.circles <- function(position.org, radius, col = "grey30") {
  for (circ in 1:nrow(position.org)) {
    # prepare "circle data"
    center_x = position.org[circ,1]
    center_y = position.org[circ,2]
    theta = seq(0, 2 * pi, length = 200) # angles for drawing points around the circle
    
    # draw the circle
    lines(x = radius * cos(theta) + center_x, 
          y = radius * sin(theta) + center_y,
          col = col)
  }
}



# Creates dir 
fold.exists = file.exists("outns")
if (!fold.exists) {
  dir.create("outns")
}


if(clean.PNG){
  #Define the file name that will be deleted
  fn <- "outns/ns*.png"
  files.to.remove = file.path("outns",list.files(path = "outns/",pattern = ".png"))
  lapply(files.to.remove, file.remove)
}