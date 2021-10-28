# Information  ------------------------------------------------------------
# Requires imagemagick if you want to make the gif at the end 
# syss.path = system("Echo $PATH", intern = TRUE)
# Sys.getenv("PATH")
# Sys.setenv(PATH = paste(Sys.getenv("PATH"),syss.path,sep = ":"))

# Load libraries and functions  -------------------------------------------
# Remove all variables 
rm(list = ls())
# Remove all plots 
dev.off(dev.list()["RStudioGD"])
# set.seed(12345)
# Creates dir 
fold.exists = file.exists("~/Desktop/outns")
if (!fold.exists) {
  dir.create("~/Desktop/outns")
}
# export plots in png? 
png.plot = TRUE 
# dev.off()
clean.PNG = TRUE
if(clean.PNG){
  #Define the file name that will be deleted
  fn <- "~/Desktop/outns/ns*.png"
  files.to.remove = file.path("~/Desktop/outns",list.files(path = "~/Desktop/outns/",pattern = ".png"))
  lapply(files.to.remove, file.remove)
}

make.gif = TRUE 


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


# System variables --------------------------------------------------------
ssleep = 0.0

# Initial plot ------------------------------------------------------------
# Background 
back = "black"
buff = 0.1 # percentage
# Define the environment boundary 
maxrange = 50
range.x = c(0-(maxrange*buff), maxrange + (maxrange*buff))
range.y = c(0-(maxrange*buff), maxrange + (maxrange*buff))

# Make initial plot 
if (png.plot) {
  png("~/Desktop/outns/ns%04d.png",width = 5, height = 5,units = "in",res = 300) 
}
resetplot(title = "Initial", 
          xlim = range.x, 
          ylim = range.y, 
          col = back)

# Init food items ---------------------------------------------------------

food.density = .04 # between 0 and 1 
n.food = food.density*maxrange^2 # number of initial food items 
col.food = "grey20" # Colour 
food.capacity = 30 # Maximum amount of food particles in the environment 
regen.food = 50 # Number of food items added if under the food.capacity of the environment

# Generate position of food items 
food.x = runif(n.food, 0, maxrange) 
food.y = runif(n.food, 0, maxrange)
food.df = data.frame(food.x, food.y)
initial.food = food.df

# Add food to plot  
points(food.df, 
       col = col.food, 
       pch = 15)


# Init bacteria individuals -----------------------------------------------
# Generate bacteria (position = (x,y) and colour = species)
n.sp = 1
n.bac.max = 50 # this has to be even to simplify programming if the number of sp is greater than 1 
nb.bac = n.bac.max
# If more than one species, you can select the colours here 
col.bac = c("yellow","blue","red",
            "green","purple","orange","grey50","cyan")

# Generate coordinates for bacteria in the environment
pos.bac = matrix(runif(n.bac.max*n.sp*2, 0, maxrange),
                 ncol = n.sp*2)
initial.pos.bac = pos.bac

# Radius of the detection range of each bacteria 
radius = 4


# Add detection range for each bacteria 
add.circles(pos.bac,radius = radius)

# Define traits -----------------------------------------------------------
# trait to evolve : Speed (movement in distance per frame = dist/frame)
# get random variable for species (get a value for speed for a bacterium species)
mean.speed.random = 7
mut.strength = .5
# Number of iterations (generations) total
it.max = 100

speed.max.sp = abs(rnorm(n.sp, 
                         mean = mean.speed.random, 
                         sd = 1))

# Get Speed phenotype for EACH bacterium of a particular species  
dist.bac = matrix(rnorm(n.bac.max, speed.max.sp, 1),
                  ncol = n.sp, 
                  byrow = TRUE)
# Make initial bacterium species database 
df.speed = data.frame(speed = dist.bac,
                      gen = 0)
# Get average initial speed  
avg.initial = mean(dist.bac)

# Add bacteria to the environment 
points(pos.bac[,1], 
       pos.bac[,2], 
       col = col.bac[n.sp],
       # col = colorRampPalette(c('yellow', 'red'))(length(dist.bac))[rank(dist.bac)], 
       pch = 15)



# Record mean -------------------------------------------------------------
# Record the mean of all generations 
rec.mean = avg.initial
if (png.plot) {
  dev.off()
}

# Model iterations --------------------------------------------------------
if (png.plot) {
  png("~/Desktop/outns/ns%04d.png",width = 5, height = 5,units = "in",res = 300) 
}
for (ittt in 1:it.max) {
  cat("\n\nGeneration #",ittt,"---------------------\n")
  cat("There are", nrow(dist.bac),"bacteria\n")
  # Reset the fitness values in this generation 
  fit.val = matrix(0,
                   nrow = nrow(dist.bac),
                   ncol = n.sp)
  run = 1 # reset run to 1 
  max.run = 4 # Number of frames-1 in which the bacteria can SEARCH for food items 
  
  # Make the frames in which bacteria are searching 
  while (run != max.run) {
    # Make initial plot 
    if (run != max.run) {
      resetplot(title = paste(ittt,run,sep = "_"))
    }
    
    cat("Set of iterations for time", run, "\n")
    cat("Number of food items: ", nrow(food.df), "\n")
    
    food.found = NULL
    for (i in 1:nrow(pos.bac)) { # refers to the bacterium
      for (j in 1:nrow(food.df)) { # refers to the food items 
        crit = (food.df[j,1]-pos.bac[i,1])^2 + (food.df[j,2] - pos.bac[i,2])^2 < radius^2
        dat.app = data.frame(i,j,crit,
                             food.x = food.df[j,1],
                             food.y = food.df[j,2],
                             posb.x = pos.bac[i,1], 
                             posb.y = pos.bac[i,2])
        if (crit) {
          # Append the food items that are found
          food.found = rbind(food.found,dat.app)
          # Remove the food items found more than once (i.e., only one bacterium can have a food item. The bacteria can't share food)
          food.found = food.found[!duplicated(food.found$j), ]
        } # End if "crit" 
      } # End for j
    } # End for i 
    
    if (!is.null(food.found)) { # If the food.found is not null, execute 
      j.delete = food.found %>%
        group_by(i) %>%
        sample_n(1) %>% # get a random piece of food within the detection range
        ungroup()%>%
        select(j) %>%
        unlist() %>%  
        as.vector()
      
      i.fit = unique(food.found$i) # Eat only 1 piece of food in the range even if there are multiple food items in reach 
      fit.val[i.fit,] <- fit.val[i.fit,] + 1 # 1 food item = 1 unit of fitness 
      
      food.df = food.df[-j.delete,] # Remove the food from the environment 
    } # End if food.found 
    
    # numb.eaten = length(j.delete)
    if (food.capacity > nrow(food.df)) { # Check if food needs to be regenerated
      # regenerate new food 
      food.x.new = runif(regen.food, 0, maxrange)
      food.y.new = runif(regen.food, 0, maxrange)
      
      food.append = data.frame(food.x = food.x.new, food.y = food.y.new)
      food.df = rbind(food.df, food.append)
    } # End if food.capacity in the environment 
    
    # Add all the food items PLUS the ones that were generated if needed 
    points(food.df, 
           col = "red", 
           pch =15)
    
    # Make bacteria move 
    # Generate random angle in degrees 
    angles.bac.rad = matrix(runif(nrow(dist.bac), 0, 360),
                            ncol = n.sp)
    # Convert angle degrees to rad 
    angles.bac.deg = deg2rad(angles.bac.rad)
    
    # Energy for distance 
    # energie.cost = apply(dist.bac,2,function(x) .5 * x + rnorm(nrow(dist.bac)))
    
    col.mat.bac = ncol(pos.bac)
    even_indexes <- seq(2,col.mat.bac,2)
    odd_indexes <- seq(1,col.mat.bac-1,2)
    
    pos.bac2 = pos.bac
    # Add a vector with some angle and the "speed" phenotype which is how much a bacterium move in the next frame 
    pos.bac2[, even_indexes] = pos.bac2[, even_indexes] + cos(angles.bac.rad) * dist.bac
    pos.bac2[, odd_indexes] = pos.bac2[, odd_indexes] + sin(angles.bac.rad) * dist.bac
    
    # Make sure the bacteria are not OUTSIDE the plotting area 
    pos.bac2[, even_indexes][which(pos.bac2[, even_indexes] > maxrange)] <- maxrange
    pos.bac2[, even_indexes][which(pos.bac2[, even_indexes] < 0)] <- 0
    pos.bac2[, odd_indexes][which(pos.bac2[, odd_indexes] > maxrange)] <- maxrange
    pos.bac2[, odd_indexes][which(pos.bac2[, odd_indexes] < 0)] <- 0
    
    # Add the bacteria position to the plot 
    points(pos.bac2[,1], 
           pos.bac2[,2], 
           col = col.bac[n.sp],
           pch =15)
    # Add the detection range 
    add.circles(pos.bac2, radius = radius)
    
    pos.bac = pos.bac2
    
    # Make time for the user to see the result 
    if (ssleep!=0) {
      Sys.sleep(ssleep)  
    }
    ### Increment the time to search 
    run = run + 1
  } # End While run 
  
  # in each generation, find the bacteria that survived 
  survival.bac = fit.val >= 1
  # sum(!survival.bac) counts the number of bacteria that died  
  cat("---------------------\n")
  cat(sum(!survival.bac),"bacteria died\n")
  
  # Get phenotypes of surviving bacteria 
  bac.s.speed = dist.bac[survival.bac]
  
  # from the bacteria that survived, the one with a highest fitness will reproduce 
  reproduce.bac = fit.val >= 2
  nb.bac.repro = sum(reproduce.bac)
  cat(nb.bac.repro,"bacteria reproduced\n")
  
  # position of the bacteria that survived 
  pos.suv = pos.bac[survival.bac,]
  # position of the bacteria that reproduce
  pos.rep = pos.bac[reproduce.bac,]
  # 
  pos.bac = rbind(pos.suv, # position of the bacteria that survived 
                  pos.rep) # position of the duplicated bacteria (daugther cell will be on top of the maternal bacteria)
  
  # Get the phenotypes of the bacteria that reproduce 
  bac.r.speed = dist.bac[reproduce.bac]
  # Generate random fluctuations in the phenotype of the bacteria that reproduce 
  mutation = rnorm(nb.bac.repro, mean = 0, sd = mut.strength)
  
  # Get new phenotype of the bacteria that reproduced 
  dist.bac = matrix(c(bac.s.speed, bac.r.speed + mutation),
                    ncol = n.sp, 
                    byrow = TRUE)
  nb.bac = c(nb.bac,nrow(dist.bac))
  # If no bacteria reproduce, stop 
  if (length(pos.bac)==0) {
    stop("No bacterium survived!\n")
  }
  
  # Record the mean phenotype in each generation 
  rec.mean = c(rec.mean, mean(dist.bac))
  # Record the individual phenotypes of each bacterium with the generation in which it is found 
  df.speed = rbind(df.speed,
                   data.frame(speed = dist.bac, 
                              gen = ittt))
}

if (png.plot) {
  dev.off()
}
if (make.gif & it.max == ittt) {
  system("convert -delay 20 -loop 0 ~/Desktop/outns/*.png ~/Desktop/outns/ns.film.gif", intern = TRUE)
}

# Summary statistics ------------------------------------------------------
for (end in 1:1) {
  cat("\nNumber of gen:",ittt+1,"---------------------\n")
  cat("#initinal bact.:", n.bac.max,"\n")
  cat("  # final bact.:", nrow(dist.bac),"\n")
  avg.final = mean(dist.bac)
  cat("Initial average:", round(avg.initial,2),"\n")
  cat("  Final average:", round(avg.final,2),"\n")
  cat("  #initial food:", n.food,"\n")
  cat("   # final food:", nrow(food.df),"\n")
}

dim(dist.bac)

if (png.plot) {
  pdf("~/Desktop/outns/popsize.pdf",width = 5, height = 5,) 
  
}
# Make a plot of the change in mean speed through the generations 
plot(x = df.speed$gen, 
     y = df.speed$speed, 
     xlab = "Generations",
     ylab = "Mean speed",
     xlim = c(1,length(rec.mean)),
     ylim = range(df.speed$speed),pch =19,
     type = "p",cex =.9, col = alpha("red",alpha = .1))

points(x = 1:length(rec.mean),
       y = rec.mean, type = "l",lwd =2,
       pch =19, cex =.9, col = alpha("black",alpha = .9))

# Add population plot of the trait evolving 
ggplot(df.speed, aes(speed)) +
geom_histogram() +
  geom_density(aes(y=1 * ..count..))+
  ggplot2::facet_wrap(~gen, ncol = 2) + 
  theme_bw()

plot(nb.bac~c(1:(ittt+1)), pch =19, type = "l")

if (png.plot) {
  dev.off()
}
