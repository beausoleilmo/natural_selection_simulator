# Information  ------------------------------------------------------------
# Requires imagemagick if you want to make the gif at the end 
# syss.path = system("Echo $PATH", intern = TRUE)
# Sys.getenv("PATH")
# Sys.setenv(PATH = paste(Sys.getenv("PATH"),syss.path,sep = ":"))
# Remove all variables 
rm(list = ls())
startTime = Sys.time()

# System variables --------------------------------------------------------
# dev.off()
# export plots in png? 
png.plot = TRUE
clean.PNG = TRUE
make.gif = FALSE 
ssleep = 0.0
# set.seed(12345)

source("scripts/0.initialize.R")

# Initial plot ------------------------------------------------------------
# Background 
back = "black"
buff = 0.1 # percentage
# Define the environment boundary 
maxrange = 50 # Range of the environment plane 
range.x = c(0-(maxrange*buff), maxrange + (maxrange*buff))
range.y = c(0-(maxrange*buff), maxrange + (maxrange*buff))

# Make initial plot 
if (png.plot) {
  png("outns/ns%04d.png", width = 5, height = 5, units = "in", res = 300) 
}
# Generate plot 
resetplot(title = "Initial", 
          xlim = range.x, 
          ylim = range.y, 
          col = back)

# Init food items ---------------------------------------------------------
# food.density = .01 # between 0 and 1 
# n.food = food.density*maxrange^2 # number of initial food items 
n.food = 50
col.food = "grey20" # Colour 
food.capacity = 30 # Maximum amount of food particles in the environment 
regen.food = 50 # Number of food items added if under the food.capacity of the environment
reset_food = TRUE 
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

# Radius of the detection range of each bacteria (Radius of detect range)
radius = 4

# Add detection range for each bacteria 
add.circles(pos.bac,radius = radius)

# Define traits -----------------------------------------------------------
# trait to evolve : Speed (movement in distance per frame = dist/frame)
# get random variable for species (get a value for speed for a bacterium species)
mean.speed.random = 7
mut.strength = .5
# Number of iterations (generations) total
nb.gen = 100

speed.max.sp = abs(rnorm(n.sp, 
                         mean = mean.speed.random, 
                         sd = 1))

# Get Speed phenotype for EACH bacterium of a particular species  
dist.bac = matrix(rnorm(n.bac.max, speed.max.sp, 1),
                  ncol = n.sp, 
                  byrow = TRUE)
# Make initial bacterium species database 
df.speed = data.frame(speed = dist.bac,
                      gen = 0,
                      org.id = 1:nrow(dist.bac))
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

# We need to keep track of 
# 1. the fitness of all traits (OK for speed)
# 2. The energy function (OK for speed)
# 3. The nb of food items, eaten and regenerated through the generations
# 4. the age of the org. and age at death 
# 5. Population size (net change, birth and death)
# 
# Model iterations --------------------------------------------------------
if (png.plot) {
  png("outns/ns%04d.png", width = 5, height = 5,units = "in",res = 300) 
}
food.eaten.in.gen = n.food
n.birth = n.bac.max
n.death = 0
# Energy function 
expected.energy = function(speed) -1/4* speed^2
for (no.gen in 1:nb.gen) {
  cat("\n\nGeneration #", no.gen, "---------------------\n")
  cat("There are", nrow(dist.bac),"bacteria\n")
  # Reset the fitness values in this generation 
  fit.val = matrix(0,
                   nrow = nrow(dist.bac),
                   ncol = n.sp)
  # Reset the energy values in this generation 
  energy.pt = 100
  energ.val = matrix(energy.pt,
                     nrow = nrow(dist.bac),
                     ncol = n.sp)
  run = 1 # reset run to 1
  # max.run = 4 # Number of frames-1 in which the bacteria can SEARCH for food items 
  sum.energy=sum(energ.val)
  
  if(reset_food & no.gen>1){
    # Generate position of food items 
    food.x = runif(n.food, 0, maxrange) 
    food.y = runif(n.food, 0, maxrange)
    food.df = data.frame(food.x, food.y)
    # Add food to plot  
    points(food.df, 
           col = "red", 
           pch = 15)
  }
  
  # Make the frames in which bacteria are searching 
  while (sum.energy > 0) {
    # Make initial plot 
    if (sum.energy > 0) {
      resetplot(title = paste(no.gen,run,sep = "_"))
    }
    
    cat("Set of iterations for time", run, "\n")
    cat("Number of food items: ", nrow(food.df), "\n")
    
    food.found = NULL
    for (i in 1:nrow(pos.bac)) { # refers to the bacterium
      
      if (nrow(food.df)==0) { # If there is no food, don't try to make the bacteria find it. juste go to the next iteration 
        next
      }
      
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
    if(!(reset_food)){
      if (food.capacity > nrow(food.df)) { # Check if food needs to be regenerated
        # regenerate new food 
        food.x.new = runif(regen.food, 0, maxrange)
        food.y.new = runif(regen.food, 0, maxrange)
        
        food.append = data.frame(food.x = food.x.new, food.y = food.y.new)
        food.df = rbind(food.df, food.append)
      } # End if food.capacity in the environment 
    } # End if regen.food == null 
    
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
    
    # Recalculate the energy 
    energ.val = energ.val + expected.energy(dist.bac)
    energ.val = ifelse(energ.val<= 0,0,energ.val)
    sum.energy = sum(energ.val)
    cat("total energy is:",round(sum.energy,2),"\n")
    
  } # End While run 
  
  # in each generation, find the bacteria that survived 
  survival.bac = fit.val >= 1
  # sum(!survival.bac) counts the number of bacteria that died  
  cat("---------------------\n")
  cat(sum(!survival.bac),"bacteria died\n")
  
  n.death = c(n.death, sum(!survival.bac))
  
  # Get phenotypes of surviving bacteria 
  bac.s.speed = dist.bac[survival.bac]
  
  # from the bacteria that survived, the one with a highest fitness will reproduce 
  reproduce.bac = fit.val >= 2
  nb.bac.repro = sum(reproduce.bac)
  cat(nb.bac.repro,"bacteria reproduced\n")
  n.birth = c(n.birth,nb.bac.repro)
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
  
  food.eaten.in.gen = c(food.eaten.in.gen,n.food-nrow(food.df))
  # If no bacteria reproduce, stop 
  if (length(pos.bac)==0) {
    stop("No bacterium survived!\n")
  }
  
  # Record the mean phenotype in each generation 
  rec.mean = c(rec.mean, mean(dist.bac))
  max.id = max(df.speed$org.id)+1
  
  nrow(dist.bac)
  current.gen=df.speed[df.speed$gen %in% (no.gen-1),"org.id"]
  
  # Record the individual phenotypes of each bacterium with the generation in which it is found 
  df.speed = rbind(df.speed,
                   data.frame(speed = dist.bac, 
                              gen = no.gen,
                              org.id = c(current.gen[survival.bac],
                                         seq(max.id,max.id+length(bac.r.speed)-1))))
} # End for loop of the NS_algorithm 

if (png.plot) {
  dev.off()
}
if (make.gif & nb.gen == no.gen) {
  fold.exists.gif = file.exists("gif")
  if (!fold.exists.gif) {
    dir.create("gif")
  }
  system("convert -delay 20 -loop 0 ~/Github_proj/natural_selection_simulator/outns/*.png ~/Github_proj/natural_selection_simulator/gif/ns.film.gif", intern = TRUE)
}

# Summary statistics ------------------------------------------------------
for (end in 1:1) {
  cat("\nNumber of gen:",no.gen+1,"---------------------\n")
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
  fold.exists.gif = file.exists("summary_stats")
  if (!fold.exists.gif) {
    dir.create("summary_stats")
  }
  pdf("summary_stats/popsize.pdf",width = 5, height = 5,) 
  
}

# Add population plot of the trait evolving 
ggplot(df.speed, aes(speed)) +
  geom_histogram() +
  geom_density(aes(y=1 * ..count..))+
  ggplot2::facet_wrap(~gen, ncol = 2) + 
  theme_bw()

par(mfrow = c(5,1))

# Make a plot of the change in mean speed through the generations 
plot(x = df.speed$gen+1, 
     y = df.speed$speed, 
     xlab = "Generations",
     ylab = "Mean speed",
     main = "Trait change over generations",
     xlim = c(1,length(rec.mean)),
     ylim = range(df.speed$speed), 
     pch =19,
     type = "n",
     cex =.9, 
     col = alpha("red",alpha = .1))

pol.dat = df.speed %>% 
  group_by(gen) %>% 
  summarise(sd = sd(speed),
            mean = mean(speed),
            min=min(speed),
            max=max(speed))

y.pol = c(pol.dat$mean+pol.dat$sd,rev(pol.dat$mean-pol.dat$sd))

polygon(x = c(1:(no.gen+1), rev(1:(no.gen+1))), 
        y = y.pol, 
        col = scales::alpha("grey50",.5), border = NA)

points(x = df.speed$gen+1, 
       y = df.speed$speed, 
       cex =.9, 
       pch =19,
       col = alpha("red",alpha = .1))
points(x = 1:length(rec.mean),
       y = rec.mean, 
       type = "l",lwd =2,
       pch =19, cex =.9, col = alpha("black",alpha = .9))
points(x = 1:length(rec.mean),
       y = pol.dat$min, 
       type = "l",lwd =2,
       pch =19, cex =.9, col = alpha("grey50",alpha = .6))
points(x = 1:length(rec.mean),
       y = pol.dat$max, 
       type = "l",lwd =2,
       pch =19, cex =.9, col = alpha("grey50",alpha = .6))
abline(h = rec.mean[1],lty = 3)
### 
plot(nb.bac~c(1:(no.gen+1)), 
     xlab = "Generations",
     ylab = "Number of bacteria",
     ylim =range(0,nb.bac),
     pch =19, type = "l", 
     main = "Population size")
plot(food.eaten.in.gen~c(1:(no.gen+1)), 
     xlab = "Generations",
     ylab = "Number of food items eaten",
     ylim =range(0,food.eaten.in.gen),
     pch =19, type = "l", 
     main = "Food eaten through time")
abline(h = n.food, lty=3)
popul.dyn.birth = data.frame(nb = n.birth,type ="birth",gen = 1:length(n.birth))
popul.dyn.death = data.frame(nb = -n.death,type ="death",gen = 1:length(n.death))
pop.dyn = rbind(popul.dyn.birth,popul.dyn.death)

bp.out = barplot(height = popul.dyn.birth$nb,
                 names.arg = popul.dyn.birth$gen,
                 col = "green",
                 main = "",
                 xlab = "Generations",
                 ylab = "Number of organisms",
                 ylim = range(pop.dyn$nb))
barplot(height = popul.dyn.death$nb,
        names.arg = popul.dyn.death$gen, 
        main = "Population change life table",
        col = "red",
        add = TRUE)
abline(h = 0)
points(x = bp.out,
       y = popul.dyn.birth$nb + popul.dyn.death$nb, 
       type = "l",
       lwd = 2,
       pch = 19)


# age summary

# make histogram of the maximum age for each individual
# hist((sort(as.vector(table(df.speed$org.id)))))
# Make a "capture history"
ch = table(df.speed$org.id,df.speed$gen)
# make a datafram ethat records the age of all individuals at all time 
age.df = NULL
for (age.ch in 1:ncol(ch)) {
  sum.up.to.gen = apply(ch[,1:age.ch, drop = FALSE],1,sum)
  age.df = cbind(age.df,sum.up.to.gen)
}
# When an organism is not found, switch it to 0 
age.df[which(ch==0)]<-0

plot(apply(age.df,2,mean),
     ylab = "Age",
     xlab = "Generations",
     ylim = range(age.df), pch =19)
lines(apply(age.df,2,mean))
lines(apply(age.df,2,min))
lines(apply(age.df,2,max))

par(mfrow = c(1,1))

if (png.plot) {
  dev.off()
}

# print elapsed time
endTime <- Sys.time() - startTime # calculate difference
print(endTime) 
