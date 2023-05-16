# Description  ------------------------------------------------------------
#### ### ### ## #### ### ### ## #### ### ### ## 
# Created by Marc-Olivier Beausoleil
# Why: 
  # example of drift (bottleneck or founder effect)
# Requires:
# NOTES: 
# Reference : 
#### ### ### ## #### ### ### ## #### ### ### ## 

par(mfrow = c(1,1), mar = c(0,0,5,0))
# set.seed(123456)

cols = c("blue", "yellow1", "pink", "chocolate2", "firebrick1", 
         "forestgreen", "orchid2", "orangered4", "skyblue1", "wheat2")

# How many marbles to select? 
(size.sel = 4)

# The greater the by, the less there are marbles 
by = 10
x = seq(1,100, by = by)
y = x
xy = expand.grid(x,y)
min.range = -10
max.range = 100
par(mfrow = c(1,2))
plot(xy[,1:2], 
     asp = 1, 
     main = "Original population",
     axes=FALSE, frame.plot=FALSE,
     xlim = range(c(min.range, xy[,"Var1"],max.range)),
     ylim = range(c(min.range, xy[,"Var2"],max.range)), 
     type = "n")


# random colours 
xy$col = sample(cols, size = nrow(xy), replace = TRUE)

for (i in 1:nrow(xy)) {
  plotrix::draw.circle(x = xy[i,1],y = xy[i,2], 
                       radius = log(by),
                       nv=1000,
                       border=NULL,
                       col=scales::alpha(xy[i,"col"], 1),
                       lty=1,lwd=.5)
}


sample.select = sample(1:nrow(xy), size = size.sel)
xy$gen = 1
# Select rows 
xy$select = ifelse((1:nrow(xy)) %in% sample.select, T, F)

select.col = xy[sort(sample.select),]

for (i in 1:nrow(select.col)) {
  plotrix::draw.circle(x = select.col[i,1],
                       y = select.col[i,2], 
                       radius = log(by),
                       nv=1000,
                       border="black",
                       col=select.col[i,"col"],
                       lty=1,lwd=5)
}

survived = xy[which(xy$select),]

tmp = 0
i = 1
while (tmp < nrow(xy)) {
  tmp.dat = survived[rep(1:nrow(survived), each = i), ]
  tmp = nrow(tmp.dat)
  i = i +1
}

new.gen = tmp.dat[sample(1:nrow(tmp.dat)),]
new.gen = new.gen[1:nrow(xy),]
new.gen[,1] <- xy$Var1
new.gen[,2] <- xy$Var2

plot(xy[,1:2], 
     asp = 1, 
     axes=FALSE, frame.plot=FALSE,
     main = "After bottleneck / founder effect",
     xlim = range(c(min.range, xy[,"Var1"],max.range)),
     ylim = range(c(min.range, xy[,"Var2"],max.range)), 
     type = "n")

for (i in 1:nrow(new.gen)) {
  plotrix::draw.circle(x = new.gen[i,1],y = new.gen[i,2], 
                       radius = log(by),
                       nv=1000,
                       border=NULL,
                       col=scales::alpha(new.gen[i,"col"], 1),
                       lty=1,lwd=.5)
}

# Original population
tab.col = table(xy$col)
# sum(tab.col)
prop.table(tab.col)
# as.data.frame(tab.col)

# Population After drift
tab.col.new.gen = table(new.gen$col)
# sum(tab.col.new.gen)
prop.table(tab.col.new.gen)
# as.data.frame(tab.col.new.gen)

par(mfrow = c(1,1), mar = c(4,4,4,4))
# plot(table(xy$col, new.gen$col))
track.allele = data.frame(cols, pop1 = NA, pop2 = NA)
for (k in names(tab.col)) {
  track.allele[track.allele$cols == k,"pop1"] <- tab.col[k]
  track.allele[track.allele$cols == k,"pop2"] <- ifelse(test = is.na(tab.col.new.gen[k]), "lost", tab.col.new.gen[k])
}

track.allele
