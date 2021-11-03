# set.seed(12345)
par(mfrow = c(2,2))
x = rnorm(100, 15, 1)
y = x + rnorm(length(x), 0, 1)

# Each dot represent a family (mean from mom and dad (x), and mean of all offsprings traits (regardless of the number of offspings). 
# BUT the evolution is shown in the MEAN of the cloud of points
plot(y~x, pch =19, main = "Without selection",
     xlab = "Mean parents (male and female)",
     ylab = "Mean offspings")
lm.out=lm(y~x)
abline(lm.out)
points(mean(x),mean(y), 
       col = scales::alpha("black",.5), 
       pch =19, cex = 4)
h.s.event = 15
x.hs = x[x>h.s.event]
y.hs = y[x>h.s.event]
plot(y~x, #pch =19, 
     col = "red",
     main = "With selection",
     type = "n",
     xlab = "Mean parents (male and female)",
     ylab = "Mean offspings")
points(y.hs~x.hs, pch =19, col = "darkgreen", bg = "darkgreen")
points(y[x<=h.s.event]~x[x<=h.s.event], pch =21, 
       col = scales::alpha("blue",.5), bg = scales::alpha("blue",.5))
abline(lm.out)
abline(v = h.s.event, lty = 3)

points(mean(x.hs),mean(y.hs), 
       col = scales::alpha("darkgreen",.5), 
       pch =19,cex = 4)

points(mean(x),mean(y), 
       col = scales::alpha("black",.5), 
       pch =19,cex = 2)

plot(density(x.hs),col ="darkgreen",
     xlim = c(12,18),
     main = "Distribution of parental phenotypes")
lines(density(x),col = "black")
abline(v = h.s.event, col = "red")

breaks = seq(floor(min(x)),ceiling(max(x)),1)
x.hist = hist(x, breaks=breaks, col = scales::alpha("blue",.5))
x.hs.hist = hist(x.hs, add=TRUE, breaks=breaks)
plot(x.hs.hist, col = scales::alpha("darkgreen",.5), add = TRUE)
abline(v = h.s.event, col = "red")


