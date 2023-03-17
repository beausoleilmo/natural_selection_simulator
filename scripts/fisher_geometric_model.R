# Fisher's geometric model 
r = .01 # random mutation of phenotypic size r
n = 1 # number of characters 
z = 1 # Distance to the optimum 
standardized.mutational.size <- function(r,n,z) {
  (r*sqrt(n))/(2*z)
}

# Get the sandardized mutational size of various distances to the optimum 
x = standardized.mutational.size(r,n,seq(0,3, by = .0001))
plot(x = x, y = 1-pnorm(x), 
     xlim = c(0, 3), ylim = c(0, 1), 
     type = "l",
     main = "Fisher's geometric model of adaptation")
abline(h = c(0,.5), lty = 3)
