x = c(seq(-15,15,by = .1))
fx <- function(x, betas = c(1,1), expo = 1, int = 0) {
  if (!(length(betas) == expo)) {
    stop("Number of coefficients should match the exponent: \n  -> Decrease the exponent or add betas")
  }
  x.par = paste0("x^", (expo):1)
  
  my.fun = paste(c(paste(betas, x.par, sep = " * "), int), collapse = " + ")
  fungen = as.expression(my.fun)
  D.sc <- D(parse(text=fungen), "x")
  
  polyout = eval(parse(text=fungen), envir=list(x=x)) # cbind(poly(x, expo), int)  %*% matrix(c(betas,1), ncol = 1)
  return(list(polyout = polyout,
              fungen = fungen,
              D.sc = D.sc))
}

nb.exp = 5
fx.out = fx(x, expo = nb.exp, betas = c(1/2,1/3,-1,2/5,  0), int = 0)
fx.out$fungen
fx.out$D.sc
plot(x = x, y = fx.out$polyout, xlim = c(-10, 10), ylim = c(-10, 10), type = "l")
# plot(x = x, y = x^2, xlim = c(-10, 10), ylim = c(-10, 10))
points(x = x, eval(fx.out$D.sc, envir=list(x=x)), 
       type = "l", 
       col = "red", lwd = 3, lty = 1)
abline(h = 0, v = 0, lty = 2)




# Example in Graham Coop's pop gen book 
x = c(seq(-3,3,by = .01))

nb.exp = 4
fx.out = fx(x, expo = nb.exp, betas = c(-1/3,-5/6,0,1), int = 0)
fx.out$fungen
fx.out$D.sc
plot(x = x, y = fx.out$polyout, xlim = c(-2, 2), ylim = c(-1, 1), type = "l")
# plot(x = x, y = x^2, xlim = c(-10, 10), ylim = c(-10, 10))
points(x = x, eval(fx.out$D.sc, envir=list(x=x)), 
       type = "l", 
       col = "red", lwd = 3, lty = 1)
abline(h = 0, v = 0, lty = 2)
f = function(x) eval(fx.out$D.sc)
abline(v = uniroot(f = f,interval = c(-2,2))$root, lty = 3)
abline(v = uniroot(f = f,interval = c(-2,-1))$root, lty = 3)
abline(v = uniroot(f = f,interval = c(-1,0))$root, lty = 3)




samp=c(12, 10, 8, 13)
100-sum(samp)

p_a.and.b = samp[2]/100
p_a = samp[1]/100
p_b = samp[3]/100

p.jab = p_a.and.b/p_b
p.jba = p_a.and.b/p_a

p.jab*p_b
p.jba*p_a
# Bayes' Rule
# p.jba = p.jab*p_b/p_a

n = 4
k = 2
choose(n,k)
factorial(n)/factorial(n-k)