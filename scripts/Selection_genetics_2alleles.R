# Visualisation de la selection su HW -------------------------------------
# see Futuyma p. 109 
# positive selection favoring the spread of a beneficial mutation
# see also https://www.apsnet.org/edcenter/disimpactmngmnt/topc/PopGenetics/Pages/NaturalSelection.aspx 

# This is a more general solution to what is found in Futuyma 
selection.cal <- function(p=.5, s=0, h = .5, t = 2*.5 
) {
  q = 1-p
  p2 = p^2
  pq2 = p*q*2
  q2 = q^2
  s = s # percentage of increase in fitness (e.g., survival) for one allele 
  before.sel = c(p2, pq2, q2)
  # This part was added and might be wrong... 
  rel.fit.hr = (1+s)
  rel.fit.ht = (1+h)
  rel.fit.hd = (1+t)
  
  w_bar = 
    p2*rel.fit.hr + 
    pq2*rel.fit.ht + 
    q2*rel.fit.hd
  after.sel = c(p2 *rel.fit.hr/w_bar, 
                pq2*rel.fit.ht/w_bar, 
                q2 *rel.fit.hd/w_bar)
  q_prime = ((q2 *rel.fit.hd)+(p*q*rel.fit.ht))/w_bar
  df.ba = data.frame(Genotype = c("p^2","2pq","1^2"),
                     before.sel,
                     after.sel = round(after.sel, 2))
  return(list(p = p,
              q = q,
              s = s,
              h = h, 
              t = t,
              w_bar = w_bar,
              p_prime = 1-q_prime,
              q_prime = q_prime,
              df.change = df.ba))
}


draw.selection <- function(generations = 200, 
                           p=.5, 
                           s = 0, 
                           h = .05, 
                           t = .1, 
                           add = FALSE, 
                           # dom = "no",
                           show.p = TRUE) {
  generations = generations
  p=p
  show.var = ifelse(show.p,"Frequency of p","Frequency of q")
  evo.mat = matrix(NA,
                   ncol = 2,
                   nrow = generations)
  for (gen in 1:generations) {
    sel.out = selection.cal(p = p,
                            s = s, 
                            h = h, 
                            t = t
                            )
    evo.mat[gen,1] <- sel.out$p
    evo.mat[gen,2] <- gen
    p = sel.out$p_prime
  }
  if (show.p) {
    y.val = evo.mat[,1]
  } else {y.val = 1-evo.mat[,1]}
  
  if (!add) {
    plot(y.val~evo.mat[,2], type = "l",
         ylim = c(0,1),
         xlab = "Generations",
         ylab = show.var)    
  } else {
    lines(y.val ~ evo.mat[,2])
  }
  
}
draw.selection(p = 0.3, s = 0, h = .05, t = 0.1, show.p = T)
draw.selection(p = 0.6, s = 0, h = .05, t = 0.1, show.p = T, add = T)
draw.selection(p = 0.9, s = 0, h = .05, t = 0.1, show.p = T, add = T)


draw.selection(p = 0.3, s = 0, h = .05, t = 0.1, show.p = F)
draw.selection(p = 0.6, s = 0, h = .05, t = 0.1, show.p = F, add = T)
draw.selection(p = 0.9, s = 0, h = .05, t = 0.1, show.p = F, add = T)

draw.selection(p = 0.001, s = -.1, h = 0, t = -.2, show.p = F)
draw.selection(p = 0.35, s = -.1, h = 0, t = -.2, show.p = F, add = T)
draw.selection(p = 0.999, s = -.1, h = 0, t = -.2, show.p = F, add = T)

draw.selection(p = 0.001, s = .1, h = 0, t = .2, show.p = F)
draw.selection(p = 0.35, s = .1, h = 0, t = .2, show.p = F, add = T)
draw.selection(p = 0.999, s = .1, h = 0, t = .2, show.p = F, add = T)

# Positive 
# Evolution "towards" fixation 
for (i in seq(0.,1,length.out = 25)) {
  if (i==0) {
    add = FALSE
  } else {add = TRUE}
  draw.selection(p = i, s = 0, h = 0.05, t = .1, show.p = F, add = add)
}

# Overdominance 
# Stable equilibrium
for (i in seq(0.,1,length.out = 25)) {
  if (i==0) {
    add = FALSE
  } else {add = TRUE}
  draw.selection(p = i, s = -.1, h = 0, t = -.2, show.p = F, add = add)
}

# Underdominance 
# Unstable equilibrium
for (i in seq(0.,1,length.out = 25)) {
  if (i==0) {
    add = FALSE
  } else {add = TRUE}
  draw.selection(p = i, s = .1, h = 0, t = .2, show.p = F, add = add)
}




