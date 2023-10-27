# Source (translated from python)
# https://gitlab.com/ole.tange/just-for-fun/-/blob/master/python/diceevolution.py?ref_type=heads

# How does evolution through natural selection work?
# What is the chance of getting 1000 sixes with 1000 dice?
# If you only throw once: 1/6^1000 =
# 7.059 × 10^-779 (see https://www.wolframalpha.com/)
# But if you keep the sixes after each throw, then after just 40 throws
# the chance is 50% that you have 1000 sixes.

# Number of replicates
tests <- 10000

# Function that will throw die
throw.die <- function(dice = 1000) {
  sixes <- 0              # Start with 0 sixes 
  throws <- 0             # Start without throwing
  while (sixes != dice) { # If nb of sixes is NOT equal to nb of die, continue. 
    throws <- throws + 1  # Add one to count the number of throws (a trhow is tossing ALL non-six dice at the same time)
    # Get a vector will all possible faces of all die that were not 0 
    throw.all = sample(x = 1:6,           # Get faces of a dice 
                       size = dice-sixes, # Nb of die minus number of sixes represent the NON-6 die that need to be tossed 
                       replace = TRUE)    # Replace as there will be more die than the number of faces 
    nb.six = length(which(throw.all==6))  # Count the number of 6s in the data 
    sixes = sixes + nb.six                # Add the nubmer of sixes to the previous counts from the new throw
  }
  return(throws)           # Return number of throws 
}

# Replicate the experiement multiple time to get the 'long-term' effect 
nb.throws.all = replicate(n = tests, 
                          expr = throw.die(dice = 1000))
# Make histogram 
hist(nb.throws.all);abline(v = median(nb.throws.all))
print(median(nb.throws.all))
# 1 => 4
# 10 => 16
# 100 => 28
# 1000 => 40-41
# 10000 => 53
# 100000 => 66
# 1000000 => 78
rec.med = c()
max.pwr = 4
for (i in 10^(seq(0,max.pwr, by = 1))) {
  print(i)
  nb.throws.all = replicate(n = tests, 
                            expr = throw.die(dice = i))
  rec.med = c(rec.med, median(nb.throws.all))
}
dat = data.frame(x = 10^(seq(0,max.pwr, by = 1)),
           y = rec.med)
dat = rbind(dat, data.frame(x = c(10^5, 10^6),
                      y = c(66,78)))
plot(dat, type = "p")
