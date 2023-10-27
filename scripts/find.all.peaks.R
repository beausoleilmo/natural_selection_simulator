# Source: https://stackoverflow.com/questions/68618852/r-how-to-get-values-related-to-last-peak-of-density-plot

data <- c(11,35,35,35,6,75,49,74,82,49,
          75,8,74,37,73,7,47,47,72,48,
          46,9,73,49,73,51,50,9,73,47)

d <- density(data)
peaks <- NULL
for (i in 2:(length(d$y)-1)) {
  if (d$y[i-1] <= d$y[i] & d$y[i] >= d$y[i+1]) {
    peaks <- cbind(peaks, c(d$x[i], d$y[i]))
  }
}
peaks
# result:
# each column is a peak point
#          [,1]        [,2]        [,3]
#[1,] 9.8433702 46.92771112 70.90705938
#[2,] 0.0074287  0.01527099  0.01276465

# show the peak in a graph
plot(d)
points(peaks[1,], peaks[2,])

# third peak, if any
peaks[,3]
#[1] 70.90705938  0.01276465
# you can use this elsewhere you want for further analysis
