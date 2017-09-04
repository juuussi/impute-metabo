xc <- rnorm(25) # complete sample
xt <- xc[xc >= 0] # left truncated sample

models <- truncgof::ks.test(xt, "pnorm",list(0, 1), H = 0)


dplot(xc, "pnorm", list(0,1), vertical = TRUE)
# df of the left truncated sample
dplot(xt, "pnorm", list(0,1), H = 0, vertical = TRUE)



xt <-na.omit(miss_data[,1])
threshold <- min(na.omit((miss_data[,1])))

models <- truncgof::ks.test(xt, "pnorm",list(  mean(simulated_data),  sd(simulated_data)), H = threshold,  alternative ="two.sided")
