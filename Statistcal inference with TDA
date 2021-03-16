seirmod2=function(t, y, parms){
  S=y[1]
  E=y[2]
  I=y[3]
  R=y[4]
  with(as.list(parms),{
    dS = mu * (N - S) - beta0 * (1+beta1 *
                                   cos(2 * pi * t)) * S * I / N
    dE = beta0 * (1 + beta1 * cos(2*pi * t)) *
      S * I / N - (mu + sigma) * E
    dI = sigma * E - (mu + gamma) * I
    dR = gamma * I - mu * R
    res=c(dS, dE, dI, dR)
    list(res)
  })
}  

require(deSolve)
times = seq(0, 5, by=1/120)
paras = c(mu = 1/50, N = 1, beta0 = 1000, beta1 = 0.2,
          sigma = 365/2, gamma = 365/5)
start = c(S=0.06, E=0, I=0.01, R = 0.939)
out = as.data.frame(ode(start, times, seirmod2, paras))
par(mfrow=c(1,2)) #Side-by-side plot
plot(times, out$I, ylab = "Infected",
     xlab = "Time", type = "l")
k = data.matrix(out, rownames.force = NA)


times = seq(0, 5, by=1/120)
paras = c(mu = 1/50, N = 1, beta0 = 1000, beta1 = 0.2,
          sigma = 365/20, gamma = 365/5)
start = c(S=0.06, E=0, I=0.001, R = 0.939)
out1 = as.data.frame(ode(start, times, seirmod2, paras))
plot(times, out1$I, ylab = "Infected",
     xlab = "Time", type = "l")
l = data.matrix(out1, rownames.force = NA)

# calculate homologies for both datasets
k.phom <- calculate_homology(k, dim = 1)
l.phom1 <- calculate_homology(l, dim = 1)

perm.test <- permutation_test(k,l, iterations = 100)
# display p-value for 0-cycles
print(perm.test[[1]]$pvalue)
# display p-value for 1-cycles
print(perm.test[[2]]$pvalue)

hist(perm.test[[1]]$permvals,
     xlab = "Wasserstein distance",
     ylab = "Counts",
     main = "Null distribution for 0-cycles",
     xlim = c(0, .4))
abline(v = perm.test[[1]]$wasserstein)

hist(perm.test[[2]]$permvals,
     xlab = "Wasserstein distance",
     ylab = "Counts",
     main = "Null distribution for 1-cycles",
     xlim = c(0, .4))
abline(v = perm.test[[2]]$wasserstein)




