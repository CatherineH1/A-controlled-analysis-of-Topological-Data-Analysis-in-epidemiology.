install.packages("epimdr")
install.packages("deSolve")
install.packages("rootSolve")
install.packages("phaseR")
install.packages("shiny")

require(deSolve)
#Sir model
sirmod = function(t, y, parms) {
  # Pull state variables from y vector
  S = y[1]
  I = y[2]
  R = y[3]
  # Pull parameter values from parms vector
  beta = parms["beta"]
  mu = parms["mu"]
  gamma = parms["gamma"]
  N = parms["N"]
  # Define equations
  dS = mu * (N - S) - beta * S * I/N
  dI = beta * S * I/N - (mu + gamma) * I
  dR = gamma * I - mu * R
  res = c(dS, dI, dR)
  # Return list of gradients
  list(res)
}
times = seq(0, 26, by = 1/10)
parms = c(mu = 0, N = 1, beta = 2, gamma = 1/2)
start = c(S = 0.999, I = 0.001, R = 0)
out=ode(y=start, times=times, func=sirmod, parms=
          parms)
out=as.data.frame(out)
head(round(out, 3))
plot(x=out$time, y=out$S, ylab="Fraction", xlab=
       "Time", type="l")
lines(x=out$time, y=out$I, col="red")
lines(x=out$time, y=out$R, col="green")

#Calculate R0
R0=parms["beta"]/(parms["gamma"]+parms["mu"])
#Adjust margins to accommodate a second right axis
par(mar = c(5,5,2,5))
#Plot state variables
plot(x=out$time, y=out$S, ylab="Fraction", xlab="Time",
     type="l")
lines(x=out$time, y=out$I, col="red")
lines(x=out$time, y=out$R, col="green")
#Add vertical line at turnover point
xx=out$time[which.max(out$I)]
lines(c(xx,xx), c(1/R0,max(out$I)), lty=3)
#prepare to superimpose 2nd plot
par(new=TRUE)
#plot effective reproductive ratio (w/o axes)
plot(x=out$time, y=R0*out$S, type="l", lty=2, lwd=2,
     col="black", axes=FALSE, xlab=NA, ylab=NA,
     ylim=c(-.5, 4.5))
lines(c(xx, 26), c(1,1), lty=3)
#Add right-hand axis for RE
axis(side = 4)
mtext(side = 4, line = 4, expression(R[E]))
#Add legend
legend("right", legend=c("S", "I", "R",
                         expression(R[E])), lty=c(1,1,1, 2),
       col=c("black", "red", "green", "black"))
require(rootSolve)
equil=runsteady(y=c(S=1-1E-5, I=1E-5, R=0),
                times=c(0,1E5), func=sirmod, parms=parms)
round(equil$y, 3)



#Candidate values for R0 and beta (not that relevent)
R0 = seq(0.1, 5, length=50) 
betas= R0 * 1/2
#Vector of NAs to be filled with numbers
f = rep(NA, 50)
#Loop over i from 1, 2, ..., 50
for(i in seq(from=1, to=50, by=1)){
  equil=runsteady(y=c(S=1-1E-5, I=1E-5,
                      R=0), times=c(0,1E5), func=sirmod,
                  parms=c(mu=0, N=1, beta=betas[i], gamma=1/2))
  f[i]=equil$y["R"]
}
plot(R0, f, type="l", xlab=expression(R[0]))
curve(1-exp(-x), from=1, to=5, add=TRUE, col="red")

#2.4 opn epidemic
times = seq(0, 52*50, by=.1)
parms = c(mu = 1/(50*52), N = 1, beta = 2,
          gamma = 1/2)
start = c(S=0.19, I=0.01, R = 0.8)
out = as.data.frame(ode(y=start, times=times,
                        func=sirmod, parms=parms))
par(mfrow=c(1,2)) #Make room for side-by-side plots
plot(times, out$I, ylab="Fraction", xlab="Time",
     type="l")
plot(out$S, out$I, type="l", xlab="Susceptible",
     ylab="Infected")




#More Realistic Infectious Periods (plotting issue)
chainSIR=function(t, logx, params){
  x=exp(logx)
  u=params["u"]
  S=x[1]
  I=x[2:(u+1)]
  R=x[u+2]
  with(as.list(params),{
    dS = mu * (N - S) - sum(beta * S * I) / N
    dI = rep(0, u)
    dI[1] = sum(beta * S * I) / N - (mu + u*gamma) * I[1]
    if(u>1){
      for(i in 2:u){
        dI[i]= u*gamma * I[i-1] - (mu+u*gamma)* I[i]
      }
    }
    dR = u*gamma * I[u] - mu * R
    res=c(dS/S, dI/I, dR/R)
    list(res)
  })
}
times = seq(0, 10, by=1/52)
paras2 = c(mu = 1/75, N = 1, beta = 625,
           gamma = 365/14, u=1)
xstart2 = log(c(S=.06, I=c(0.001, rep(0.0001,
                                      paras2["u"]-1)), R = 0.0001))
out = as.data.frame(ode(xstart2, times, chainSIR,
                        paras2))
plot(times, exp(out[,3]), ylab="Infected", xlab=
       "Time", ylim=c(0, 0.01), type=’l’)
paras2["u"] =2
xstart2 = log(c(S=.06, I=c(0.001, rep(0.0001/
                                        paras2["u"], paras2["u"]-1)), R = 0.0001))
out2 = as.data.frame(ode(xstart2, times, chainSIR,
                         paras2))
lines(times, apply(exp(out2[,-c(1:2,length(out2))]),
                   1 ,sum), col=’blue’)
paras2["u"] =73
xstart2 = log(c(S=.06, I=c(0.001, rep(0.0001/
                                        paras2["u"], paras2["u"]-1)), R = 0.0001))
out3 = as.data.frame(ode(xstart2, times, chainSIR,
                         paras2))
lines(times, apply(exp(out3[,-c(1:2,length(out3))]),
                   1, sum), col=’red’, lwd=2, lty=2)
paras2["u"] =500
xstart2 = log(c(S=.06, I=c(0.001, rep(0.0001/
                                        paras2["u"], paras2["u"]-1)), R = 0.0001))
out4 = as.data.frame(ode(xstart2, times, chainSIR,
                         paras2))
lines(times, apply(exp(out4[,-c(1:2,length(out4))]),
                   1,sum, na.rm=TRUE), col=’green’)
legend("topright", legend=c("SIR", "u=2", "u=500",
                            "u=73 (H-S)"), lty=c(1,1,1,2), lwd=c(1,1,1, 2),
       col=c("black", "blue", "green", "red"))


#3.5 Stochastic Simulation (Plotting issue)
sim.cb=function(S0, beta, I0){
  I=I0
  S=S0
  i=1
  while(!any(I==0)){
    i=i+1
    I[i]=rbinom(1, size=S[i-1], prob=1-
                  exp(-beta*I[i-1]/S0))
    S[i]=S[i-1]-I[i]
  }
  out=data.frame(S=S, I=I)
  return(out)
}
plot(y, type="n", xlim=c(1,18),
     ylab="Predicted/observed", xlab="Week")
for(i in 1:100){
  sim=sim.cb(S0=floor(coef(fit)["S0"]),
             beta=coef(fit)["beta"], I0=11)
  lines(sim$I, col=grey(.5))
}
points(y, type="b", col=2)


#5.2 The Seasonally Forced SEIR Model
seirmod = function(t, y, parms) {
  S = y[1]
  E = y[2]
  I = y[3]
  R = y[4]
  mu = parms["mu"]
  N = parms["N"]
  beta = parms["beta"]
  sigma = parms["sigma"]
  gamma = parms["gamma"]
  dS = mu * (N - S) - beta * S * I/N
  dE = beta * S * I/N - (mu + sigma) * E
  dI = sigma * E - (mu + gamma) * I
  dR = gamma * I - mu * R
  res = c(dS, dE, dI, dR)
  list(res)
}
require(deSolve)
times = seq(0, 10, by=1/120)
paras = c(mu = 1/50, N = 1, beta = 1000,
          sigma = 365/8, gamma = 365/5)
start = c(S=0.06, E=0, I=0.001, R = 0.939)

R0 = expression(sigma/(sigma + mu) * beta/(gamma + mu))
with(as.list(paras), eval(R0))
out = as.data.frame(ode(start, times, seirmod, paras))
par(mfrow = c(1,2)) #Two plots side by side
plot(times, out$I, ylab = "Prevalence",
     xlab = "Time", type = "l")
plot(out$S, out$I, ylab = "Prevalence",
     xlab = "Susceptible", type = "l")

