library(epimdr)
library(deSolve)
library(rootSolve)
library(phaseR)
library(shiny)

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
plot(times, exp(out[,3]), ylab="Infected", xlab="Time", ylim=c(0, 0.01), type='l')

paras2["u"] =2
xstart2 = log(c(S=.06, I=c(0.001, rep(0.0001/
                                        paras2["u"], paras2["u"]-1)), R = 0.0001))
out2 = as.data.frame(ode(xstart2, times, chainSIR,
                         paras2))
lines(times, apply(exp(out2[,-c(1:2,length(out2))]),
                   1 ,sum), col='blue')
paras2["u"] =73
xstart2 = log(c(S=.06, I=c(0.001, rep(0.0001/
                                        paras2["u"], paras2["u"]-1)), R = 0.0001))
out3 = as.data.frame(ode(xstart2, times, chainSIR,
                         paras2))
lines(times, apply(exp(out3[,-c(1:2,length(out3))]),
                   1, sum), col='red', lwd=2, lty=2)
paras2["u"] =500
xstart2 = log(c(S=.06, I=c(0.001, rep(0.0001/
                                        paras2["u"], paras2["u"]-1)), R = 0.0001))
out4 = as.data.frame(ode(xstart2, times, chainSIR,
                         paras2))
lines(times, apply(exp(out4[,-c(1:2,length(out4))]),
                   1,sum, na.rm=TRUE), col='green')
legend("topright", legend=c("SIR", "u=2", "u=500",
                            "u=73 (H-S)"), lty=c(1,1,1,2), lwd=c(1,1,1, 2),
       col=c("black", "blue", "green", "red"))



Amount <- as.data.frame(exp(out[,3]))
transposeAmount<-t(Amount)

Amount2 <- as.data.frame(exp(out2[,-c(1:2,length(out2))]))   
transposeAmount2<-t(Amount2)

Amount3 <- as.data.frame(exp(out3[,-c(1:2,length(out3))])) 
transposeAmount3<-t(Amount3)

Amount4 <- as.data.frame(exp(out4[,-c(1:2,length(out4))])) 
transposeAmount4<-t(Amount4)

x<- data.matrix(transposeAmount3)
a<- buildTakens(x,2,3)
plot(a)

