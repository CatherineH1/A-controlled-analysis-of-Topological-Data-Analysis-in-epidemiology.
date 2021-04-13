#Load required libraries

library(ggplot2)
library(deSolve)
library(reshape2)

#ODE 
library(deSolve)

#Graphs and Data Manipulation
library(tidyverse)
library(data.table)

#Confusion Matrix
library(caret) 

#Persistent Homology
library(nonlinearTseries)
library(TDAstats)

#Betti Numbers
library(matrixStats)

#databinding
library(dplyr)

#Undriven
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
out = as.data.frame(ode(start, times, seirmod, paras))
par(mfrow = c(1,2)) #Two plots side by side
#plot(times, out$I, ylab = "Infected",
     xlab = "Time", type = "l")
#Create data frame consisting of row vectors of output from SIR model
Dat<- data.frame()

#change SIR in select as desired
Amount <- as.data.frame(out$I)
transposeAmount<-t(Amount)
Dat <- rbind.data.frame(transposeAmount, Dat)


#Apply Takens Embedding, calculate homology and determine if persistent for One Peak
homology <- data.frame()
for (i in 1:1){
  x<- data.matrix(Dat[i,])
  a<- buildTakens(x,2,10)
  plot(a, xlim =c(0,1e-03), ylim =c(0, 11e-04))
  hom <- calculate_homology(a,return_df = TRUE)
  hom <- hom %>%
    mutate(persistence = death-birth) %>%
    mutate(persistent = ifelse(persistence > max(persistence)-0.00001, 1,0))
  hom_matrix <- data_frame(hom) %>% select(dimension, persistent)
  hom_matrix <- as.data.frame(hom_matrix) 
  p1 <- hom_matrix[hom_matrix$persistent == '1',] 
  homology <- rbind.data.frame(homology,p1) 
}

pers <- calculate_homology(a)
par(mfrow = c(1,2))
plot_barcode(pers)
plot_persist(pers)










#Seasonal
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
times = seq(0, 10, by=1/120)
paras = c(mu = 1/50, N = 1, beta0 = 1000, beta1 = 0.2,
          sigma = 365/8, gamma = 365/5)
start = c(S=0.06, E=0, I=0.001, R = 0.939)
out = as.data.frame(ode(start, times, seirmod2, paras))
par(mfrow=c(1,2)) #Side-by-side plot
plot(times, out$I, ylab="Infected", xlab="Time")

#Create data frame consisting of row vectors of output from SIR model
OnePeakData <- data.frame()

#change SIR in select as desired
Amount <- as.data.frame(out$I)
transposeAmount<-t(Amount)
OnePeakData <- rbind.data.frame(transposeAmount,OnePeakData)

#Apply Takens Embedding, calculate homology and determine if persistent for One Peak
homology <- data.frame()

x<- data.matrix(OnePeakData[1,])
a<- buildTakens(x,2,3)
plot(a)
hom <- calculate_homology(a,return_df = TRUE)
hom <- hom %>%
  mutate(persistence = death-birth) %>%
  mutate(persistent = ifelse(persistence > max(persistence)-0.00001, 1,0))
hom_matrix <- data_frame(hom) %>% select(dimension, persistent)
hom_matrix <- as.data.frame(hom_matrix) 
p1 <- hom_matrix[hom_matrix$persistent == '1',] 
homology <- rbind.data.frame(homology,p1) 

pers <- calculate_homology(a)
par(mfrow = c(1,2))
plot_barcode(pers)
plot_persist(pers)

