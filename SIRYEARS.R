#Changing value of gamma(the recovery rate) and showing the difference
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

FirstVariant <- data.frame()

for (i in 1:2){
# Model input

initial_values=c(S=999999,I=1,R=0)
parameters=c(gamma=0.2*365,beta=0.45*365,sigma=1/(2))

# Time points

time=seq(from=1,to=5,by=1/365)

#SIR model function
sir_model4 <- function(time,state,parameters){
  with(as.list(c(state,parameters)),{
    N=S+I+R
    lambda=beta*(I/N)
    dS=-lambda*S+sigma*R
    dI=lambda*S-gamma*I
    dR=gamma*I-sigma*R
    
    return(list(c(dS,dI,dR)))
  }
  )
}



# Solving the differential equations:
output<-as.data.table(ode(y=initial_values,func = sir_model4,parms=parameters,times = time))

out_long=melt(output,id="time")
#Prevalence plot
plotg <- ggplot(data = out_long,          
       aes(x = time, y = value/1000000, colour = variable, group = variable)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                         
  ylab("Prevalence") +scale_color_discrete(name="State")

plot(plotg)

#Create data frame consisting of row vectors of output from SIR model

Amount <- as.data.frame(output) %>% select(I)  #change SIR in select as desired
transposeAmount<-t(Amount)
FirstVariant<- rbind.data.frame(transposeAmount,FirstVariant)
}

#Apply Takens Embedding, calculate homology and determine if persistent for One Peak
homology <- data.frame()
for (i in 1:2){
  x<- data.matrix(FirstVariant[i,])
  a<- buildTakens(x,2,15)
  plot(a)
  hom <- calculate_homology(a,return_df = TRUE)
  hom <- hom %>%
    mutate(persistence = death-birth) %>%
    mutate(persistent = ifelse(persistence > max(persistence)-0.00001, 1,0))
  hom_matrix <- data_frame(hom) %>% select(dimension, persistent)
  hom_matrix <- as.data.frame(hom_matrix) 
  p1 <- hom_matrix[hom_matrix$persistent == '1',] 
  homology <- rbind.data.frame(homology,p1) 
}

#here we go again
SecondVariant <- data.frame()

for (i in 1:2){
  # Model input
  
  initial_values=c(S=999999,I=1,R=0)
  parameters=c(gamma=0.3*365,beta=0.45*365,sigma=1/(2))
  
  # Time points
  
  time=seq(from=1,to=5,by=1/365)
  
  #SIR model function
  sir_model4 <- function(time,state,parameters){
    with(as.list(c(state,parameters)),{
      N=S+I+R
      lambda=beta*(I/N)
      dS=-lambda*S+sigma*R
      dI=lambda*S-gamma*I
      dR=gamma*I-sigma*R
      
      return(list(c(dS,dI,dR)))
    }
    )
  }
  
  
  
  # Solving the differential equations:
  output<-as.data.table(ode(y=initial_values,func = sir_model4,parms=parameters,times = time))
  
  out_long=melt(output,id="time")
  #Prevalence plot
  plotg <- ggplot(data = out_long,          
                  aes(x = time, y = value/1000000, colour = variable, group = variable)) +  
    geom_line() +                                                          
    xlab("Time (years)")+                         
    ylab("Prevalence") +scale_color_discrete(name="State")
  
  plot(plotg)
  
  #Create data frame consisting of row vectors of output from SIR model
  
  Amount <- as.data.frame(output) %>% select(I)  #change SIR in select as desired
  transposeAmount<-t(Amount)
  SecondVariant<- rbind.data.frame(transposeAmount,SecondVariant)
}

#Apply Takens Embedding, calculate homology and determine if persistent for One Peak
homology <- data.frame()
for (i in 1:2){
  x<- data.matrix(SecondVariant[i,])
  a<- buildTakens(x,2,15)
  plot(a)
  hom <- calculate_homology(a,return_df = TRUE)
  hom <- hom %>%
    mutate(persistence = death-birth) %>%
    mutate(persistent = ifelse(persistence > max(persistence)-0.00001, 1,0))
  hom_matrix <- data_frame(hom) %>% select(dimension, persistent)
  hom_matrix <- as.data.frame(hom_matrix) 
  p1 <- hom_matrix[hom_matrix$persistent == '1',] 
  homology <- rbind.data.frame(homology,p1) 
}

