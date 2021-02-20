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

#First variant.
FirstVariant <- data.frame()

for (i in 1:2){
  
  # Model inputs
  
  initial_state_values=c(S=9999,I=1,R=0) #change as desired
  parameters=c(gamma=1/14,beta=0.5) #change as desired
  S <- initial_state_values[1]
  I <- initial_state_values[2]
  R <- initial_state_values[3]
  
  # Time points
  
  time=seq(from=1,to=75,by=1) #change as desired
  
  # SIR model function 
  
  sir_model <- function(time,state,parameters){
    with(as.list(c(state,parameters)),{
      N=S+I+R
      lambda=beta*(I/N) 
      dS=-lambda*S
      dI=lambda*S-gamma*I
      dR=gamma*I
      
      return(list(c(dS,dI,dR)))
    }
    )
  }
  
  
  
  
  #Solving the differential equations
  output<-as.data.table(ode(y=initial_state_values,func = sir_model,parms=parameters,times = time))
  
  out_long=melt(output,id="time")
  
  # To plot the proportion of susceptible, infected and recovered individuals over time
  graph <- ggplot(data = out_long,          
                  aes(x = time, y = value/(S+I+R), colour = variable, group = variable)) +  
    geom_line() +xlab("Time (days)")+ylab("Proportion of the population")+scale_color_discrete(name="State")
  plot(graph)
  
  #Create data frame consisting of row vectors of output from SIR model
  
  Amount <- as.data.frame(output) %>% select(S)  #change SIR in select as desired
  transposeAmount<-t(Amount)
  FirstVariant<- rbind.data.frame(transposeAmount,FirstVariant)
}

#Apply Takens Embedding, calculate homology and determine if persistent for One Peak
homology <- data.frame()
for (i in 1:2){
  x<- data.matrix(FirstVariant[i,])
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
}


#Second Variant
SecondVariant <- data.frame()

for (i in 1:2){
  
  # Model inputs
  
  initial_state_values=c(S=9999,I=1,R=0) #change as desired
  parameters=c(gamma=0.0714,beta=1) #change as desired
  S <- initial_state_values[1]
  I <- initial_state_values[2]
  R <- initial_state_values[3]
  
  # Time points
  
  time=seq(from=1,to=75,by=1) #change as desired
  
  # SIR model function 
  
  sir_model <- function(time,state,parameters){
    with(as.list(c(state,parameters)),{
      N=S+I+R
      lambda=beta*(I/N) 
      dS=-lambda*S
      dI=lambda*S-gamma*I
      dR=gamma*I
      
      return(list(c(dS,dI,dR)))
    }
    )
  }
  
  
  
  
  #Solving the differential equations
  output<-as.data.table(ode(y=initial_state_values,func = sir_model,parms=parameters,times = time))
  
  out_long=melt(output,id="time")
  
  # To plot the proportion of susceptible, infected and recovered individuals over time
  graph <- ggplot(data = out_long,          
                  aes(x = time, y = value/(S+I+R), colour = variable, group = variable)) +  
    geom_line() +xlab("Time (days)")+ylab("Proportion of the population")+scale_color_discrete(name="State")
  plot(graph)
  
  #Create data frame consisting of row vectors of output from SIR model
  
  Amount <- as.data.frame(output) %>% select(S)  #change SIR in select as desired
  transposeAmount<-t(Amount)
  SecondVariant <- rbind.data.frame(transposeAmount,SecondVariant)
}

#Apply Takens Embedding, calculate homology and determine if persistent for One Peak
homology <- data.frame()
for (i in 1:2){
  x<- data.matrix(SecondVariant[i,])
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
}
