#10/03/2021
#classification attempt - needs to be finished
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
FullData <- data.frame()
FirstVariant <- data.frame()

for (i in 1:10){
  # Model input
  Inft = sample(0:100,1)
  initial_values=c(S=(1000000-Inft),I=Inft,R=0)
  S <- initial_values[1]
  I <- initial_values[2]
  R <- initial_values[3]
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
                  aes(x = time, y = value/(S+I+R), colour = variable, group = variable)) +  
    geom_line() +                                                          
    xlab("Time (years)")+                         
    ylab("Prevalence") +scale_color_discrete(name="State")
  
  #plot(plotg)
  
  #Create data frame consisting of row vectors of output from SIR model
  
  Amount <- as.data.frame(output) %>% select(I)  #change SIR in select as desired
  transposeAmount<-t(Amount)
  FirstVariant<- rbind.data.frame(transposeAmount,FirstVariant)
}
plot(plotg)
#get average betti number for flute
average_list_flute <- c()
for (i in 1:10) {
  
  
  Flute_homology <- data.frame()
  
  flute_A4_matrix <- data.matrix(FirstVariant[i,])
  tak <- buildTakens(flute_A4_matrix,2,10)
  #plot(tak)
  hom <- calculate_homology(tak,return_df = TRUE) 
  hom_matrix <- calculate_homology(tak,return_df = FALSE) 
  hom <- hom %>%
    mutate(persistence = death-birth) %>%
    mutate(persistent = ifelse(persistence > 5000, 1,0))
  hom_matrix <- data.frame(hom) %>% select(dimension, persistent)
  hom_matrix <- as.data.frame(hom_matrix) 
  p1 <- hom_matrix[hom_matrix$persistent == '1',] 
  Flute_homology<-rbind.data.frame(Flute_homology,p1)
  avg <- nrow(Flute_homology)
  average_list_flute <- append(average_list_flute,avg)
}
plot(tak)
avg_betti_number_flute <- sum(average_list_flute)/length(average_list_flute)
avg_betti_number_flute

