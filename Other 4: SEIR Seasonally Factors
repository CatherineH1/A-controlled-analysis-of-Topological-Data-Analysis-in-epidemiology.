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
          sigma = 365/8, gamma = 365/5)
start = c(S=0.06, E=0, I=0.001, R = 0.939)
out = as.data.frame(ode(start, times, seirmod2, paras))
par(mfrow=c(1,2)) #Side-by-side plot
plot(times, out$I, ylab = "Infected",
     xlab = "Time", type = "l")

#Create data frame consisting of row vectors of output from SIR model
Data <- data.frame()

Amount <- as.data.frame(out) %>% select(I)  #change SIR in select as desired
transposeAmount<-t(Amount)
Data <- rbind.data.frame(transposeAmount,Data)

#Apply Takens Embedding, calculate homology and determine if persistent for One Peak
homology <- data.frame()
for (i in 1:1){
  x<- data.matrix(Data[i,])
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

average_list_flute <- c()
for (i in 1:30) {
  Flute_homology <- data.frame()
  
  flute_A4_matrix <- data.matrix(Data[i,])
  tak <- buildTakens(flute_A4_matrix,2,10)
  hom <- calculate_homology(tak,return_df = TRUE) 
  hom_matrix <- calculate_homology(tak,return_df = FALSE) 
  hom <- hom %>%
    mutate(persistence = death-birth) %>%
    mutate(persistent = ifelse(persistence > 19000, 1,0))
  hom_matrix <- data.frame(hom) %>% select(dimension, persistent)
  hom_matrix <- as.data.frame(hom_matrix) 
  p1 <- hom_matrix[hom_matrix$persistent == '1',] 
  Flute_homology<-rbind.data.frame(Flute_homology,p1)
  avg <- nrow(Flute_homology)
  average_list_flute <- append(average_list_flute,avg)
}
avg_betti_number_flute <- sum(average_list_flute)/length(average_list_flute)
avg_betti_number_flute
