#8.5 SEIR with Error
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

times = seq(0, 10, by=1/52)
paras = c(mu = 1/50, N = 1, beta = 1000,
          sigma = 365/8, gamma = 365/5)
start = c(S=0.08, E=0, I=0.001, R = 0.939)
out = as.data.frame(ode(start, times, seirmod, paras))
datay = jitter(out$I, amount = 1e-04) #Add a small amount of noise to a numeric vector
plot(times, datay, ylab = "Infected", xlab = "Time")
lines(times, out$I, col = 2)

#Create data frame consisting of row vectors of output from SIR model
OnePeakData <- data.frame()

Amount <- as.data.frame(datay) #change SIR in select as desired
transposeAmount<-t(Amount)
OnePeakData <- rbind.data.frame(transposeAmount,OnePeakData)

#Apply Takens Embedding, calculate homology and determine if persistent for One Peak
homology <- data.frame()
for (i in 1:1){
  x<- data.matrix(OnePeakData[i,])
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
  
  flute_A4_matrix <- data.matrix(OnePeakData[i,])
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

