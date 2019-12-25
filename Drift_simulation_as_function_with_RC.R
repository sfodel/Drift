Drift_simulation_RC <- function(num.species=2000, num.coms=100, num.years=1000, selection.strength, disp.rate, disp.lim=TRUE) {
  
  ##This script grows communities with and without drift just like the main script and calculates the proportion of communities under dispersal limitation or homogenizing dispersal every 100 years
  ##To run this version of the script, you need first to load the modified RC script from James Stegen where you can replace "dist" with "bcdist" from package ecodist to make it run a bit faster. Please make sure your script is named <<raup_crick_abundance>> so that it can be called properly
  ##Moreover, you have to load the simple script <<check.integer>> that you can find here: https://stackoverflow.com/questions/3476782/check-if-the-number-is-integer
  ##If disp.lim is TRUE, the script returns the proportion of dispersal-limited community pairs, otherwise it returns the proportion of community pairs under homogenizing dispersal
  
  library(vegan)
  library(ecodist)
  year <- 1
  
  freq.1.mat <- matrix(nrow = num.species, ncol = num.coms*num.years) ##We will record growth with drift in this matrix##
  
  start.abundances <- ceiling(rlnorm(num.species, meanlog=4, sdlog=1.1)) #draw the means of the populations across the metacommunity
  for (i in 1:num.species) {
    freq.1.mat[i,1:num.coms] <- rnorm(num.coms, mean=start.abundances[i], sd=selection.strength)
  }
  
  freq.2.mat <- freq.1.mat ##duplicate matrix for growth without drift
  

  avg.DL.drift <- vector(length = num.years/100)
  avg.DL.nodrift <- vector(length = num.years/100)
  
  avg.DL.drift[1] <- 0
  avg.DL.nodrift[1] <- 0
  
  
  for (j in 2:num.years) {
    
    ##record initial data into the community matrices##
    
    
    year <- year +1 
    
    ## specify drift based on its distribution according to experimental data. Numbers here are CVs##
    drift <- rlnorm(num.species, meanlog = 0.625, sdlog=0.889) 
    rand <- sample(num.species,0.806*num.species) #create random number sequences selecting from a string of length=num.species. Drift hits randomly the remaining 16% of the communities.
    drift[rand] <- 0 #replace the drift values with 0s at the random positions
    
    ##specify the growth rates at each generation##
    growth.mat.nodrift <- matrix(nrow = length(start.abundances), ncol = num.coms) ##empty matrix to start with##
    growth.mat.drift <- matrix(nrow = length(start.abundances), ncol = num.coms)  ##empty matrix to start with##
    
    for (i in 1:num.coms) {
      growth.mat.nodrift[,i] <- rnorm(num.species, mean =1, sd=selection.strength)  ## growth rates are different for each community and are normally distributed. The community is not growing on average (mean=1). The environmental heterogeneity depends on selection.strength.
    }
    
    hit.drift <- growth.mat.drift #replicate the empty matrix
    
    for (i in 1:num.coms) { 
      hit.drift[,i] <- rnorm(num.species, mean=0, sd=growth.mat.nodrift[,i]*drift/100) #drift hits the same populations across communities
    }
    
    ##Grow the communities according to the growth rates with drift
    growth.mat.drift <- growth.mat.nodrift*(1-hit.drift)
    
    #Grow drift-impacted communities
    freq.1.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)] <- freq.1.mat[, (num.coms*(year-2)+1):(num.coms*(year-1))]*growth.mat.drift
    
    #Dispersal in drift-impacted communities
    freq.1.disp <- disp.rate*freq.1.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)] #A (species * (comms-1)) matrix with the immigrant numbers. We will fix the first community in the end.
    freq.1.mat[, (num.coms*year-(num.coms -2)):(num.coms*year)] <- (1-disp.rate)*freq.1.mat[, (num.coms*year-(num.coms -2)):(num.coms*year)]+freq.1.disp[,1:(num.coms -1)] #Communities from the second to the last receive immigrants from the first to the one-before-last
    freq.1.mat[, (num.coms*year-(num.coms -1))] <- (1-disp.rate)*freq.1.mat[, (num.coms*year-(num.coms -1))]+freq.1.disp[,num.coms] ##First community gets immigrants from the last community
    
    #Extinction of populations with less than 0.5 individuals
    freq.1.mat[freq.1.mat[, (num.coms*year-(num.coms-1)):(num.coms*year)] <= 0.5] <- 0 #species with 0.5 or less individuals get extinct
    
    ##Record drift-impacted communities
    if(check.integer(j/100) == TRUE){
    RC_temp_d <- raup_crick_abundance(spXsite = (log(t(freq.1.mat[, (num.coms*year-(num.coms-1)):(num.coms*year)])+1)), plot_names_in_col1 = FALSE, reps =99, as.distance.matrix = TRUE)
    if (disp.lim ==TRUE) {
      avg.DL.drift[j/100] <- sum(RC_temp_d > 0.95) / length(RC_temp_d) #proportion of dispersal-limited community pairs
    }
    if (disp.lim ==FALSE) {
      avg.DL.drift[j/100] <- sum(RC_temp_d < -0.95) / length(RC_temp_d) #proportion of dispersal-homogenized community pairs
      }
      rm(RC_temp_d)}
    
    #Grow drift-free communities
    freq.2.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)] <- freq.2.mat[, (num.coms*(year-2)+1):(num.coms*(year-1))]*growth.mat.nodrift
    
    #Dispersal in drift-free communities
    freq.2.disp <- disp.rate*freq.2.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)] #A (species * (comms-1)) matrix with the immigrant numbers. We will fix the first community in the end.
    freq.2.mat[, (num.coms*year-(num.coms -2)):(num.coms*year)] <- (1-disp.rate)*freq.2.mat[, (num.coms*year-(num.coms -2)):(num.coms*year)]+freq.2.disp[,1:(num.coms -1)] #Communities from the second to the last receive immigrants from the first to the one-before-last
    freq.2.mat[, (num.coms*year-(num.coms -1))] <- (1-disp.rate)*freq.2.mat[, (num.coms*year-(num.coms -1))]+freq.2.disp[,num.coms] ##First community gets immigrants from the last community
    
    #Extinction of populations with less than 0.5 individuals
    freq.2.mat[freq.2.mat[, (num.coms*year-(num.coms-1)):(num.coms*year)] <= 0.5] <- 0 #species with 0.5 or less individuals get extinct
    
    #Record drift-free communities
    if(check.integer(j/100) == TRUE){
    RC_temp_nd <- raup_crick_abundance(spXsite = (log(t(freq.2.mat[, (num.coms*year-(num.coms-1)):(num.coms*year)])+1)), plot_names_in_col1 = FALSE, reps =99, as.distance.matrix = TRUE)
    if (disp.lim ==TRUE) {
      avg.DL.nodrift[j/100] <- sum(RC_temp_nd > 0.95) / length(RC_temp_nd) #proportion of dispersal-limited community pairs
    }
    if (disp.lim ==FALSE) {
      avg.DL.nodrift[j/100] <- sum(RC_temp_nd < -0.95) / length(RC_temp_nd) #proportion of dispersal-homogenized community pairs
    }
    rm(RC_temp_nd)}
    
    print(year)
    timestamp()
    
  }
  
  rm(freq.1.mat)
  rm(freq.2.mat)
  
  generations <- c(100,200,300,400,500,600,700,800,900,1000)
  
  all_data <- data.frame(generations, avg.DL.drift, avg.DL.nodrift)
  
  return(all_data)
}
