Drift_simulation <- function(num.species=2000, num.coms=100, num.years=1000, selection.strength, disp.rate) {
  library(vegan)
  
  year <- 1
  
  freq.1.mat <- matrix(nrow = num.species, ncol = num.coms*num.years) ##We will record growth with drift in this matrix##
  
  start.abundances <- ceiling(rlnorm(num.species, meanlog=4, sdlog=1.1)) #draw the means of the populations across the metacommunity
  for (i in 1:num.species) {
    freq.1.mat[i,1:num.coms] <- rnorm(num.coms, mean=start.abundances[i], sd=selection.strength)
  }
  
  freq.2.mat <- freq.1.mat ##duplicate matrix for growth without drift
  
  avg.J.nodrift <- vector(length = num.years)
  avg.J.drift <- vector(length = num.years)
  avg.S.drift <- vector(length = num.years)
  avg.S.nodrift <- vector(length = num.years)
  avg.even.drift <- vector(length = num.years)
  avg.even.nodrift <- vector(length = num.years)
  avg.BC.drift <- vector(length = num.years)
  avg.BC.nodrift <- vector(length = num.years)
  
  avg.BC.across <- vector(length=num.years)
  
  avg.J.nodrift[1] <- mean(colSums(freq.2.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]))
  avg.J.drift[1] <- mean(colSums(freq.1.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]))
  avg.S.drift[1] <- mean(colSums(freq.1.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]!=0))
  avg.S.nodrift[1] <- mean(colSums(freq.2.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]!=0))
  avg.even.drift[1] <-  mean(diversity(t(freq.1.mat)[(num.coms*year-(num.coms -1)):(num.coms*year),])/log(specnumber(t(freq.1.mat)[(num.coms*year-(num.coms -1)):(num.coms*year),])))
  avg.even.nodrift[1] <-  mean(diversity(t(freq.2.mat)[(num.coms*year-(num.coms -1)):(num.coms*year),])/log(specnumber(t(freq.2.mat)[(num.coms*year-(num.coms -1)):(num.coms*year),])))
  avg.BC.drift[1] <- mean(vegdist(t(log(freq.1.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]+1))))
  avg.BC.nodrift[1] <- mean(vegdist(t(log(freq.2.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]+1))))
  
  avg.BC.across[1] <- 0
  
  stdev.J.drift <- vector(length = num.years)
  stdev.S.drift <- vector(length = num.years)
  stdev.even.drift <- vector(length = num.years)
  stdev.J.nodrift <- vector(length = num.years)
  stdev.S.nodrift <- vector(length = num.years)
  stdev.even.nodrift <- vector(length = num.years)
  stdev.BC.drift <- vector(length = num.years)
  stdev.BC.nodrift <- vector(length = num.years)
  
  stdev.BC.across <- vector(length=num.years)
  
  stdev.J.drift[1] <- sd(colSums(freq.1.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]))
  stdev.S.drift[1] <- sd(colSums(freq.1.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]!=0))
  stdev.even.drift[1] <- sd(diversity(t(freq.1.mat)[(num.coms*year-(num.coms -1)):(num.coms*year),])/log(specnumber(t(freq.1.mat)[(num.coms*year-(num.coms -1)):(num.coms*year),])))
  stdev.BC.drift[1] <- sd(vegdist(t(log(freq.1.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]+1))))
  stdev.J.nodrift[1] <- sd(colSums(freq.2.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]))
  stdev.S.nodrift[1] <- sd(colSums(freq.2.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]!=0))
  stdev.even.nodrift[1] <- sd(diversity(t(freq.2.mat)[(num.coms*year-(num.coms -1)):(num.coms*year),])/log(specnumber(t(freq.2.mat)[(num.coms*year-(num.coms -1)):(num.coms*year),])))
  stdev.BC.nodrift[1] <- sd(vegdist(t(log(freq.2.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]+1))))
  
  stdev.BC.across[1] <- 0
  
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
    freq.1.mat[freq.1.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)] <= 0.5] <- 0 #species with 0.5 or less individuals get extinct
    
    ##Record drift-impacted communities
    avg.J.drift[j] <- mean(colSums(freq.1.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]))
    avg.S.drift[j] <- mean(colSums(freq.1.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]!=0))
    avg.even.drift[j] <- mean(diversity(t(freq.1.mat)[(num.coms*year-(num.coms -1)):(num.coms*year),])/log(specnumber(t(freq.1.mat)[(num.coms*year-(num.coms -1)):(num.coms*year),])))
    stdev.J.drift[j] <- sd(colSums(freq.1.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]))
    stdev.S.drift[j] <- sd(colSums(freq.1.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]!=0))
    stdev.even.drift[j] <- sd(diversity(t(freq.1.mat)[(num.coms*year-(num.coms -1)):(num.coms*year),])/log(specnumber(t(freq.1.mat)[(num.coms*year-(num.coms -1)):(num.coms*year),])))
    avg.BC.drift[j] <- mean(vegdist(t(log(freq.1.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]+1))))
    stdev.BC.drift[j] <- sd(vegdist(t(log(freq.1.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]+1))))
    
    #Grow drift-free communities
    freq.2.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)] <- freq.2.mat[, (num.coms*(year-2)+1):(num.coms*(year-1))]*growth.mat.nodrift
    
    #Dispersal in drift-free communities
    freq.2.disp <- disp.rate*freq.2.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)] #A (species * (comms-1)) matrix with the immigrant numbers. We will fix the first community in the end.
    freq.2.mat[, (num.coms*year-(num.coms -2)):(num.coms*year)] <- (1-disp.rate)*freq.2.mat[, (num.coms*year-(num.coms -2)):(num.coms*year)]+freq.2.disp[,1:(num.coms -1)] #Communities from the second to the last receive immigrants from the first to the one-before-last
    freq.2.mat[, (num.coms*year-(num.coms -1))] <- (1-disp.rate)*freq.2.mat[, (num.coms*year-(num.coms -1))]+freq.2.disp[,num.coms] ##First community gets immigrants from the last community
    
    #Extinction of populations with less than 0.5 individuals
    freq.2.mat[freq.2.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)] <= 0.5] <- 0 #species with 0.5 or less individuals get extinct
    
    #Record drift-free communities
    avg.J.nodrift[j] <- mean(colSums(freq.2.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]))
    avg.S.nodrift[j] <- mean(colSums(freq.2.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]!=0))
    avg.even.nodrift[j] <- diversity(t(freq.2.mat)[(num.coms*year-(num.coms -1)),])/log(specnumber(t(freq.2.mat)[(num.coms*year-(num.coms -1)),]))
    stdev.J.nodrift[j] <- sd(colSums(freq.2.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]))
    stdev.S.nodrift[j] <- sd(colSums(freq.2.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]!=0))
    stdev.even.nodrift[j] <- sd(diversity(t(freq.2.mat)[(num.coms*year-(num.coms -1)):(num.coms*year),])/log(specnumber(t(freq.2.mat)[(num.coms*year-(num.coms -1)):(num.coms*year),])))
    avg.BC.nodrift[j] <- mean(vegdist(t(log(freq.2.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]+1))))
    stdev.BC.nodrift[j] <- sd(vegdist(t(log(freq.2.mat[, (num.coms*year-(num.coms -1)):(num.coms*year)]+1))))
    
    #Trying to calculate BC_across
    temp.BC.across <- vector(length=num.coms)
    for (i in 1:num.coms) {
      temp.BC.across[i] <- mean(vegdist(t(log(data.frame((freq.1.mat[,(year-1)*100+i]+1),(freq.2.mat[,(year-1)*100+i]+1))))))
    }
    avg.BC.across[j] <- mean(temp.BC.across)
    stdev.BC.across[j] <- sd(temp.BC.across)
    
    rm(temp.BC.across)
    
    print(year)
    timestamp()
    
  }
  
  rm(freq.1.mat)
  rm(freq.2.mat)
  
  generations <- c(1:num.years)
  
  all_data <- data.frame(generations,avg.J.drift,avg.S.drift,avg.even.drift,stdev.J.drift,stdev.S.drift,stdev.even.drift,avg.BC.drift,stdev.BC.drift,avg.J.nodrift,avg.S.nodrift,avg.even.nodrift,stdev.J.nodrift,stdev.S.nodrift,stdev.even.nodrift,avg.BC.nodrift,stdev.BC.nodrift, avg.BC.across, stdev.BC.across)
  
  return(all_data)
}
