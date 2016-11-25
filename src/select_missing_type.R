## create a function that simulates the three types of missing data(mcar,mnar,mar)
## choice of different percenatges

select_missing_type <- function(data,type,alphaP){
  N <- length(simulated_data)
  if (type == 1){
  # mcar
    mcar <- rbinom(N,1,prob=alphaP)
    simulation <-data*(1-mcar)+mcar*99999  
    simulation[simulation==99999]=NA                  #change 99999 to NA (R's notation for missing)
    
  }
  else if (type == 2){
    #mar 
    # create a depedence
    

    mean_mar <- mean(x=data)
    sd_mar  <- sd(x=data)
    cut   <- mean_mar+sd_mar*qnorm(alphaP)
    simulation <-  data
    simulation[simulation<cut] <- NA
  
  }

  # else if (type == 3){
  #   simulation
  # }
  # 
  # else if (type = 4){
  #   simulation
  # }
  
  return  (simulation)
  
  
}

