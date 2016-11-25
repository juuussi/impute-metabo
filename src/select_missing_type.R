#' Title: select_missing_type
# simulates the three types of missing data (mcar,mnar,mar)
# and the choice of different percentages of missingness 
#' @param data matrix to be used as reference
#' @param type indicates the type of missing data 1 == MCAR ,2 == MAR,3 == MNAR 
#' @param alphaP percentage of missingness
#'
#' @return data matrix with the missing values 
#' @export
#'
#' @examples
select_missing_type <- function(data,type,alphaP){
  
  ## a function that simulates the three types of missing data (mcar,mnar,mar)
  ## choice of different percentages of missingness 
  
  
  N <- length(simulated_data)
  if (type == 1){
  # mcar(missingness is random,flip of a coin)
    mcar <- rbinom(N,1,prob=alphaP)
    #simulate MCAR data
    simulation <-data*(1-mcar)+mcar*99999  
    simulation[simulation==99999] <- NA   #change 99999 to NA (R's notation for missing)
    
  }
  else if (type == 2){
    #mar : missing at random
    # create a depedence (missing is dependent on the observed values)
    mean_mar <- mean(x=data)
    sd_mar  <- sd(x=data)
    threshold   <- mean_mar+sd_mar*qnorm(alphaP)
    # simulate MAR data
    simulation <-  data
    simulation[simulation >threshold] <- NA
  
  }

  else if (type == 3){
    # mnar: missing not at random
    # dependence of missingness in both observed and unobserved values(left censored)
    mean_mnar <- colMeans(x=data)
    sd_mnar  <- apply(x=data, 2, sd)
    cut1   <- mean_mnar+sd_mnar*qnorm(alphaP)
    cut2   <- mean_mnar-sd_mnar*qnorm(alphaP)
    threshold <- cbind(cut1,cut2)
    # simulate MNAR
    simulation <-  data
    simulation[simulation > min(threshold) |simulation > max(threshold) ] <- NA
    }
  
  return  (simulation)
  
  
}

