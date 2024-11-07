##### Checking for sampling Uncertainities
##########Using random simulation of n=30 samples from the fitted GEV
# and using them to estimate parmaeters of 10,000.
# Then evaluating the RP100 10,000 tmes and later estimating the ratio for other durations
rm(list = ls(all = TRUE))
set.seed(123)
library(mev)
library(doParallel)
library(parallel)
library(foreach)
library(dplyr)
# comparing it with df.hot

################
# defining constants
n_simulations = 1000
RP = c(2,5,10,25,50,100,250,500)
irp = 1
prob = 1-(1./RP)
###################
# passing arguments from array job
input <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
for (input in 1:8){
scenario.label <- c('rcp45cooler_mid','rcp45hotter_mid','rcp85cooler_mid','rcp85hotter_mid','rcp45cooler_end','rcp45hotter_end','rcp85cooler_end','rcp85hotter_end')
scenario <- scenario.label[input]
#path = '/Users/njayadevan/nanditha/IDF/'
path = 'H:/nanditha/IDF/'
#*******************************
#* creating directories
dir.create(paste0(path,'ratios/',scenario), showWarnings = FALSE)
dir.create(paste0(path,'ratios/',scenario,'/RP',RP[irp]), showWarnings = FALSE)
#reading the latlons within US
latlon <- read.csv(paste0(path,'/latlon_clipped2.csv'))

# reading files
# loading daily GEV parameters
fut.daily1 <- read.csv(paste0(path,'parameters/',scenario,'/',scenario,'_hour_24.csv'))
hist.daily1 <- read.csv(paste0(path,'parameters/historical/historical_hour_24.csv'))

################################

### Saving the ratio of point estimate of RP100 precep every hour

subhourly_AF <- NULL
for (ihr in 1:24){
  fut1 <- read.csv(paste0(path,'parameters/',scenario,'/',scenario,'_hour_',ihr,'.csv'))
  hist1 <- read.csv(paste0(path,'parameters/historical/historical_hour_',ihr,'.csv'))
  # comparing the latitude and longitude values and picking only land grids
  fut <- inner_join(latlon,fut1, by = c('lat','lon'))
  hist <- inner_join(latlon,hist1, by = c('lat','lon'))
  temp.out <- NULL
  count1 <- count2 <- count3 <- count4 <- 0
  for (irow in 1:nrow(fut)){
    # Amplification factor for sub daily
    point.intensity.fut = mev::qgev(prob[irp],loc = fut$loc[irow], scale = fut$scale[irow],shape = fut$shape[irow],lower.tail = TRUE, log.p = FALSE )
    point.intensity.hist = mev::qgev(prob[irp],loc = hist$loc[irow], scale = hist$scale[irow],shape = hist$shape[irow],lower.tail = TRUE, log.p = FALSE )
    
    A.hourly = point.intensity.fut/point.intensity.hist
    temp.out <- rbind(temp.out,A.hourly)
  }
  subhourly_AF = cbind(subhourly_AF,temp.out)
  print(ihr)
}
colnames(subhourly_AF) = paste0('AF_',c(1:24))
subhourly_AF <- cbind(latlon[,1:2],subhourly_AF)
out.path = paste0(path,'ratios/',scenario,'/',scenario,'_AmplificationFactor_subdaily_PE_RP_',RP[irp],'.csv')
write.csv(subhourly_AF,out.path,row.names = FALSE)
print(input)
}
###########################################

hourly_ratios = read.csv(paste0(path,'/ratios/',scenario,'/',scenario,'_AmplificationFactor_subdaily_PE_RP_',RP[irp],'.csv'))
# considering only grids within CONUS

fut.daily <- inner_join(latlon,fut.daily1, by = c('lat','lon'))

hist.daily <- inner_join(latlon,hist.daily1, by = c('lat','lon'))


count1 <- count2 <- count3 <- count4 <- 0
####################
n.cores = 56
ind = round(nrow(fut.daily)/n.cores)
cl <- parallel::makeCluster(n.cores)
doParallel::registerDoParallel(cl)
###################

foreach::foreach(icore = c(1:n.cores)) %dopar% {
  strt = 1+(icore-1)*ind
  stop = ifelse(icore<n.cores, icore*ind,nrow(fut.daily))
  df.reject <- NULL
  df.reject.fs <- NULL
  df.ci <- NULL
  #rows = c(26  , 133   ,186   ,207  , 266  , 356 ,  389 ,  405  , 406,   442,   476  , 517 ,  519,   533 ,  545  , 554,   567   ,585  , 586  , 608)
  for (irow in strt:stop){
    # Generating Monte Carlo Simulations 1000 samples from the data of 30 year duration
    # estimating 1000 return levels for RP100
    rl.fut <- NULL
    for (iter in 1:n_simulations){
      repeat{
        ts = mev::rgev(30,loc = fut.daily$loc[irow], scale = fut.daily$scale[irow],shape = fut.daily$shape[irow])
        param <- tryCatch({
          mev::fit.gev(ts)
        }, error = function(e){
          NULL  # Return NULL if there's an error
        }) 
        if (!is.null(param)) {
          rl.fut <- rbind(rl.fut,mev::qgev(prob[irp],loc = param$estimate[1], scale = param$estimate[2],shape = param$estimate[3],lower.tail = TRUE, log.p = FALSE ))
          break
        }
        
      }
      #print(iter)
    }
    
    rl.hist <- NULL
    for (iter in 1:n_simulations){
      repeat {
        # Generate the time series
        ts <- mev::rgev(30, loc = hist.daily$loc[irow], scale = hist.daily$scale[irow], shape = hist.daily$shape[irow])
        
        # Try to estimate the parameters
        param <- tryCatch({
          mev::fit.gev(ts)
        }, error = function(e) {
          NULL  # Return NULL if there's an error
        })
        
        # Break the loop if the fitting was successful (param is not NULL)
        if (!is.null(param)) {
          rl.hist <- rbind(rl.hist,mev::qgev(prob[irp],loc = param$estimate[1], scale = param$estimate[2],shape = param$estimate[3],lower.tail = TRUE, log.p = FALSE ))
          break
        }
      }
    }
    
    
    
    
    ############################
    # test for comparison with daily precipitation
    
    # estimate the ci of the amplification factor and compare it with the point estimate
    # estimating the amplification factor for daily
    A.daily = rl.fut/rl.hist
    # Amplification factor for sub daily
    is.reject <- NULL
    for (icol in 1:24){
      A.hourly =  hourly_ratios[irow,icol+2]
      # CI for the sampling distribution and not the mean
      A.ratio = A.hourly/A.daily
      ratio.CI = quantile(A.ratio,c(0.025,0.975))
      # Check if all elements of ratio.CI are either < 1 or > 1
      if (all(ratio.CI < 1) | all(ratio.CI > 1)) {
        is.reject <- cbind(is.reject,1)
      } else {
        is.reject <- cbind(is.reject,0)
      }
    }
    
    # Considering the field significance, therefore 5% significance level becomes 0.2% when a buffer of 2500 grids considered
    is.reject.fs <- NULL
    for (icol in 1:24){
      A.hourly =  hourly_ratios[irow,icol+2]
      # CI for the sampling distribution and not the mean
      A.ratio = A.hourly/A.daily
      # considering 100 surrounding grids
      # 0.05% which is 0.025% and 99.975%
      ratio.CI = quantile(A.ratio,c(0.00025,0.99975))
      #ratio.CI = quantile(A.ratio,c(0.002,0.998))
      # Check if all elements of ratio.CI are either < 1 or > 1
      if (all(ratio.CI < 1) | all(ratio.CI > 1)) {
        is.reject.fs <- cbind(is.reject.fs,1)
      } else {
        is.reject.fs <- cbind(is.reject.fs,0)
      }
      
    }
    out <- cbind(mean(A.daily,na.rm = TRUE),sd(A.daily,na.rm = TRUE),t(quantile(A.daily,c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))))
    
    out.ratio.AF <- cbind(mean(A.ratio,na.rm = TRUE),sd(A.ratio,na.rm = TRUE),t(quantile(A.ratio,c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))))
    
    df.reject = rbind.data.frame(df.reject,cbind(lat = fut.daily$lat[irow],lon = fut.daily$lon[irow],
                                                 is.reject))
    df.reject.fs = rbind.data.frame(df.reject.fs,cbind(lat = fut.daily$lat[irow],lon = fut.daily$lon[irow],
                                                       is.reject.fs))
    
    df.ci = rbind.data.frame(df.ci,cbind(lat = fut.daily$lat[irow],lon = fut.daily$lon[irow],
                                         out))
    
    print(irow)
  }
  colnames(df.reject)<- c('lat','lon',paste0('hour',c(1:24)))
  colnames(df.reject.fs)<- c('lat','lon',paste0('hour',c(1:24)))
  out.path = paste0(path,'ratios/',scenario,'/RP',RP[irp],'/',scenario,'_rejection_indicator_RP',RP[irp],'_',icore,'.csv')
  write.csv(df.reject,out.path,row.names = FALSE)
  out.path = paste0(path,'ratios/',scenario,'/RP',RP[irp],'/',scenario,'_rejection_indicator_FS_05_RP',RP[irp],'_',icore,'.csv')
  write.csv(df.reject.fs,out.path,row.names = FALSE)
  
  colnames(df.ci) <- c('lat','lon','mean','sd','lowCI','5perc','25perc','median','75perc','95perc','upCI')
  # daily amplification factor
  out.path = paste0(path,'ratios/',scenario,'/RP',RP[irp],'/',scenario,'_ci_daily_AF_RP',RP[irp],'_',icore,'.csv') # wrote ci_ratio_AF initially
  write.csv(df.ci,out.path,row.names = FALSE)
  
}
parallel::stopCluster(cl)
