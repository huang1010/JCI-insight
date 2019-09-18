library('lubridate')
#Functions:

# Main function:cosinor_2hmcs: fit the data with two component (with period hmc1 and hmc2) cosinor model  
#' @param obs and time: data to fit with time. NOTE: time should be the absulute time of obs (NOT CLOCKTIME) with no decreasing values
#' @param hmc1 & hmc2: periods in hours
#' @param uncertainty: if TRUE, compute the bootstrapped parameter estimates  by performing pseduo bootstrap trails (resampling the fitting residuals for n.b times) 
#' @param n.b: number of bootstrap trails applied if uncertainty=T
#' @return datafits: time, obs and corresponding curve fitting
#' @return params: paramters, i.e. mesor, amplitude, 1st and 2nd (if exist) acrophase and bathyphase 
#' @return params_boots: bootstrapped paramters
#NOTE: acrophase and bathyphase in params and params_boots are the timing based on time input. NEED TO BE CONVERTED TO CLOCKHOUR if hour_day_start is not 0am

#' 
#Other functions 
#find_day_start: find index in clockhour time where clocktime[index]=hour_day_start
#obs_in_days: cuts the data in to days and compute the average day oscillation with starting day time =  hour_day_start
#curve_by_pars:generate two component cosinor fitting curves by parameters estimated in cosinor_2hmcs
#extract_pars: extract amplitude and phase informaion from the curve fitting
#boots_cosinor_2hmcs: compute bootstrapped  cosinor_2hmcs paramater estimations 
find_day_start<-function(hour_day_start,clocktime,sf){
  hour_time<-hour(clocktime)
  min_time<-minute(clocktime)
  A<-which(hour_time==hour_day_start)
  B<-which(min_time<60*sf) #
  one_day_start<-A[min(which(A%in%B, arr.ind = TRUE))]  
  return(one_day_start)
}

obs_in_days<-function(obs,clocktime,hour_day_start=12,sf=1/12){
  day_length<-24/sf
  n_days<-floor(length(obs)/day_length)+1
  obs_in_day<-matrix(NA,nrow=day_length,ncol=(n_days+1))
  
  col_name<-c()
  for (i in 1:n_days){
    day_start<-(i-1)*day_length+1
    day_end<-min((day_start+day_length-1),length(obs))
    size_day<-length(1:(day_end-day_start+1))
    obs_in_day[,(i+1)][1:size_day]<-obs[day_start:day_end]
    col_name[i]<-paste0('day',i)
  }
  obs_in_day[,1]<-rowMeans(obs_in_day[,2:(n_days+1)], na.rm = TRUE)
  obs_in_day<-data.frame(obs_in_day)
  col_name<-c('day average',col_name)
  colnames(obs_in_day)<-col_name
  obs_in_day$clocktime=clocktime[1:day_length]
  
  
  #select slot by hour_day_start
  start_hour_index<-find_day_start(hour_day_start=hour_day_start,clocktime=clocktime,sf=sf)
  obs_twoday<-c(obs_in_day[,1],obs_in_day[,1])
  
  profile<-obs_twoday[start_hour_index:(start_hour_index+day_length-1)]
  time_index_oneday=seq(hour_day_start,(hour_day_start+24-sf),sf) #here we use decimal hours
  day_average=data.frame(time_oneday=time_index_oneday,obs=profile) 
  return(list(day_average=day_average,obs_in_day=obs_in_day))
}


curve_by_pars<-function(pars,hmc1,hmc2,time){
  #generate two component cosinor fitting curves by pars estimated in cosinor_2hmcs
  pars<-as.matrix(unname(pars))
  dim(pars)<-c(5,1)
  X <- cbind(rep(1,length(time)), sin(time*2*pi/hmc1), cos(time*2*pi/hmc1),
             sin(time*2*pi/hmc2), cos(time*2*pi/hmc2))
  opt_fit <- as.vector(X%*%pars)
  return(opt_fit)
}

extract_pars<-function(time,curve){
  acrophase1<-acrophase2<-bathyphase1<-bathyphase2<-NA
  amp<-(max(curve,na.rm = T)-min(curve,na.rm = T))/2
  acrophase1<-time[which.max(curve)]
  bathyphase1<-time[which.min(curve)]
  
  peaks.locs <- pracma::findpeaks(curve, sortstr=T,npeaks=2)
  
  if(length(nrow(peaks.locs))>0){
    if(nrow(peaks.locs)==2){
      acrophase2<-time[peaks.locs[2,2]]  }  
  }
  
  min.locs <- pracma::findpeaks(-curve, sortstr=T,npeaks=2)
  if(length(nrow(min.locs))>0){
    if(nrow(min.locs)==2){
   bathyphase2<-time[min.locs[2,2]]}
  }
  

  summarypars=data.frame(amp=amp,acrophase1=acrophase1,acrophase2=acrophase2,
                         bathyphase1=bathyphase1,bathyphase2=bathyphase2)
  return(summarypars)
}

boots_cosinor_2hmcs<-function(fittings,residuals, time,hmc1,hmc2,n.b){
  mesor.boots<-amp.boots<-bathyphase1.boots<-bathyphase2.boots<-acrophase1.boots<-acrophase2.boots<-rep(NA,n.b)
  
  fitted_index=which(!is.na(fittings))
  L<-length(fitted_index)

  for (b in 1:n.b){
    residulas_b<-rep(NA,length(fittings))
    index<-sample(1:L, L, replace = TRUE)
    residulas_b[fitted_index]<-residuals[index]
    obs_b<-fittings+residulas_b

    fit.b<-lm(obs_b~sin(time*2*pi/hmc1)+cos(time*2*pi/hmc1)+sin(time*2*pi/hmc2)+cos(time*2*pi/hmc2),na.action = na.exclude)
    pars.b<-fit.b$coefficients
    mesor.boots[b]<-pars.b[1]
    curve.b<-curve_by_pars(pars.b,hmc1,hmc2,time)
    
    summary_pars.b=extract_pars(time,curve.b)
    amp.boots[b]<-summary_pars.b$amp
    acrophase1.boots[b]<-summary_pars.b$acrophase1
    acrophase2.boots[b]<-summary_pars.b$acrophase2
    bathyphase1.boots[b]<-summary_pars.b$bathyphase1
    bathyphase2.boots[b]<-summary_pars.b$bathyphase2
    
  }
  
  boots_pars<-data.frame(mesor=mesor.boots,amp=amp.boots, acrophase1=acrophase1.boots,acrophase2=acrophase2.boots,
                         bathyphase1=bathyphase1.boots,bathyphase2=bathyphase2.boots)
  return(boots_pars)
}

cosinor_2hmcs<-function(obs,time, hmc1=24,hmc2=12,uncertainty=T,n.b=1000){

  fit<-lm(obs~sin(time*2*pi/hmc1)+cos(time*2*pi/hmc1)+sin(time*2*pi/hmc2)+cos(time*2*pi/hmc2),na.action = na.exclude)
  pars<-fit$coefficients #fitting coefficients

  curve<-curve_by_pars(pars,hmc1,hmc2,time) #generate fitting curves by pars
  params<-extract_pars(time,curve)
  params$mesor<-pars[1]
  params<-params[,c('mesor',"amp","acrophase1","acrophase2","bathyphase1",'bathyphase2')]

  
  pars_boot=NA
  if(uncertainty){
  fitted_index<-which(!is.na(obs))
  fittings<-residuals<-rep(NA,length(obs))
  fittings[fitted_index]<-fit$fitted.values
  residuals[fitted_index]<-fit$residuals
  params_boots<-boots_cosinor_2hmcs(fittings,residuals, time,hmc1,hmc2,n.b)
  }
  
  datafits<-data.frame(time=time,obs=obs,fitting=curve)
  return(list(datafits=datafits,params=params,params_boots=params_boots))
}



