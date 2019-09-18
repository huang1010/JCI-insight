#Set Working Directory first


rm(list=ls())

library('MASS')
library('lubridate')

# Here we show 2 example individual data for chest temperature and rest-activity
# which were used to present averaged day chest surface temperature profile in Figure 3C

exp1.chest<-read.csv('chest.temp.RA.exp1.csv')
exp2.chest<-read.csv('chest.temp.RA.exp2.csv')

# and another 2 example individual data for core body temperature
# which were used to present  averaged day core body temperature profile in Figure 2C
exp1.core<-read.csv('core.temp.exp1.csv')
exp2.core<-read.csv('core.temp.exp2.csv')

#In the paper, we fitted the two-component cosinor model (with periods 12-h and 24-h) to describe the average day oscillation of both chest surface and core body temperature

#The following give an example how to perform cosinor analysis to the example 5-min data

#########################take exp1.chest as example#########################

info<-exp1.chest

#datetime: clock time of observations.  
datetime<-info$datetime
datetime<- as.POSIXct(datetime, format="%Y-%m-%d %H:%M:%S") #get appropriate time format 
#NOTE: format='' should follow your time format in the data recording

obs=info$chest.temp #here obs is chest surface temperature, may exsit NA values (allowed) 

#5-min resolution is used.
sf=1/12 # obs in 5-min data hence with sampling frequency sf=1/12 (hour)

#########################
#step 1: smooth the raw data
smooth_window=1 #in hours, smooth the raw data by hourly window

ma <- function(x, n = 1){stats::filter(x, rep(1 / n, n), sides = 2)}
obs_smoothed=ma(obs,n=smooth_window/sf)

##check the data 
#plot(obs); lines(obs_smoothed,col='red',lwd=2)



#########################
source('functions-two-component-cosinor.R') #all functions used in the following steps, see details in  functions-two-component-cosinor.R
#########################
#step 2: compute average day oscillation (24-h span) based on obs_smoothed
hour_day_start=12 #  should be in clockhour from 0 to 23, here we consider one day intervel start at 12pm, i.e. noon-noon(day+1)

cut_days=obs_in_days(obs=obs_smoothed,clocktime=datetime,hour_day_start=hour_day_start,sf=sf) 
# obs_in_days cuts the obs_smoothed in to days and compute the average day oscillation with starting day time =  hour_day_start
# for details, see in functions-two-component-cosinor.R 

day_average<-cut_days$day_average  #average day oscillation with one day starts at hour_day_start
#NOTE:time_oneday in day_average ranged from hour_day_start to hour_day_start+24-sf (NOT CLOCKTIME)

#plot(day_average$time_oneday,day_average$obs) #have a check



#########################
#step 3: two-component cosinor model (with period hmc1 and hmc2) analysis of averaged day oscillation
output<-cosinor_2hmcs(obs=day_average$obs,time=day_average$time_oneday,
                      hmc1=24,hmc2=12,uncertainty=T,n.b=100)
datafits<-output$datafits 
#$time:  time of a day, ranged from hour_day_start to hour_day_start+24-sf
#$obs:   average day oscillation input
#$fitting:   cosinor fitting for average day oscillation input

params<-output$params
#parameter estimation from two component cosinor model, i.e. mesor (mean level), amp (amplitude), phase1&2 (1st and 2nd acrophase if all exist), bathyphase1&2 (1st and 2nd bathyphase if all exist)

#have a check
plot(datafits$time,datafits$obs,col=gray(0.5)) 
lines(datafits$time,datafits$fitting,col=gray(0.2),lwd=2)  #cosinor fits
#abline(v=c(params$phase1,params$phase2),col='red') #the timing of curve fitting maxiaml 
#abline(v=c(params$bathyphase1,params$bathyphase1),col='blue') #the timing of curve fitting minimal 



params_boots<-output$params_boots
#bootstrapped parameters, which can be used to compute confidence interval of parameters.


#NOTE: phase and bathyphase in params and params_boots are the timing based on day_average$time_oneday. NEED TO BE CONVERTED TO CLOCKTIME if hour_day_start is not 0am
#for example, with hour_day_start=12, phase1=27 means acrophase=27-24 = 3 am