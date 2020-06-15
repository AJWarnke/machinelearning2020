###R Code###
rm(list=ls())
##Start with loading relevant Datasets
Kaub_latest=read.table(file="Tagesmittelwerte_Pegel_Kaub_20200607.txt",header=TRUE,dec=",")  ##Latest Dataset for the Pegel
attach(Kaub_latest)
Kaub_latest
level_t=ts(level, frequency = 365, start=c(2016,1))  ###Even though we have two Leap years (2016, 2020), we keep that in mind for the forecasts

setwd("C:/Users/PhiGe/Desktop/R")  ##Dataset of the rainfall in the catchment area of Rhine, ##I use the simple arithmetic mean of the measurement stations
library(readxl)
RAIN <- read_excel("Rain_all.xlsx")
attach(RAIN)
rainfall=ts(Wert, freq=365, start=c(2016,1))

############################################
###########Plots (Section 2)################
############################################
plot(level_t,ylim=c(0,700),main="Rhine Water Level",xlab="Year",ylab="Level in cm",col = 4)   ###This creates the Rhine water level plot
acf(level_t, main="ACF plot Rhine Water Level", lag.max = 30)   #This creates the acf plot
pacf(level_t, main="PACF plot Rhine Water Level",lag.max = 30) #This creates the pacf plot
lag.plot(level_t, lags=4, main="Lag Plots Rhine Water Level") #This creates the lag plot up to lag 4
plot(rainfall,main="Rainfall in the Rhine Catchment area",xlab="Year",ylab="Rainfall in mm",col = 4)   ###This creates the Rainfall plot

#############################################
#########BIC for AR(p) models################
#############################################
##Start by assessing a simple AR(p) model via the Bayes Information Criterion
lag_0=window(level_t, start=c(2016,11), end=c(2020,160))
lag_1=window(level_t, start=c(2016,10), end=c(2020,160-1))
lag_2=window(level_t, start=c(2016,9), end=c(2020,160-2))
lag_3=window(level_t, start=c(2016,8), end=c(2020,160-3))
lag_4=window(level_t, start=c(2016,7), end=c(2020,160-4))
lag_5=window(level_t, start=c(2016,6), end=c(2020,160-5))
lag_6=window(level_t, start=c(2016,5), end=c(2020,160-6))
lag_7=window(level_t, start=c(2016,4), end=c(2020,160-7))
lag_8=window(level_t, start=c(2016,3), end=c(2020,160-8))
lag_9=window(level_t, start=c(2016,2), end=c(2020,160-9))
lag_10=window(level_t, start=c(2016,1), end=c(2020,160-10))

AR_1=lm(lag_0~lag_1)
AR_2=lm(lag_0~lag_1+lag_2)
AR_3=lm(lag_0~lag_1+lag_2+lag_3)
AR_4=lm(lag_0~lag_1+lag_2+lag_3+lag_4)
AR_5=lm(lag_0~lag_1+lag_2+lag_3+lag_4+lag_5)
AR_6=lm(lag_0~lag_1+lag_2+lag_3+lag_4+lag_5+lag_6)
AR_7=lm(lag_0~lag_1+lag_2+lag_3+lag_4+lag_5+lag_6+lag_7)
AR_8=lm(lag_0~lag_1+lag_2+lag_3+lag_4+lag_5+lag_6+lag_7+lag_8)
AR_9=lm(lag_0~lag_1+lag_2+lag_3+lag_4+lag_5+lag_6+lag_7+lag_8+lag_9)
AR_10=lm(lag_0~lag_1+lag_2+lag_3+lag_4+lag_5+lag_6+lag_7+lag_8+lag_9+lag_10)

SRR_1=sum((AR_1$residuals)^2)
SRR_2=sum((AR_2$residuals)^2)
SRR_3=sum((AR_3$residuals)^2)
SRR_4=sum((AR_4$residuals)^2)
SRR_5=sum((AR_5$residuals)^2)
SRR_6=sum((AR_6$residuals)^2)
SRR_7=sum((AR_7$residuals)^2)
SRR_8=sum((AR_8$residuals)^2)
SRR_9=sum((AR_9$residuals)^2)
SRR_10=sum((AR_10$residuals)^2)
n=length(lag_0)

AR_BIC_1=log(SRR_1/n)+(1+1)*log(n)/n
AR_BIC_2=log(SRR_2/n)+(2+1)*log(n)/n
AR_BIC_3=log(SRR_3/n)+(3+1)*log(n)/n
AR_BIC_4=log(SRR_4/n)+(4+1)*log(n)/n
AR_BIC_5=log(SRR_5/n)+(5+1)*log(n)/n
AR_BIC_6=log(SRR_6/n)+(6+1)*log(n)/n
AR_BIC_7=log(SRR_7/n)+(7+1)*log(n)/n
AR_BIC_8=log(SRR_8/n)+(8+1)*log(n)/n
AR_BIC_9=log(SRR_9/n)+(9+1)*log(n)/n
AR_BIC_10=log(SRR_10/n)+(10+1)*log(n)/n

AR_BIC=c(AR_BIC_1,AR_BIC_2,AR_BIC_3,AR_BIC_4,AR_BIC_5,AR_BIC_6,AR_BIC_7,AR_BIC_8,AR_BIC_9,AR_BIC_10)
###The BIC for the different AR(p) models are 
AR_BIC
which.min(AR_BIC) ##The model with p=4 lags minimizes the Bayes Information criterion. For the test for stationary I use the AR(4) model

###Summary of the AR_4 model
summary(AR_4)
##adj.R^2=0.986  
library(sandwich)
SE_AR4=sqrt(diag(vcovHC(AR_4,type="HC1")))  ###Heteroskedastic Robust standard errors
SE_AR4

#############################################
#########Test for Stationarity###############
#############################################
#To perform the Augmented Dickey Fuller Test for Stationarity I use the transformed model (see paper)
change_t=diff(level_t)
change_0=window(change_t, start=c(2016,11), end=c(2020,160))
change_1=window(change_t, start=c(2016,10), end=c(2020,160-1))
change_2=window(change_t, start=c(2016,9), end=c(2020,160-2))
change_3=window(change_t, start=c(2016,8), end=c(2020,160-3))
ADF_model=lm(change_0~lag_1+change_1+change_2+change_3)
library(sandwich)
ADF_model_tstat=ADF_model$coefficients/sqrt(diag(vcovHC(ADF_model,type="const")))
ADF_model_tstat[2] ## The test statistic for the stationarity test is -5.774601. The critical value at the 1% level is -3.43, we thus reject the H_0 of nonstationarity
#I further perform the test for nonstationarity around a deterministic time trend (linear)
t=seq(1:n)
ADF_model_trend=lm(change_0~lag_1+t+change_1+change_2+change_3)
ADF_model_tstat_trend=ADF_model_trend$coefficients/sqrt(diag(vcovHC(ADF_model_trend,type="const")))
ADF_model_tstat_trend[2] ## The test statistic for the stationarity test is -5.83112  The critical value at the 1% level is -3.96, we thus reject the H_0 of nonstationarity


#############################################
###########Testing for Breaks################
#############################################
0.7*length(lag_0) ##1127 observations; I test the inner 70% of the observations for breaks 
0.15*length(lag_0) ##241.5 observations
QLR=rep(0, 1127)

for (i in 1:1127){
  D=c(rep(0,241+i),rep(1,1127+242-i))
  Dlag_1=D*lag_1
  Dlag_2=D*lag_2
  Dlag_3=D*lag_3
  Dlag_4=D*lag_4
  QLR_model=lm(lag_0~lag_1+lag_2+lag_3+lag_4+D+Dlag_1+Dlag_2+Dlag_3+Dlag_4)
  QLR_coef=QLR_model$coefficients
  library(sandwich)
  covhat=vcovHC(QLR_model, type="HC1")  ##Usage of heteroskedastic robust covariance
  R=cbind(matrix(0,5,5),diag(5))
  QLR[i]=t(R%*%QLR_coef)%*%solve(R%*%covhat%*%t(R))%*%(R%*%QLR_coef)/5
}
QLR_t=ts(QLR, freq=365, start=c(2016,243)) ##First break date that is testes is 2016,243 (30.08.2016); last break date is 2019,274 (2.Oktober 2019)(Schaltjahr beachten)
max(QLR_t)  ###QLR-Statistic is the maximum of all F test values, i.e. 2.185099. The critical value at the 5% level is 3.26, thus we do not reject the H_0 that there is no break at the 5% level
which.max(QLR_t) ###Maximum Value at 2018,328, i.e. 24.November 2018
plot(QLR_t,main="QLR Statistic (AR(4))", xlab="Year (August 2016 - Oktober 2019)", ylab="QLR Statistic", col=4)
###Maximum value occurs in 2018,328 24.November 2018

##To ensure that there is no break in the relationship at the very end of the sample (i.e. the last 15%), I use poos to assess the model accuracy
#First estimate the model for the 85% of the data.
lag85_0=window(level_t, start=c(2016,11), end=c(2019,284))
lag85_1=window(level_t, start=c(2016,10), end=c(2019,283))
lag85_2=window(level_t, start=c(2016,9), end=c(2019,282))
lag85_3=window(level_t, start=c(2016,8), end=c(2019,281))
lag85_4=window(level_t, start=c(2016,7), end=c(2019,280))
AR4_model85=lm(lag85_0~lag85_1+lag85_2+lag85_3+lag85_4)
summary(AR4_model85) ##The residual standard error is 12.35
AR4_model85_SE=sqrt(diag(vcovHC(AR4_model85, type="HC1")))
AR4_model85_SE

####Pseudo Out of Sample Forecast 
AR4_15percent=rep(0,242)
for (i in 0:241){
  lag15_0=window(level_t, start=c(2016,11), end=c(2019,284+i))
  lag15_1=window(level_t, start=c(2016,10), end=c(2019,283+i))
  lag15_2=window(level_t, start=c(2016,9), end=c(2019,282+i))
  lag15_3=window(level_t, start=c(2016,8), end=c(2019,281+i))
  lag15_4=window(level_t, start=c(2016,7), end=c(2019,280+i))
  AR4_model15=lm(lag15_0~lag15_1+lag15_2+lag15_3+lag15_4)
  n=length(lag15_0)
  AR4_15percent[i+1]=sum(AR4_model15$coefficients*c(1,lag15_0[n:(n-3)]))
}
AR5_15percent_t=ts(AR4_15percent, freq=365, start=c(2019,285))
actual_15=window(level_t, start=c(2019,285), end=c(2020,160))
RMSFE_15=sqrt(mean((actual_15-AR5_15percent_t)^2))
RMSFE_15  ###The RMSFE for the last 15 observations is 12.22573 and does not differ substantially  from the residual standard error
# I thus conclude that there is no change in the relationship at the end of the sample 


#############################################
###########Model Choice######################
#############################################

##Plotting the acf of the AR(4)-residuals to assess  whether the residuals are White noise/i.i.d.
acf(AR_4$residuals,lag.max=20, main="ACF plot of the AR(4) residuals")
##I perform the Box-Ljung Test for the first 20 lags
Box.test(AR_4$residuals, type="Ljung-Box", lag=20)  ##Test statistic of Ljung-Box test is 32.376 
qchisq(0.95,20) ##The Ljung-Box test statistic follows an chisquared distribution; 
#The critical value at the 5% level is 31.41043, we thus reject the hypothesis of i.i.d. at the 5% level
##We can, thus, conclude that the AR4 model is not sufficient to explain the variation in the residuals, therefore I test ADL model with the additional regressor of rainfall

###Estimation of Different ADL(p,p) models and computation of the BIC
rain_0=window(rainfall, start=c(2016,11), end=c(2020,160))
rain_1=window(rainfall, start=c(2016,10), end=c(2020,160-1))
rain_2=window(rainfall, start=c(2016,9), end=c(2020,160-2))
rain_3=window(rainfall, start=c(2016,8), end=c(2020,160-3))
rain_4=window(rainfall, start=c(2016,7), end=c(2020,160-4))
rain_5=window(rainfall, start=c(2016,6), end=c(2020,160-5))
rain_6=window(rainfall, start=c(2016,5), end=c(2020,160-6))
rain_7=window(rainfall, start=c(2016,4), end=c(2020,160-7))
rain_8=window(rainfall, start=c(2016,3), end=c(2020,160-8))
rain_9=window(rainfall, start=c(2016,2), end=c(2020,160-9))
rain_10=window(rainfall, start=c(2016,1), end=c(2020,160-10))


ADL_11=lm(lag_0~lag_1+rain_1)
ADL_22=lm(lag_0~lag_1+lag_2+rain_1+rain_2)
ADL_33=lm(lag_0~lag_1+lag_2+lag_3+rain_1+rain_2+rain_3)
ADL_44=lm(lag_0~lag_1+lag_2+lag_3+lag_4+rain_1+rain_2+rain_3+rain_4)
ADL_55=lm(lag_0~lag_1+lag_2+lag_3+lag_4+lag_5+rain_1+rain_2+rain_3+rain_4+rain_5)
ADL_66=lm(lag_0~lag_1+lag_2+lag_3+lag_4+lag_5+lag_6+rain_1+rain_2+rain_3+rain_4+rain_5+rain_6)
ADL_77=lm(lag_0~lag_1+lag_2+lag_3+lag_4+lag_5+lag_6+lag_7+rain_1+rain_2+rain_3+rain_4+rain_5+rain_6+rain_7)
ADL_88=lm(lag_0~lag_1+lag_2+lag_3+lag_4+lag_5+lag_6+lag_7+lag_8+rain_1+rain_2+rain_3+rain_4+rain_5+rain_6+rain_7+rain_8)
ADL_99=lm(lag_0~lag_1+lag_2+lag_3+lag_4+lag_5+lag_6+lag_7+lag_8++lag_9+rain_1+rain_2+rain_3+rain_4+rain_5+rain_6+rain_7+rain_8+rain_9)
ADL_1010=lm(lag_0~lag_1+lag_2+lag_3+lag_4+lag_5+lag_6+lag_7+lag_8+lag_9+lag_10+rain_1+rain_2+rain_3+rain_4+rain_5+rain_6+rain_7+rain_8+rain_9+rain_10)

ADL_SRR_1=sum((ADL_11$residuals)^2)
ADL_SRR_2=sum((ADL_22$residuals)^2)
ADL_SRR_3=sum((ADL_33$residuals)^2)
ADL_SRR_4=sum((ADL_44$residuals)^2)
ADL_SRR_5=sum((ADL_55$residuals)^2)
ADL_SRR_6=sum((ADL_66$residuals)^2)
ADL_SRR_7=sum((ADL_77$residuals)^2)
ADL_SRR_8=sum((ADL_88$residuals)^2)
ADL_SRR_9=sum((ADL_99$residuals)^2)
ADL_SRR_10=sum((ADL_1010$residuals)^2)
n=length(lag_0)

ADL_BIC_11=log(ADL_SRR_1/n)+2*log(n)/n
ADL_BIC_22=log(ADL_SRR_2/n)+4*log(n)/n
ADL_BIC_33=log(ADL_SRR_3/n)+6*log(n)/n
ADL_BIC_44=log(ADL_SRR_4/n)+8*log(n)/n
ADL_BIC_55=log(ADL_SRR_5/n)+10*log(n)/n
ADL_BIC_66=log(ADL_SRR_6/n)+12*log(n)/n
ADL_BIC_77=log(ADL_SRR_7/n)+14*log(n)/n
ADL_BIC_88=log(ADL_SRR_8/n)+16*log(n)/n
ADL_BIC_99=log(ADL_SRR_9/n)+18*log(n)/n
ADL_BIC_1010=log(ADL_SRR_10/n)+20*log(n)/n

ADL_BIC=c(ADL_BIC_11,ADL_BIC_22,ADL_BIC_33,ADL_BIC_44,ADL_BIC_55,ADL_BIC_66,ADL_BIC_77,ADL_BIC_88,ADL_BIC_99,ADL_BIC_1010)
which.min(ADL_BIC) #The model with 5 lags minimizes the BIC criterion
ADL_BIC


##Plotting the acf of the ADL(5,5)-residuals to assess whether the residuals are White noise/i.i.d.
acf(ADL_55$residuals,lag.max=20, main="ACF plot of the ADL(5,5) residuals")
##Not many autocorrelations are outside of the Bartlett bounds, which is an indication that the residuals are i.i.d. 
##I perform the Box-Ljung Test for the first 20 lags
Box.test(ADL_55$residuals, type="Ljung-Box", lag=20)  ##Test statistic of Ljung-Box test is 17.432 
qchisq(0.95,20) ##The Ljung-Box test statistic follows an chisquared distribution; 
#The critical value at the 5% level is 31.41043, we thus do not reject the hypothesis of i.i.d. at the 5% level
summary(ADL_55)
round(ADL_55$coefficients,5)
##Compared to the AR(4) model the adjusted R^2 increases from  0.986 to  0.9905
ADL_55_SE=sqrt(diag(vcovHC(ADL_55, type="HC1"))) ##Heteroskedastic Robust Standard Errors
ADL_55_SE

########################################################################
###########Test for Breaks - rainfall coefficients######################
########################################################################
QLR_ADL_r=rep(0,1127)
for (i in 1:1127){
  D=c(rep(0,241+i),rep(1,1127+242-i))
  Dlag_1=D*lag_1
  Dlag_2=D*lag_2
  Dlag_3=D*lag_3
  Dlag_4=D*lag_4
  Dlag_5=D*lag_5
  Drain_1=D*rain_1
  Drain_2=D*rain_2
  Drain_3=D*rain_3
  Drain_4=D*rain_4
  Drain_5=D*rain_5
  QLR_model=lm(lag_0~lag_1+lag_2+lag_3+lag_4+lag_5+rain_1+rain_2+rain_3+rain_4+rain_5+D+Dlag_1+Dlag_2+Dlag_3+Dlag_4+Dlag_5+Drain_1+Drain_2+Drain_3+Drain_4+Drain_5)
  QLR_coef=QLR_model$coefficients
  library(sandwich)
  covhat=vcovHC(QLR_model, type="HC1")  ##Usage of heteroskedastic robust covariance
  R=cbind(matrix(0,5,17),diag(5))
  QLR_ADL_r[i]=t(R%*%QLR_coef)%*%solve(R%*%covhat%*%t(R))%*%(R%*%QLR_coef)/5
}
QLR_ADL_r_t=ts(QLR_ADL_r, freq=365, start=c(2016,243)) ##First break date that is testes is 2016,243 (30.08.2016); last break date is 2019,274 (2.Oktober 2019)(Schaltjahr beachten)
max(QLR_ADL_r_t)  ###QLR-Statistic is the maximum of all F test values, i.e. 1.919126. The critical value at the 5% level is 3.26, thus we do not reject the H_0 that there is no break at the 5% level
which.max(QLR_ADL_r_t)
plot(QLR_ADL_r_t,main="QLR Statistic (ADL) Rainfall predictors", xlab="Year (August 2016 - Oktober 2019)", ylab="QLR Statistic", col=4)


###############################################################################
##Comparison of the models according the RMSFE criterion#######################
###############################################################################

#########################################################
#################RMSFE AR(4) model#######################
#########################################################
#This reproduces the RMSFE for the AR(4) model
AR_4_1=rep(0,171)
AR_4_2=rep(0,171)
AR_4_3=rep(0,171)
AR_4_4=rep(0,171)
AR_4_5=rep(0,171)
AR_4_6=rep(0,171)
AR_4_7=rep(0,171)
AR_4_8=rep(0,171)
AR_4_9=rep(0,171)
AR_4_10=rep(0,171)

length(window(lag_0, start=c(2019,355), end=c(2020,160))) ##171 observations in this time window

for (i in 0:170){
  lag_0=window(level_t, start=c(2016,30), end=c(2019,355+i))
  lag_1=window(level_t, start=c(2016,29), end=c(2019,354+i))
  lag_2=window(level_t, start=c(2016,28), end=c(2019,353+i))
  lag_3=window(level_t, start=c(2016,27), end=c(2019,352+i))
  lag_4=window(level_t, start=c(2016,26), end=c(2019,351+i))
  AR4_model=lm(lag_0~lag_1+lag_2+lag_3+lag_4)
  AR4_coef=AR4_model$coefficients
  n=length(lag_0)
  AR_4_1[i+1]=sum(AR4_coef*c(1,lag_0[n:(n-3)]))
  AR_4_2[i+1]=sum(AR4_coef*c(1,AR_4_1[i+1],lag_0[n:(n-2)]))
  AR_4_3[i+1]=sum(AR4_coef*c(1,AR_4_2[i+1],AR_4_1[i+1],lag_0[n:(n-1)]))
  AR_4_4[i+1]=sum(AR4_coef*c(1,AR_4_3[i+1],AR_4_2[i+1],AR_4_1[i+1],lag_0[n]))
  AR_4_5[i+1]=sum(AR4_coef*c(1,AR_4_4[i+1],AR_4_3[i+1],AR_4_2[i+1],AR_4_1[i+1]))
  AR_4_6[i+1]=sum(AR4_coef*c(1,AR_4_5[i+1],AR_4_4[i+1],AR_4_3[i+1],AR_4_2[i+1]))
  AR_4_7[i+1]=sum(AR4_coef*c(1,AR_4_6[i+1],AR_4_5[i+1],AR_4_4[i+1],AR_4_3[i+1]))
  AR_4_8[i+1]=sum(AR4_coef*c(1,AR_4_7[i+1],AR_4_6[i+1],AR_4_5[i+1],AR_4_4[i+1]))
  AR_4_9[i+1]=sum(AR4_coef*c(1,AR_4_8[i+1],AR_4_7[i+1],AR_4_6[i+1],AR_4_5[i+1]))
  AR_4_10[i+1]=sum(AR4_coef*c(1,AR_4_9[i+1],AR_4_8[i+1],AR_4_7[i+1],AR_4_6[i+1]))
}

AR_4_1_t=ts(AR_4_1, freq=365, start=c(2019,356))
AR_4_2_t=ts(AR_4_2, freq=365, start=c(2019,357))
AR_4_3_t=ts(AR_4_3, freq=365, start=c(2019,358))
AR_4_4_t=ts(AR_4_4, freq=365, start=c(2019,359))
AR_4_5_t=ts(AR_4_5, freq=365, start=c(2019,360))
AR_4_6_t=ts(AR_4_6, freq=365, start=c(2019,361))
AR_4_7_t=ts(AR_4_7, freq=365, start=c(2019,362))
AR_4_8_t=ts(AR_4_8, freq=365, start=c(2019,363))
AR_4_9_t=ts(AR_4_9, freq=365, start=c(2019,364))
AR_4_10_t=ts(AR_4_10, freq=365, start=c(2019,365))
actual=window(level_t, start=c(2019,356), end=c(2020,160))

AR_4_RMSFE_1=sqrt(mean((actual-AR_4_1_t)^2))
AR_4_RMSFE_2=sqrt(mean((actual-AR_4_2_t)^2))
AR_4_RMSFE_3=sqrt(mean((actual-AR_4_3_t)^2))
AR_4_RMSFE_4=sqrt(mean((actual-AR_4_4_t)^2))
AR_4_RMSFE_5=sqrt(mean((actual-AR_4_5_t)^2))
AR_4_RMSFE_6=sqrt(mean((actual-AR_4_6_t)^2))
AR_4_RMSFE_7=sqrt(mean((actual-AR_4_7_t)^2))
AR_4_RMSFE_8=sqrt(mean((actual-AR_4_8_t)^2))
AR_4_RMSFE_9=sqrt(mean((actual-AR_4_9_t)^2))
AR_4_RMSFE_10=sqrt(mean((actual-AR_4_10_t)^2))

AR_4_RMSFE=c(AR_4_RMSFE_1,AR_4_RMSFE_2,AR_4_RMSFE_3,AR_4_RMSFE_4,AR_4_RMSFE_5,AR_4_RMSFE_6,AR_4_RMSFE_7,AR_4_RMSFE_8,AR_4_RMSFE_9,AR_4_RMSFE_10)
AR_4_RMSFE



##########################################################
#################RMSFE VAR(5) model#######################
##########################################################
#This reproduce the RMSFE for the VAR(5) model
VAR_5_level_1=rep(0,171)
VAR_5_level_2=rep(0,171)
VAR_5_level_3=rep(0,171)
VAR_5_level_4=rep(0,171)
VAR_5_level_5=rep(0,171)
VAR_5_level_6=rep(0,171)
VAR_5_level_7=rep(0,171)
VAR_5_level_8=rep(0,171)
VAR_5_level_9=rep(0,171)
VAR_5_level_10=rep(0,171)
VAR_5_rain_1=rep(0,171)
VAR_5_rain_2=rep(0,171)
VAR_5_rain_3=rep(0,171)
VAR_5_rain_4=rep(0,171)
VAR_5_rain_5=rep(0,171)
VAR_5_rain_6=rep(0,171)
VAR_5_rain_7=rep(0,171)
VAR_5_rain_8=rep(0,171)
VAR_5_rain_9=rep(0,171)
VAR_5_rain_10=rep(0,171)

for (i in 0:170){
  lag_0=window(level_t, start=c(2016,30), end=c(2019,355+i))
  lag_1=window(level_t, start=c(2016,29), end=c(2019,354+i))
  lag_2=window(level_t, start=c(2016,28), end=c(2019,353+i))
  lag_3=window(level_t, start=c(2016,27), end=c(2019,352+i))
  lag_4=window(level_t, start=c(2016,26), end=c(2019,351+i))
  lag_5=window(level_t, start=c(2016,25), end=c(2019,350+i))
  rain_0=window(rainfall, start=c(2016,30), end=c(2019,355+i))
  rain_1=window(rainfall, start=c(2016,29), end=c(2019,354+i))
  rain_2=window(rainfall, start=c(2016,28), end=c(2019,353+i))
  rain_3=window(rainfall, start=c(2016,27), end=c(2019,352+i))
  rain_4=window(rainfall, start=c(2016,26), end=c(2019,351+i))
  rain_5=window(rainfall, start=c(2016,25), end=c(2019,350+i))
  VAR_5_level_model=lm(lag_0~lag_1+lag_2+lag_3+lag_4+lag_5+rain_1+rain_2+rain_3+rain_4+rain_5)
  VAR_5_rain_model=lm(rain_0~lag_1+lag_2+lag_3+lag_4+lag_5+rain_1+rain_2+rain_3+rain_4+rain_5)
  level_coef=VAR_5_level_model$coefficient
  rain_coef=VAR_5_rain_model$coefficient
  n=length(lag_0)
  VAR_5_level_1[i+1]=sum(level_coef[1:6]*c(1,lag_0[n:(n-4)]))+sum(level_coef[7:11]*c(rain_0[n:(n-4)]))
  VAR_5_rain_1[i+1]=sum(rain_coef[1:6]*c(1,lag_0[n:(n-4)]))+sum(rain_coef[7:11]*c(rain_0[n:(n-4)]))
  VAR_5_level_2[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_1[i+1],lag_0[n:(n-3)]))+sum(level_coef[7:11]*c(VAR_5_rain_1[i+1],rain_0[n:(n-3)]))
  VAR_5_rain_2[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_1[i+1],lag_0[n:(n-3)]))+sum(rain_coef[7:11]*c(VAR_5_rain_1[i+1],rain_0[n:(n-3)]))
  VAR_5_level_3[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_2[i+1],VAR_5_level_1[i+1],lag_0[n:(n-2)]))+sum(level_coef[7:11]*c(VAR_5_rain_2[i+1],VAR_5_rain_1[i+1],rain_0[n:(n-2)]))
  VAR_5_rain_3[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_2[i+1],VAR_5_level_1[i+1],lag_0[n:(n-2)]))+sum(rain_coef[7:11]*c(VAR_5_rain_2[i+1],VAR_5_rain_1[i+1],rain_0[n:(n-2)]))
  VAR_5_level_4[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_3[i+1],VAR_5_level_2[i+1],VAR_5_level_1[i+1],lag_0[n:(n-1)]))+sum(level_coef[7:11]*c(VAR_5_rain_3[i+1],VAR_5_rain_2[i+1],VAR_5_rain_1[i+1],rain_0[n:(n-1)]))
  VAR_5_rain_4[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_3[i+1],VAR_5_level_2[i+1],VAR_5_level_1[i+1],lag_0[n:(n-1)]))+sum(rain_coef[7:11]*c(VAR_5_rain_3[i+1],VAR_5_rain_2[i+1],VAR_5_rain_1[i+1],rain_0[n:(n-1)]))
  VAR_5_level_5[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_4[i+1],VAR_5_level_3[i+1],VAR_5_level_2[i+1],VAR_5_level_1[i+1],lag_0[n]))+sum(level_coef[7:11]*c(VAR_5_rain_4[i+1],VAR_5_rain_3[i+1],VAR_5_rain_2[i+1],VAR_5_rain_1[i+1],rain_0[n]))
  VAR_5_rain_5[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_4[i+1],VAR_5_level_3[i+1],VAR_5_level_2[i+1],VAR_5_level_1[i+1],lag_0[n]))+sum(rain_coef[7:11]*c(VAR_5_rain_4[i+1],VAR_5_rain_3[i+1],VAR_5_rain_2[i+1],VAR_5_rain_1[i+1],rain_0[n]))
  VAR_5_level_6[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_5[i+1],VAR_5_level_4[i+1],VAR_5_level_3[i+1],VAR_5_level_2[i+1],VAR_5_level_1[i+1]))+sum(level_coef[7:11]*c(VAR_5_rain_5[i+1],VAR_5_rain_4[i+1],VAR_5_rain_3[i+1],VAR_5_rain_2[i+1],VAR_5_rain_1[i+1]))
  VAR_5_rain_6[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_5[i+1],VAR_5_level_4[i+1],VAR_5_level_3[i+1],VAR_5_level_2[i+1],VAR_5_level_1[i+1]))+sum(rain_coef[7:11]*c(VAR_5_rain_5[i+1],VAR_5_rain_4[i+1],VAR_5_rain_3[i+1],VAR_5_rain_2[i+1],VAR_5_rain_1[i+1]))
  VAR_5_level_7[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_6[i+1],VAR_5_level_5[i+1],VAR_5_level_4[i+1],VAR_5_level_3[i+1],VAR_5_level_2[i+1]))+sum(level_coef[7:11]*c(VAR_5_rain_6[i+1],VAR_5_rain_5[i+1],VAR_5_rain_4[i+1],VAR_5_rain_3[i+1],VAR_5_rain_2[i+1]))
  VAR_5_rain_7[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_6[i+1],VAR_5_level_5[i+1],VAR_5_level_4[i+1],VAR_5_level_3[i+1],VAR_5_level_2[i+1]))+sum(rain_coef[7:11]*c(VAR_5_rain_6[i+1],VAR_5_rain_5[i+1],VAR_5_rain_4[i+1],VAR_5_rain_3[i+1],VAR_5_rain_2[i+1]))
  VAR_5_level_8[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_7[i+1],VAR_5_level_6[i+1],VAR_5_level_5[i+1],VAR_5_level_4[i+1],VAR_5_level_3[i+1]))+sum(level_coef[7:11]*c(VAR_5_rain_7[i+1],VAR_5_rain_6[i+1],VAR_5_rain_5[i+1],VAR_5_rain_4[i+1],VAR_5_rain_3[i+1]))
  VAR_5_rain_8[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_7[i+1],VAR_5_level_6[i+1],VAR_5_level_5[i+1],VAR_5_level_4[i+1],VAR_5_level_3[i+1]))+sum(rain_coef[7:11]*c(VAR_5_rain_7[i+1],VAR_5_rain_6[i+1],VAR_5_rain_5[i+1],VAR_5_rain_4[i+1],VAR_5_rain_3[i+1]))
  VAR_5_level_9[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_8[i+1],VAR_5_level_7[i+1],VAR_5_level_6[i+1],VAR_5_level_5[i+1],VAR_5_level_4[i+1]))+sum(level_coef[7:11]*c(VAR_5_rain_8[i+1],VAR_5_rain_7[i+1],VAR_5_rain_6[i+1],VAR_5_rain_5[i+1],VAR_5_rain_4[i+1]))
  VAR_5_rain_9[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_8[i+1],VAR_5_level_7[i+1],VAR_5_level_6[i+1],VAR_5_level_5[i+1],VAR_5_level_4[i+1]))+sum(rain_coef[7:11]*c(VAR_5_rain_8[i+1],VAR_5_rain_7[i+1],VAR_5_rain_6[i+1],VAR_5_rain_5[i+1],VAR_5_rain_4[i+1]))
  VAR_5_level_10[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_9[i+1],VAR_5_level_8[i+1],VAR_5_level_7[i+1],VAR_5_level_6[i+1],VAR_5_level_5[i+1]))+sum(level_coef[7:11]*c(VAR_5_rain_9[i+1],VAR_5_rain_8[i+1],VAR_5_rain_7[i+1],VAR_5_rain_6[i+1],VAR_5_rain_5[i+1]))
  VAR_5_rain_10[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_9[i+1],VAR_5_level_8[i+1],VAR_5_level_7[i+1],VAR_5_level_6[i+1],VAR_5_level_5[i+1]))+sum(rain_coef[7:11]*c(VAR_5_rain_9[i+1],VAR_5_rain_8[i+1],VAR_5_rain_7[i+1],VAR_5_rain_6[i+1],VAR_5_rain_5[i+1]))
  
}
actual=window(level_t, start=c(2019,356), end=c(2020,160))
VAR_5_level_1_t=ts(VAR_5_level_1, freq=365, start=c(2019,356))
VAR_5_level_2_t=ts(VAR_5_level_2, freq=365, start=c(2019,357))
VAR_5_level_3_t=ts(VAR_5_level_3, freq=365, start=c(2019,358))
VAR_5_level_4_t=ts(VAR_5_level_4, freq=365, start=c(2019,359))
VAR_5_level_5_t=ts(VAR_5_level_5, freq=365, start=c(2019,360))
VAR_5_level_6_t=ts(VAR_5_level_6, freq=365, start=c(2019,361))
VAR_5_level_7_t=ts(VAR_5_level_7, freq=365, start=c(2019,362))
VAR_5_level_8_t=ts(VAR_5_level_8, freq=365, start=c(2019,363))
VAR_5_level_9_t=ts(VAR_5_level_9, freq=365, start=c(2019,364))
VAR_5_level_10_t=ts(VAR_5_level_10, freq=365, start=c(2019,365))

VAR_5_RMSFE_1=sqrt(mean((actual-VAR_5_level_1_t)^2))
VAR_5_RMSFE_2=sqrt(mean((actual-VAR_5_level_2_t)^2))
VAR_5_RMSFE_3=sqrt(mean((actual-VAR_5_level_3_t)^2))
VAR_5_RMSFE_4=sqrt(mean((actual-VAR_5_level_4_t)^2))
VAR_5_RMSFE_5=sqrt(mean((actual-VAR_5_level_5_t)^2))
VAR_5_RMSFE_6=sqrt(mean((actual-VAR_5_level_6_t)^2))
VAR_5_RMSFE_7=sqrt(mean((actual-VAR_5_level_7_t)^2))
VAR_5_RMSFE_8=sqrt(mean((actual-VAR_5_level_8_t)^2))
VAR_5_RMSFE_9=sqrt(mean((actual-VAR_5_level_9_t)^2))
VAR_5_RMSFE_10=sqrt(mean((actual-VAR_5_level_10_t)^2))
VAR_5_RMSFE=c(VAR_5_RMSFE_1,VAR_5_RMSFE_2,VAR_5_RMSFE_3,VAR_5_RMSFE_4,VAR_5_RMSFE_5,VAR_5_RMSFE_6,VAR_5_RMSFE_7,VAR_5_RMSFE_8,VAR_5_RMSFE_9,VAR_5_RMSFE_10)
VAR_5_RMSFE
AR_4_RMSFE
###The RMSFE decreases substantially when the additional regressor of rainfall is included

##########################################################
#################RMSFE Auto ARIMA#########################
##########################################################
#This reproduce the RMSFE for the Auto.arima model

ARIMA_1=rep(0,171)
ARIMA_2=rep(0,171)
ARIMA_3=rep(0,171)
ARIMA_4=rep(0,171)
ARIMA_5=rep(0,171)
ARIMA_6=rep(0,171)
ARIMA_7=rep(0,171)
ARIMA_8=rep(0,171)
ARIMA_9=rep(0,171)
ARIMA_10=rep(0,171)

library(forecast)
for (i in 0:170){
  lag_0=window(level_t, start=c(2016,1), end=c(2019,355+i))
  res=auto.arima(lag_0) ##No additional Arguments used
  ARIMA_1[i+1]=predict(res, n.ahead=1)$pred
  ARIMA_2[i+1]=predict(res, n.ahead=2)$pred[2]
  ARIMA_3[i+1]=predict(res, n.ahead=3)$pred[3]
  ARIMA_4[i+1]=predict(res, n.ahead=4)$pred[4]
  ARIMA_5[i+1]=predict(res, n.ahead=5)$pred[5]
  ARIMA_6[i+1]=predict(res, n.ahead=6)$pred[6]
  ARIMA_7[i+1]=predict(res, n.ahead=7)$pred[7]
  ARIMA_8[i+1]=predict(res, n.ahead=8)$pred[8]
  ARIMA_9[i+1]=predict(res, n.ahead=9)$pred[9]
  ARIMA_10[i+1]=predict(res, n.ahead=10)$pred[10]
}

ARIMA_1_t=ts(ARIMA_1, freq=365, start=c(2019,356))
ARIMA_2_t=ts(ARIMA_2, freq=365, start=c(2019,356+1))
ARIMA_3_t=ts(ARIMA_3, freq=365, start=c(2019,356+2))
ARIMA_4_t=ts(ARIMA_4, freq=365, start=c(2019,356+3))
ARIMA_5_t=ts(ARIMA_5, freq=365, start=c(2019,356+4))
ARIMA_6_t=ts(ARIMA_6, freq=365, start=c(2019,356+5))
ARIMA_7_t=ts(ARIMA_7, freq=365, start=c(2019,356+6))
ARIMA_8_t=ts(ARIMA_8, freq=365, start=c(2019,356+7))
ARIMA_9_t=ts(ARIMA_9, freq=365, start=c(2019,356+8))
ARIMA_10_t=ts(ARIMA_10, freq=365, start=c(2019,356+9))
actual_t=window(level_t,start=c(2019,356), end=c(2019,356+170))

ARIMA_RMSFE_1=sqrt(mean((ARIMA_1_t-actual_t)^2))
ARIMA_RMSFE_2=sqrt(mean((ARIMA_2_t-actual_t)^2))
ARIMA_RMSFE_3=sqrt(mean((ARIMA_3_t-actual_t)^2))
ARIMA_RMSFE_4=sqrt(mean((ARIMA_4_t-actual_t)^2))
ARIMA_RMSFE_5=sqrt(mean((ARIMA_5_t-actual_t)^2))
ARIMA_RMSFE_6=sqrt(mean((ARIMA_6_t-actual_t)^2))
ARIMA_RMSFE_7=sqrt(mean((ARIMA_7_t-actual_t)^2))
ARIMA_RMSFE_8=sqrt(mean((ARIMA_8_t-actual_t)^2))
ARIMA_RMSFE_9=sqrt(mean((ARIMA_9_t-actual_t)^2))
ARIMA_RMSFE_10=sqrt(mean((ARIMA_10_t-actual_t)^2))

ARIMA_RMSFE=c(ARIMA_RMSFE_1,ARIMA_RMSFE_2,ARIMA_RMSFE_3,ARIMA_RMSFE_4,ARIMA_RMSFE_5,ARIMA_RMSFE_6,ARIMA_RMSFE_7,ARIMA_RMSFE_8,ARIMA_RMSFE_9,ARIMA_RMSFE_10)
ARIMA_RMSFE


###################################################
#################Forecast###########################
####################################################
lag_0=window(level_t, start=c(2016,11), end=c(2019,355+170))
lag_1=window(level_t, start=c(2016,10), end=c(2019,354+170))
lag_2=window(level_t, start=c(2016,9), end=c(2019,353+170))
lag_3=window(level_t, start=c(2016,8), end=c(2019,352+170))
lag_4=window(level_t, start=c(2016,7), end=c(2019,351+170))
lag_5=window(level_t, start=c(2016,6), end=c(2019,350+170))

rain_0=window(rainfall, start=c(2016,11), end=c(2019,355+170))
rain_1=window(rainfall, start=c(2016,10), end=c(2019,354+170))
rain_2=window(rainfall, start=c(2016,9), end=c(2019,353+170))
rain_3=window(rainfall, start=c(2016,8), end=c(2019,352+170))
rain_4=window(rainfall, start=c(2016,7), end=c(2019,351+170))
rain_5=window(rainfall, start=c(2016,6), end=c(2019,350+170))


VAR_55=lm(lag_0~lag_1+lag_2+lag_3+lag_4+lag_5+rain_1+rain_2+rain_3+rain_4+rain_5) ##I use the VAR(5) model since it has the minimized RMSFE
Rainy_55=lm(rain_0~lag_1+lag_2+lag_3+lag_4+lag_5+rain_1+rain_2+rain_3+rain_4+rain_5)
ADL_55_coef=VAR_55$coefficients
Rainy_55_coef=Rainy_55$coefficients
n=length(lag_0)
h_1=sum(ADL_55_coef[1:6]*c(1,lag_0[n:(n-4)]))+sum(ADL_55_coef[7:11]*c(rain_0[n:(n-4)]))
h_1_r=sum(Rainy_55_coef[1:6]*c(1,lag_0[n:(n-4)]))+sum(Rainy_55_coef[7:11]*c(rain_0[n:(n-4)]))
h_2=sum(ADL_55_coef[1:6]*c(1,h_1,lag_0[n:(n-3)]))+sum(ADL_55_coef[7:11]*c(h_1_r,rain_0[n:(n-3)]))
h_2_r=sum(Rainy_55_coef[1:6]*c(1,h_1,lag_0[n:(n-3)]))+sum(Rainy_55_coef[7:11]*c(h_1_r,rain_0[n:(n-3)]))
h_3=sum(ADL_55_coef[1:6]*c(1,h_2,h_1,lag_0[n:(n-2)]))+sum(ADL_55_coef[7:11]*c(h_2_r,h_1_r,rain_0[n:(n-2)]))
h_3_r=sum(Rainy_55_coef[1:6]*c(1,h_2,h_1,lag_0[n:(n-2)]))+sum(Rainy_55_coef[7:11]*c(h_2_r,h_1_r,rain_0[n:(n-2)]))
h_4=sum(ADL_55_coef[1:6]*c(1,h_3,h_2,h_1,lag_0[n:(n-1)]))+sum(ADL_55_coef[7:11]*c(h_3_r,h_2_r,h_1_r,rain_0[n:(n-1)]))
h_4_r=sum(Rainy_55_coef[1:6]*c(1,h_3,h_2,h_1,lag_0[n:(n-1)]))+sum(Rainy_55_coef[7:11]*c(h_3_r,h_2_r,h_1_r,rain_0[n:(n-1)]))
h_5=sum(ADL_55_coef[1:6]*c(1,h_4,h_3,h_2,h_1,lag_0[n]))+sum(ADL_55_coef[7:11]*c(h_4_r,h_3_r,h_2_r,h_1_r,rain_0[n]))
h_5_r=sum(Rainy_55_coef[1:6]*c(1,h_4,h_3,h_2,h_1,lag_0[n]))+sum(Rainy_55_coef[7:11]*c(h_4_r,h_3_r,h_2_r,h_1_r,rain_0[n]))
h_6=sum(ADL_55_coef[1:6]*c(1,h_5,h_4,h_3,h_2,h_1))+sum(ADL_55_coef[7:11]*c(h_5_r,h_4_r,h_3_r,h_2_r,h_1_r))
h_6_r=sum(Rainy_55_coef[1:6]*c(1,h_5,h_4,h_3,h_2,h_1))+sum(Rainy_55_coef[7:11]*c(h_5_r,h_4_r,h_3_r,h_2_r,h_1_r))
h_7=sum(ADL_55_coef[1:6]*c(1,h_6,h_5,h_4,h_3,h_2))+sum(ADL_55_coef[7:11]*c(h_6_r,h_5_r,h_4_r,h_3_r,h_2_r))
h_7_r=sum(Rainy_55_coef[1:6]*c(1,h_6,h_5,h_4,h_3,h_2))+sum(Rainy_55_coef[7:11]*c(h_6_r,h_5_r,h_4_r,h_3_r,h_2_r))
h_8=sum(ADL_55_coef[1:6]*c(1,h_7,h_6,h_5,h_4,h_3))+sum(ADL_55_coef[7:11]*c(h_7_r,h_6_r,h_5_r,h_4_r,h_3_r))
h_8_r=sum(Rainy_55_coef[1:6]*c(1,h_7,h_6,h_5,h_4,h_3))+sum(Rainy_55_coef[7:11]*c(h_7_r,h_6_r,h_5_r,h_4_r,h_3_r))
h_9=sum(ADL_55_coef[1:6]*c(1,h_8,h_7,h_6,h_5,h_4))+sum(ADL_55_coef[7:11]*c(h_8_r,h_7_r,h_6_r,h_5_r,h_4_r))
h_9_r=sum(Rainy_55_coef[1:6]*c(1,h_8,h_7,h_6,h_5,h_4))+sum(Rainy_55_coef[7:11]*c(h_8_r,h_7_r,h_6_r,h_5_r,h_4_r))
h_10=sum(ADL_55_coef[1:6]*c(1,h_9,h_8,h_7,h_6,h_5))+sum(ADL_55_coef[7:11]*c(h_9_r,h_8_r,h_7_r,h_6_r,h_5_r))
h_10_r=sum(Rainy_55_coef[1:6]*c(1,h_9,h_8,h_7,h_6,h_5))+sum(Rainy_55_coef[7:11]*c(h_9_r,h_8_r,h_7_r,h_6_r,h_5_r))
h=c(h_1,h_2,h_3,h_4,h_5,h_6,h_7,h_8,h_9,h_10)  ###Rhine water level forecasts
h
##Constructing the 90%,95%, 99% confidence intervals 

upper_95_h=h+qnorm(0.975)*VAR_5_RMSFE
lower_95_h=h-qnorm(0.975)*VAR_5_RMSFE
upper_90_h=h+qnorm(0.9)*VAR_5_RMSFE
lower_90_h=h-qnorm(0.9)*VAR_5_RMSFE
upper_99_h=h+qnorm(0.99)*VAR_5_RMSFE
lower_99_h=h-qnorm(0.99)*VAR_5_RMSFE

upper_95_h_t=ts(upper_95_h, freq=365, start=c(2020,161))
lower_95_h_t=ts(lower_95_h, freq=365, start=c(2020,161))
upper_90_h_t=ts(upper_90_h, freq=365, start=c(2020,161))
lower_90_h_t=ts(lower_90_h, freq=365, start=c(2020,161))
upper_99_h_t=ts(upper_99_h, freq=365, start=c(2020,161))
lower_99_h_t=ts(lower_99_h, freq=365, start=c(2020,161))



##############################################################
####################Fanchart Plot#############################
##############################################################
level_h=c(level_t, h)
level_h_t=ts(level_h, freq=365, start=c(2016,1))
x=seq(2020.161,2020.170,by=0.001)
z=seq(2020.140,2020.170,by=0.001)
level_new=window(level_h_t, start=c(2020,140), end=c(2020,170))
h_t=ts(h,freq=365,start=c(2020,161))
upper_95_h_t=ts(upper_95_h, freq=365, start=c(2020,161))
plot(z,level_new,type="n",yaxt="n",xaxt="n",ylim=c(-100,500),main="Fanchart Rhine Water Level", ylab="Rhin Water Level", xlab="Time")
lines(z,level_new)
length(level_new)
length(z)
axis(1,at=c(seq(2020.140,2020.170,by=0.004)),labels=c(seq(2020.140,2020.170,by=0.004)))
axis(2,at=c(seq(-100,500,by=50)),labels=c(seq(-100,500,by=50)))
max(level_new)
xx=c(2020.161,x,2020.170,2020.161,x,2020.170)
yy_u_95=c(0,upper_90_h_t,0,0,upper_95_h_t,0)
yy_l_95=c(0,upper_90_h_t,0,0,lower_95_h_t,0)
yy_u_90=c(0,h_t,0,0,upper_90_h_t,0)
yy_l_90=c(0,h_t,0,0,lower_90_h_t,0)
yy_u_99=c(0,upper_95_h_t,0,0,upper_99_h_t,0)
yy_l_99=c(0,lower_95_h_t,0,0,lower_99_h_t,0)
polygon(xx,yy_u_95,col=c("#3399FF"),border=F)
polygon(xx,yy_l_95,col=c("#3399FF"),border=F)
polygon(xx,yy_u_90,col=c("#0066CC"),border=F)
polygon(xx,yy_l_90,col=c("#0066CC"),border=F)
polygon(xx,yy_u_99,col=c("#99CCFF"),border=F)
polygon(xx,yy_l_99,col=c("#99CCFF"),border=F)
lines(x,h_t)

##############################################################
#########Results for different catchment areas################
##############################################################
#Note: I overwrite the existing variables for rain and copy the code before


#####################################################################
#################RMSFE VAR(5) Downstream model#######################
#####################################################################
#RMSFE for VAR(5) of the catchment area without downstream measurement stations
RAIN <- read_excel("Rain_reduced.xlsx")
attach(RAIN)
rainfall=ts(Wert, freq=365, start=c(2016,1))

VAR_5_level_1=rep(0,171)
VAR_5_level_2=rep(0,171)
VAR_5_level_3=rep(0,171)
VAR_5_level_4=rep(0,171)
VAR_5_level_5=rep(0,171)
VAR_5_level_6=rep(0,171)
VAR_5_level_7=rep(0,171)
VAR_5_level_8=rep(0,171)
VAR_5_level_9=rep(0,171)
VAR_5_level_10=rep(0,171)
VAR_5_rain_1=rep(0,171)
VAR_5_rain_2=rep(0,171)
VAR_5_rain_3=rep(0,171)
VAR_5_rain_4=rep(0,171)
VAR_5_rain_5=rep(0,171)
VAR_5_rain_6=rep(0,171)
VAR_5_rain_7=rep(0,171)
VAR_5_rain_8=rep(0,171)
VAR_5_rain_9=rep(0,171)
VAR_5_rain_10=rep(0,171)

for (i in 0:170){
  lag_0=window(level_t, start=c(2016,30), end=c(2019,355+i))
  lag_1=window(level_t, start=c(2016,29), end=c(2019,354+i))
  lag_2=window(level_t, start=c(2016,28), end=c(2019,353+i))
  lag_3=window(level_t, start=c(2016,27), end=c(2019,352+i))
  lag_4=window(level_t, start=c(2016,26), end=c(2019,351+i))
  lag_5=window(level_t, start=c(2016,25), end=c(2019,350+i))
  rain_0=window(rainfall, start=c(2016,30), end=c(2019,355+i))
  rain_1=window(rainfall, start=c(2016,29), end=c(2019,354+i))
  rain_2=window(rainfall, start=c(2016,28), end=c(2019,353+i))
  rain_3=window(rainfall, start=c(2016,27), end=c(2019,352+i))
  rain_4=window(rainfall, start=c(2016,26), end=c(2019,351+i))
  rain_5=window(rainfall, start=c(2016,25), end=c(2019,350+i))
  VAR_5_level_model=lm(lag_0~lag_1+lag_2+lag_3+lag_4+lag_5+rain_1+rain_2+rain_3+rain_4+rain_5)
  VAR_5_rain_model=lm(rain_0~lag_1+lag_2+lag_3+lag_4+lag_5+rain_1+rain_2+rain_3+rain_4+rain_5)
  level_coef=VAR_5_level_model$coefficient
  rain_coef=VAR_5_rain_model$coefficient
  n=length(lag_0)
  VAR_5_level_1[i+1]=sum(level_coef[1:6]*c(1,lag_0[n:(n-4)]))+sum(level_coef[7:11]*c(rain_0[n:(n-4)]))
  VAR_5_rain_1[i+1]=sum(rain_coef[1:6]*c(1,lag_0[n:(n-4)]))+sum(rain_coef[7:11]*c(rain_0[n:(n-4)]))
  VAR_5_level_2[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_1[i+1],lag_0[n:(n-3)]))+sum(level_coef[7:11]*c(VAR_5_rain_1[i+1],rain_0[n:(n-3)]))
  VAR_5_rain_2[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_1[i+1],lag_0[n:(n-3)]))+sum(rain_coef[7:11]*c(VAR_5_rain_1[i+1],rain_0[n:(n-3)]))
  VAR_5_level_3[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_2[i+1],VAR_5_level_1[i+1],lag_0[n:(n-2)]))+sum(level_coef[7:11]*c(VAR_5_rain_2[i+1],VAR_5_rain_1[i+1],rain_0[n:(n-2)]))
  VAR_5_rain_3[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_2[i+1],VAR_5_level_1[i+1],lag_0[n:(n-2)]))+sum(rain_coef[7:11]*c(VAR_5_rain_2[i+1],VAR_5_rain_1[i+1],rain_0[n:(n-2)]))
  VAR_5_level_4[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_3[i+1],VAR_5_level_2[i+1],VAR_5_level_1[i+1],lag_0[n:(n-1)]))+sum(level_coef[7:11]*c(VAR_5_rain_3[i+1],VAR_5_rain_2[i+1],VAR_5_rain_1[i+1],rain_0[n:(n-1)]))
  VAR_5_rain_4[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_3[i+1],VAR_5_level_2[i+1],VAR_5_level_1[i+1],lag_0[n:(n-1)]))+sum(rain_coef[7:11]*c(VAR_5_rain_3[i+1],VAR_5_rain_2[i+1],VAR_5_rain_1[i+1],rain_0[n:(n-1)]))
  VAR_5_level_5[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_4[i+1],VAR_5_level_3[i+1],VAR_5_level_2[i+1],VAR_5_level_1[i+1],lag_0[n]))+sum(level_coef[7:11]*c(VAR_5_rain_4[i+1],VAR_5_rain_3[i+1],VAR_5_rain_2[i+1],VAR_5_rain_1[i+1],rain_0[n]))
  VAR_5_rain_5[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_4[i+1],VAR_5_level_3[i+1],VAR_5_level_2[i+1],VAR_5_level_1[i+1],lag_0[n]))+sum(rain_coef[7:11]*c(VAR_5_rain_4[i+1],VAR_5_rain_3[i+1],VAR_5_rain_2[i+1],VAR_5_rain_1[i+1],rain_0[n]))
  VAR_5_level_6[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_5[i+1],VAR_5_level_4[i+1],VAR_5_level_3[i+1],VAR_5_level_2[i+1],VAR_5_level_1[i+1]))+sum(level_coef[7:11]*c(VAR_5_rain_5[i+1],VAR_5_rain_4[i+1],VAR_5_rain_3[i+1],VAR_5_rain_2[i+1],VAR_5_rain_1[i+1]))
  VAR_5_rain_6[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_5[i+1],VAR_5_level_4[i+1],VAR_5_level_3[i+1],VAR_5_level_2[i+1],VAR_5_level_1[i+1]))+sum(rain_coef[7:11]*c(VAR_5_rain_5[i+1],VAR_5_rain_4[i+1],VAR_5_rain_3[i+1],VAR_5_rain_2[i+1],VAR_5_rain_1[i+1]))
  VAR_5_level_7[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_6[i+1],VAR_5_level_5[i+1],VAR_5_level_4[i+1],VAR_5_level_3[i+1],VAR_5_level_2[i+1]))+sum(level_coef[7:11]*c(VAR_5_rain_6[i+1],VAR_5_rain_5[i+1],VAR_5_rain_4[i+1],VAR_5_rain_3[i+1],VAR_5_rain_2[i+1]))
  VAR_5_rain_7[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_6[i+1],VAR_5_level_5[i+1],VAR_5_level_4[i+1],VAR_5_level_3[i+1],VAR_5_level_2[i+1]))+sum(rain_coef[7:11]*c(VAR_5_rain_6[i+1],VAR_5_rain_5[i+1],VAR_5_rain_4[i+1],VAR_5_rain_3[i+1],VAR_5_rain_2[i+1]))
  VAR_5_level_8[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_7[i+1],VAR_5_level_6[i+1],VAR_5_level_5[i+1],VAR_5_level_4[i+1],VAR_5_level_3[i+1]))+sum(level_coef[7:11]*c(VAR_5_rain_7[i+1],VAR_5_rain_6[i+1],VAR_5_rain_5[i+1],VAR_5_rain_4[i+1],VAR_5_rain_3[i+1]))
  VAR_5_rain_8[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_7[i+1],VAR_5_level_6[i+1],VAR_5_level_5[i+1],VAR_5_level_4[i+1],VAR_5_level_3[i+1]))+sum(rain_coef[7:11]*c(VAR_5_rain_7[i+1],VAR_5_rain_6[i+1],VAR_5_rain_5[i+1],VAR_5_rain_4[i+1],VAR_5_rain_3[i+1]))
  VAR_5_level_9[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_8[i+1],VAR_5_level_7[i+1],VAR_5_level_6[i+1],VAR_5_level_5[i+1],VAR_5_level_4[i+1]))+sum(level_coef[7:11]*c(VAR_5_rain_8[i+1],VAR_5_rain_7[i+1],VAR_5_rain_6[i+1],VAR_5_rain_5[i+1],VAR_5_rain_4[i+1]))
  VAR_5_rain_9[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_8[i+1],VAR_5_level_7[i+1],VAR_5_level_6[i+1],VAR_5_level_5[i+1],VAR_5_level_4[i+1]))+sum(rain_coef[7:11]*c(VAR_5_rain_8[i+1],VAR_5_rain_7[i+1],VAR_5_rain_6[i+1],VAR_5_rain_5[i+1],VAR_5_rain_4[i+1]))
  VAR_5_level_10[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_9[i+1],VAR_5_level_8[i+1],VAR_5_level_7[i+1],VAR_5_level_6[i+1],VAR_5_level_5[i+1]))+sum(level_coef[7:11]*c(VAR_5_rain_9[i+1],VAR_5_rain_8[i+1],VAR_5_rain_7[i+1],VAR_5_rain_6[i+1],VAR_5_rain_5[i+1]))
  VAR_5_rain_10[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_9[i+1],VAR_5_level_8[i+1],VAR_5_level_7[i+1],VAR_5_level_6[i+1],VAR_5_level_5[i+1]))+sum(rain_coef[7:11]*c(VAR_5_rain_9[i+1],VAR_5_rain_8[i+1],VAR_5_rain_7[i+1],VAR_5_rain_6[i+1],VAR_5_rain_5[i+1]))
  
}
acutal=window(level_t, start=c(2019,356), end=c(2020,160))
VAR_5_level_1_t=ts(VAR_5_level_1, freq=365, start=c(2019,356))
VAR_5_level_2_t=ts(VAR_5_level_2, freq=365, start=c(2019,357))
VAR_5_level_3_t=ts(VAR_5_level_3, freq=365, start=c(2019,358))
VAR_5_level_4_t=ts(VAR_5_level_4, freq=365, start=c(2019,359))
VAR_5_level_5_t=ts(VAR_5_level_5, freq=365, start=c(2019,360))
VAR_5_level_6_t=ts(VAR_5_level_6, freq=365, start=c(2019,361))
VAR_5_level_7_t=ts(VAR_5_level_7, freq=365, start=c(2019,362))
VAR_5_level_8_t=ts(VAR_5_level_8, freq=365, start=c(2019,363))
VAR_5_level_9_t=ts(VAR_5_level_9, freq=365, start=c(2019,364))
VAR_5_level_10_t=ts(VAR_5_level_10, freq=365, start=c(2019,365))

VAR_5_RMSFE_1=sqrt(mean((actual-VAR_5_level_1_t)^2))
VAR_5_RMSFE_2=sqrt(mean((actual-VAR_5_level_2_t)^2))
VAR_5_RMSFE_3=sqrt(mean((actual-VAR_5_level_3_t)^2))
VAR_5_RMSFE_4=sqrt(mean((actual-VAR_5_level_4_t)^2))
VAR_5_RMSFE_5=sqrt(mean((actual-VAR_5_level_5_t)^2))
VAR_5_RMSFE_6=sqrt(mean((actual-VAR_5_level_6_t)^2))
VAR_5_RMSFE_7=sqrt(mean((actual-VAR_5_level_7_t)^2))
VAR_5_RMSFE_8=sqrt(mean((actual-VAR_5_level_8_t)^2))
VAR_5_RMSFE_9=sqrt(mean((actual-VAR_5_level_9_t)^2))
VAR_5_RMSFE_10=sqrt(mean((actual-VAR_5_level_10_t)^2))
VAR_5_RMSFE_reduced=c(VAR_5_RMSFE_1,VAR_5_RMSFE_2,VAR_5_RMSFE_3,VAR_5_RMSFE_4,VAR_5_RMSFE_5,VAR_5_RMSFE_6,VAR_5_RMSFE_7,VAR_5_RMSFE_8,VAR_5_RMSFE_9,VAR_5_RMSFE_10)
VAR_5_RMSFE_reduced #Column 6 in Table 5.4


##########################################################
#################RMSFE VAR(5) BAYERN model################
##########################################################
#RMSFE for VAR(5) of the catchment area without Bayern measurement stations
RAIN <- read_excel("Rain_bayern.xlsx")
attach(RAIN)
rainfall=ts(Wert, freq=365, start=c(2016,1))

VAR_5_level_1=rep(0,171)
VAR_5_level_2=rep(0,171)
VAR_5_level_3=rep(0,171)
VAR_5_level_4=rep(0,171)
VAR_5_level_5=rep(0,171)
VAR_5_level_6=rep(0,171)
VAR_5_level_7=rep(0,171)
VAR_5_level_8=rep(0,171)
VAR_5_level_9=rep(0,171)
VAR_5_level_10=rep(0,171)
VAR_5_rain_1=rep(0,171)
VAR_5_rain_2=rep(0,171)
VAR_5_rain_3=rep(0,171)
VAR_5_rain_4=rep(0,171)
VAR_5_rain_5=rep(0,171)
VAR_5_rain_6=rep(0,171)
VAR_5_rain_7=rep(0,171)
VAR_5_rain_8=rep(0,171)
VAR_5_rain_9=rep(0,171)
VAR_5_rain_10=rep(0,171)

for (i in 0:170){
  lag_0=window(level_t, start=c(2016,30), end=c(2019,355+i))
  lag_1=window(level_t, start=c(2016,29), end=c(2019,354+i))
  lag_2=window(level_t, start=c(2016,28), end=c(2019,353+i))
  lag_3=window(level_t, start=c(2016,27), end=c(2019,352+i))
  lag_4=window(level_t, start=c(2016,26), end=c(2019,351+i))
  lag_5=window(level_t, start=c(2016,25), end=c(2019,350+i))
  rain_0=window(rainfall, start=c(2016,30), end=c(2019,355+i))
  rain_1=window(rainfall, start=c(2016,29), end=c(2019,354+i))
  rain_2=window(rainfall, start=c(2016,28), end=c(2019,353+i))
  rain_3=window(rainfall, start=c(2016,27), end=c(2019,352+i))
  rain_4=window(rainfall, start=c(2016,26), end=c(2019,351+i))
  rain_5=window(rainfall, start=c(2016,25), end=c(2019,350+i))
  VAR_5_level_model=lm(lag_0~lag_1+lag_2+lag_3+lag_4+lag_5+rain_1+rain_2+rain_3+rain_4+rain_5)
  VAR_5_rain_model=lm(rain_0~lag_1+lag_2+lag_3+lag_4+lag_5+rain_1+rain_2+rain_3+rain_4+rain_5)
  level_coef=VAR_5_level_model$coefficient
  rain_coef=VAR_5_rain_model$coefficient
  n=length(lag_0)
  VAR_5_level_1[i+1]=sum(level_coef[1:6]*c(1,lag_0[n:(n-4)]))+sum(level_coef[7:11]*c(rain_0[n:(n-4)]))
  VAR_5_rain_1[i+1]=sum(rain_coef[1:6]*c(1,lag_0[n:(n-4)]))+sum(rain_coef[7:11]*c(rain_0[n:(n-4)]))
  VAR_5_level_2[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_1[i+1],lag_0[n:(n-3)]))+sum(level_coef[7:11]*c(VAR_5_rain_1[i+1],rain_0[n:(n-3)]))
  VAR_5_rain_2[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_1[i+1],lag_0[n:(n-3)]))+sum(rain_coef[7:11]*c(VAR_5_rain_1[i+1],rain_0[n:(n-3)]))
  VAR_5_level_3[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_2[i+1],VAR_5_level_1[i+1],lag_0[n:(n-2)]))+sum(level_coef[7:11]*c(VAR_5_rain_2[i+1],VAR_5_rain_1[i+1],rain_0[n:(n-2)]))
  VAR_5_rain_3[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_2[i+1],VAR_5_level_1[i+1],lag_0[n:(n-2)]))+sum(rain_coef[7:11]*c(VAR_5_rain_2[i+1],VAR_5_rain_1[i+1],rain_0[n:(n-2)]))
  VAR_5_level_4[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_3[i+1],VAR_5_level_2[i+1],VAR_5_level_1[i+1],lag_0[n:(n-1)]))+sum(level_coef[7:11]*c(VAR_5_rain_3[i+1],VAR_5_rain_2[i+1],VAR_5_rain_1[i+1],rain_0[n:(n-1)]))
  VAR_5_rain_4[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_3[i+1],VAR_5_level_2[i+1],VAR_5_level_1[i+1],lag_0[n:(n-1)]))+sum(rain_coef[7:11]*c(VAR_5_rain_3[i+1],VAR_5_rain_2[i+1],VAR_5_rain_1[i+1],rain_0[n:(n-1)]))
  VAR_5_level_5[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_4[i+1],VAR_5_level_3[i+1],VAR_5_level_2[i+1],VAR_5_level_1[i+1],lag_0[n]))+sum(level_coef[7:11]*c(VAR_5_rain_4[i+1],VAR_5_rain_3[i+1],VAR_5_rain_2[i+1],VAR_5_rain_1[i+1],rain_0[n]))
  VAR_5_rain_5[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_4[i+1],VAR_5_level_3[i+1],VAR_5_level_2[i+1],VAR_5_level_1[i+1],lag_0[n]))+sum(rain_coef[7:11]*c(VAR_5_rain_4[i+1],VAR_5_rain_3[i+1],VAR_5_rain_2[i+1],VAR_5_rain_1[i+1],rain_0[n]))
  VAR_5_level_6[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_5[i+1],VAR_5_level_4[i+1],VAR_5_level_3[i+1],VAR_5_level_2[i+1],VAR_5_level_1[i+1]))+sum(level_coef[7:11]*c(VAR_5_rain_5[i+1],VAR_5_rain_4[i+1],VAR_5_rain_3[i+1],VAR_5_rain_2[i+1],VAR_5_rain_1[i+1]))
  VAR_5_rain_6[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_5[i+1],VAR_5_level_4[i+1],VAR_5_level_3[i+1],VAR_5_level_2[i+1],VAR_5_level_1[i+1]))+sum(rain_coef[7:11]*c(VAR_5_rain_5[i+1],VAR_5_rain_4[i+1],VAR_5_rain_3[i+1],VAR_5_rain_2[i+1],VAR_5_rain_1[i+1]))
  VAR_5_level_7[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_6[i+1],VAR_5_level_5[i+1],VAR_5_level_4[i+1],VAR_5_level_3[i+1],VAR_5_level_2[i+1]))+sum(level_coef[7:11]*c(VAR_5_rain_6[i+1],VAR_5_rain_5[i+1],VAR_5_rain_4[i+1],VAR_5_rain_3[i+1],VAR_5_rain_2[i+1]))
  VAR_5_rain_7[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_6[i+1],VAR_5_level_5[i+1],VAR_5_level_4[i+1],VAR_5_level_3[i+1],VAR_5_level_2[i+1]))+sum(rain_coef[7:11]*c(VAR_5_rain_6[i+1],VAR_5_rain_5[i+1],VAR_5_rain_4[i+1],VAR_5_rain_3[i+1],VAR_5_rain_2[i+1]))
  VAR_5_level_8[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_7[i+1],VAR_5_level_6[i+1],VAR_5_level_5[i+1],VAR_5_level_4[i+1],VAR_5_level_3[i+1]))+sum(level_coef[7:11]*c(VAR_5_rain_7[i+1],VAR_5_rain_6[i+1],VAR_5_rain_5[i+1],VAR_5_rain_4[i+1],VAR_5_rain_3[i+1]))
  VAR_5_rain_8[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_7[i+1],VAR_5_level_6[i+1],VAR_5_level_5[i+1],VAR_5_level_4[i+1],VAR_5_level_3[i+1]))+sum(rain_coef[7:11]*c(VAR_5_rain_7[i+1],VAR_5_rain_6[i+1],VAR_5_rain_5[i+1],VAR_5_rain_4[i+1],VAR_5_rain_3[i+1]))
  VAR_5_level_9[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_8[i+1],VAR_5_level_7[i+1],VAR_5_level_6[i+1],VAR_5_level_5[i+1],VAR_5_level_4[i+1]))+sum(level_coef[7:11]*c(VAR_5_rain_8[i+1],VAR_5_rain_7[i+1],VAR_5_rain_6[i+1],VAR_5_rain_5[i+1],VAR_5_rain_4[i+1]))
  VAR_5_rain_9[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_8[i+1],VAR_5_level_7[i+1],VAR_5_level_6[i+1],VAR_5_level_5[i+1],VAR_5_level_4[i+1]))+sum(rain_coef[7:11]*c(VAR_5_rain_8[i+1],VAR_5_rain_7[i+1],VAR_5_rain_6[i+1],VAR_5_rain_5[i+1],VAR_5_rain_4[i+1]))
  VAR_5_level_10[i+1]=sum(level_coef[1:6]*c(1,VAR_5_level_9[i+1],VAR_5_level_8[i+1],VAR_5_level_7[i+1],VAR_5_level_6[i+1],VAR_5_level_5[i+1]))+sum(level_coef[7:11]*c(VAR_5_rain_9[i+1],VAR_5_rain_8[i+1],VAR_5_rain_7[i+1],VAR_5_rain_6[i+1],VAR_5_rain_5[i+1]))
  VAR_5_rain_10[i+1]=sum(rain_coef[1:6]*c(1,VAR_5_level_9[i+1],VAR_5_level_8[i+1],VAR_5_level_7[i+1],VAR_5_level_6[i+1],VAR_5_level_5[i+1]))+sum(rain_coef[7:11]*c(VAR_5_rain_9[i+1],VAR_5_rain_8[i+1],VAR_5_rain_7[i+1],VAR_5_rain_6[i+1],VAR_5_rain_5[i+1]))
  
}
acutal=window(level_t, start=c(2019,356), end=c(2020,160))
VAR_5_level_1_t=ts(VAR_5_level_1, freq=365, start=c(2019,356))
VAR_5_level_2_t=ts(VAR_5_level_2, freq=365, start=c(2019,357))
VAR_5_level_3_t=ts(VAR_5_level_3, freq=365, start=c(2019,358))
VAR_5_level_4_t=ts(VAR_5_level_4, freq=365, start=c(2019,359))
VAR_5_level_5_t=ts(VAR_5_level_5, freq=365, start=c(2019,360))
VAR_5_level_6_t=ts(VAR_5_level_6, freq=365, start=c(2019,361))
VAR_5_level_7_t=ts(VAR_5_level_7, freq=365, start=c(2019,362))
VAR_5_level_8_t=ts(VAR_5_level_8, freq=365, start=c(2019,363))
VAR_5_level_9_t=ts(VAR_5_level_9, freq=365, start=c(2019,364))
VAR_5_level_10_t=ts(VAR_5_level_10, freq=365, start=c(2019,365))

VAR_5_RMSFE_1=sqrt(mean((actual-VAR_5_level_1_t)^2))
VAR_5_RMSFE_2=sqrt(mean((actual-VAR_5_level_2_t)^2))
VAR_5_RMSFE_3=sqrt(mean((actual-VAR_5_level_3_t)^2))
VAR_5_RMSFE_4=sqrt(mean((actual-VAR_5_level_4_t)^2))
VAR_5_RMSFE_5=sqrt(mean((actual-VAR_5_level_5_t)^2))
VAR_5_RMSFE_6=sqrt(mean((actual-VAR_5_level_6_t)^2))
VAR_5_RMSFE_7=sqrt(mean((actual-VAR_5_level_7_t)^2))
VAR_5_RMSFE_8=sqrt(mean((actual-VAR_5_level_8_t)^2))
VAR_5_RMSFE_9=sqrt(mean((actual-VAR_5_level_9_t)^2))
VAR_5_RMSFE_10=sqrt(mean((actual-VAR_5_level_10_t)^2))
VAR_5_RMSFE_bayern=c(VAR_5_RMSFE_1,VAR_5_RMSFE_2,VAR_5_RMSFE_3,VAR_5_RMSFE_4,VAR_5_RMSFE_5,VAR_5_RMSFE_6,VAR_5_RMSFE_7,VAR_5_RMSFE_8,VAR_5_RMSFE_9,VAR_5_RMSFE_10)
VAR_5_RMSFE_bayern  #Column 7 in Table 5.4





