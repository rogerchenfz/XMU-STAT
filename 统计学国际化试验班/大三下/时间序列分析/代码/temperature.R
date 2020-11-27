## Global temperature
rm(list=ls())
library(forecast)
library(zoo)
Gt=read.table("D:/teaching/TSA/TSA2019ÑÐ¾¿Éú/Rcode2019/m-GLBTs.txt")
Gtemp=ts(Gt$V3,frequency=12,start=c(1880,1))
plot(Gtemp,xlab='year',ylab='temperature',type='l') # Plot the data
acf(Gtemp,lag=36)
boxplot(Gtemp~cycle(Gtemp),xlab="Date", ylab = "Average global temperature")

win.graph(width=5,height =6,pointsize = 10)
par(mfcol=c(2,1))
acf(diff(Gtemp),lag=36)
pacf(diff(Gtemp),lag=36)
acf(m1$residuals,main="residual",lag=36)
m1=arima(Gtemp,order=c(1,1,2),seasonal=list(order=c(0,0,1),period=24))
m1
win.graph(width=5,height =6,pointsize = 10)
tsdiag(m1,gof=36)

timex=c(1569:1688)
pre<-predict(m1,newxreg=timex,n.head=120)
plot(pre)



## trend-stationarity
time=c(1:1568) # time index
m2=lm(Gtemp~time)
summary(m2)
win.graph(width=5,height =6,pointsize = 10)
par(mfcol=c(2,1))
acf(m2$residuals,main="residual",lag=36)
pacf(m2$residuals,main="residual",lag=36) 
m2=arima(Gtemp,order=c(2,0,1),xreg=time)
m2
#tsdiag(m2,gof=36)  # Significant ACF at lag 24.
win.graph(width=5,height =5,pointsize = 10)
acf(m2$residuals,main="residual",lag=36)
m2=arima(Gt,order=c(2,0,1),seasonal=list(order=c(0,0,1),period=24),xreg=time)
m2
win.graph(width=5,height =5,pointsize = 10)
tsdiag(m2,gof=36) # model checking
timex=c(1569:1688)
pre<-predict(m2,newxreg=timex,n.head=120)
plot(pre)
### Comparison 
source("D:/teaching/TSA2019UG/Rcode2019/backtest.R")
pm1=backtest(m1,Gtemp,1368,1)
time=as.matrix(time)
pm2=backtest(m2,Gtemp,1368,1,xre=time)
###
Gt=scan(file='m-GLBTs.txt')
time=c(1:1568)
time1=c(rep(0,1212),time[1213:1568])
mm1=lm(Gt~time+time1)
summary(mm1)
x1=cbind(time,time1)
mm1=arima(Gt,order=c(2,0,1),seasonal=list(order=c(0,0,1),period=24),xreg=x1)
mm1
tsdiag(mm1,gof=36)
Box.test(mm1$residuals,lag=8,type='Ljung')
####
da=read.table("m-ncdc-noaa-glbtemp.txt")
head(da)
tail(da)
da=da[1:1568,]
temp=da[,3]
m3=arima(temp,order=c(1,1,2),seasonal=list(order=c(0,0,1),period=24))
m3
tsdiag(m3,gof=36)
m4=arima(temp,order=c(2,0,1),seasonal=list(order=c(0,0,1),period=24),xreg=time)
m4
m4$coef
sqrt(diag(m4$var.coef))
m4$coef/sqrt(diag(m4$var.coef))  # Compute t-ratios 
tsdiag(m4,gof=36)
### Backtesting 
pm3=backtest(m3,temp,1368,1)
pm4=backtest(m4,temp,1368,1,xre=time)