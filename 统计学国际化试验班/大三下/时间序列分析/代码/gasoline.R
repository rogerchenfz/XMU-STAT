library(TSA)
library(fUnitRoots)

rm(list=ls())

da=read.table("datasets/w-petroprice.txt",header=T)
da1=read.table("datasets/w-gasoline.txt")

pgs=log(da1[,1])
pus=log(da$US)

tdx=c(1:717)/52+1997  # calendar time

par(mfcol=c(2,1))
plot(tdx,pgs,xlab='year',ylab='ln(price)',type='l')
title(main='(a) Gasoline')
plot(tdx,pus,xlab='year',ylab='ln(price)',type='l')
title(main='(b) Crude oil')
#
dpgs=diff(pgs)
plot(tdx[-1],dpgs,main="gasoline",xlab='year',ylab='Growth rate',type='l')
#
win.graph(width=6,height =6,pointsize = 10)
par(mfcol=c(2,1))
acf(dpgs,main="diff(pgs)",lag=20)
pacf(dpgs,main="diff(pgs)",lag=20)

m1=ar(diff(pgs),method='mle')
m1$order

t.test(dpgs)#检验均值是否为0

m1=arima(dpgs,order=c(5,0,0),include.mean=F)
m1#系数显著性，4阶系数不显著
#summary(m1)

m1=arima(dpgs,order=c(5,0,0),include.mean=F,fixed=c(NA,NA,NA,0,NA))
m1

#m1=arima(dpgs,order=c(5,0,0),include.mean=F,fixed=c(NA,0,NA,0,NA))
#m1

#m1=arima(dpgs,order=c(0,0,5),include.mean=F) MA(5)
#m1
#sigma^2 estimated as 0.000327:  log likelihood = 1856.89,  aic = -3701.77

win.graph(width=5,height =6,pointsize = 10)
tsdiag(m1,gof=20)#diagnostic checking

dpus=diff(pus)
m3=lm(dpgs~-1+dpus)#lm
summary(m3)
win.graph(width=5,height =6,pointsize = 10)
par(mfcol=c(2,1))
acf(m3$residuals,main="Residual of model (3)",lag=20)
pacf(m3$residuals,main="Residual of model (3)",lag=20)

m4=ar(m3$residuals,method='mle')
m4$order
m4=arima(dpgs,order=c(6,0,0),include.mean=F,xreg=dpus)
m4
m4=arima(dpgs,order=c(5,0,0),include.mean=F,xreg=dpus)
m4
m5=arima(dpgs,order=c(6,0,0),include.mean=F,xreg=dpus,fixed=c(NA,NA,NA,0,NA,0,NA))
m5
win.graph(width=5,height =6,pointsize = 10)
tsdiag(m5,gof=20)

m6=lm(dpgs[-1]~-1+dpus[-length(dpus)])
summary(m6)
win.graph(width=5,height =6,pointsize = 10)
par(mfcol=c(2,1))
acf(m6$residuals,main="Residual of model (5)",lag=20)
pacf(m6$residuals,main="Residual of model (5)",lag=20)

c1=c(NA,NA,NA,0,NA)
pm1=backtest(m1,dpgs,316,1,fixed=c1,inc.mean=F)
c4=c(NA,NA,NA,0,NA,NA)
pm4=backtest(m4,dpgs,316,1,xre=dpus,inc.mean=F,fixed=c4)
tdx=tdx[2:717]
pm4fit=dpgs[317:716]-pm4$error
pm1fit=dpgs[317:716]-pm1$error
plot(tdx[317:716],dpgs[317:716],xlab='year',ylab='growth',type='l')
points(tdx[317:716],pm1fit,pch='*')
plot(tdx[317:716],dpgs[317:716],xlab='year',ylab='growth',type='l')
points(tdx[317:716],pm4fit,pch='*')
m6=lm(dpgs[2:716]~-1+dpus[1:715])
summary(m6)
acf(m6$residuals,lag=20)
pacf(m6$residuals,lag=20)
m7=ar(m6$residuals,method='mle')
m7$order
m7=arima(dpgs[2:716],order=c(9,0,0),include.mean=F,xreg=dpus[1:715])
m7
m7=arima(dpgs[2:716],order=c(9,0,0),include.mean=F,xreg=dpus[1:715],fixed=c(NA,NA,NA,0,NA,0,0,0,NA,NA))
m7
tsdiag(m7,gof=20)
c7=c(NA,NA,NA,0,NA,0,0,0,NA,NA)
pm7=backtest(m7,dpgs[2:716],315,1,xre=dpus[1:715],inc.mean=F,fixed=c7)

