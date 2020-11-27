
########### monthly log stock returns of Intel Corporation from 1973.01---2009.12
setwd("C:/Users/陈方舟/Desktop/时间序列分析/chap10volatility")
da=read.table("m-intcsp7309.txt",header=T)
head(da)
intc=log(da$intc+1)### transfer simple returns to log-returns
rtn=ts(intc,frequency=12,start=c(1973,1))
plot(rtn,type='l',xlab='year',ylab='ln-rtn') # time plot
t.test(intc)  # testing the mean of returns
Box.test(intc,lag=12,type='Ljung')
win.graph(width=6,height =6)
par(mfcol=c(3,1))
acf(intc,lag=24) # ACF plots
acf(intc^2.,lag=24)
acf(abs(intc),lag=24) 


## ARCH test
y=intc-mean(intc)
Box.test(y^2,lag=12,type='Ljung')
source("C:/Users/陈方舟/Desktop/时间序列分析/chap10volatility/archTest.r")  # R script available on the book web site.
archTest(y,12)   # output edited.

###########the daily exchange rate between US dollar and Euro from 1999/01/04-2010/08/20


######## Illustrate ARCH modeling by log returns of Intel stock

library(fGarch) # Load package 

m1=garchFit(~1+garch(3,0),data=intc,trace=F) # Fit an ARCH(3) model
summary(m1)
m2=garchFit(~1+garch(1,0),data=intc,trace=F)
summary(m2)

resi=residuals(m2,standardize=T)
tdx=c(1:444)/12+1973
win.graph(width=6,height =6)
par(mfcol=c(3,1))
plot(tdx,resi,xlab='year',ylab='stand-resi',type='l')
acf(resi,lag=20)
pacf(resi^2,lag=20) 
plot(m2)
#### Student t innovations
m3=garchFit(~1+garch(1,0),data=intc,trace=F,cond.dist="std")
summary(m3)
##
  # Obtain volatility
resi=residuals(m4,standardize=T) # Standardized residuals
vol=ts(v1,frequency=12,start=c(1973,1))
res=ts(resi,frequency=12,start=c(1973,1))
win.graph(width=6,height =6)
par(mfcol=c(2,1))  # Show volatility and residuals
plot(vol,xlab='year',ylab='volatility',type='l')
plot(res,xlab='year',ylab='st. resi',type='l') 
par(mfcol=c(2,2)) # Obtain ACF & PACF
acf(resi,lag=24)
pacf(resi,lag=24)
acf(resi^2,lag=24)
pacf(resi^2,lag=24) 
# Obtain plot of predictive intervals
par(mfcol=c(1,1))
upp=0.0113+2*v1
low=0.0113-2*v1
tdx=c(1:444)/12+1973
plot(tdx,intc,xlab='year',ylab='series',type='l',ylim=c(-0.6,0.6))
lines(tdx,upp,lty=2,col='red')
lines(tdx,low,lty=2,col='red')
abline(h=c(0.0113)) ############# unconditional volatility
# Student-t innovations
m5=garchFit(~1+garch(1,1),data=intc,trace=F,cond.dist="std")
summary(m5)
v2=volatility(m5)
m6=garchFit(~1+garch(1,1),data=intc,trace=F,cond.dist='sstd')
summary(m6)
v3=volatility(m6)
par(mfcol=c(3,1))
plot(tdx,v1,xlab='year',ylab='volatility',type='l',ylim=c(0.06,0.3))
title(main='(a) Gaussian')
plot(tdx,v2,xlab='year',ylab='volatility',type='l',ylim=c(0.06,0.3))
title(main='(b) Student-t')
plot(tdx,v3,xlab='year',ylab='volatility',type='l',ylim=c(0.06,0.3))
title(main='(c) Skew Student-t') 
cor(cbind(v1,v2,v3))


library(fBasics)
basicStats(intc)
tt=-0.5526/sqrt(6/444) # Testing skewness of the data
tt
tt=(0.8717-1)/0.0629 # Testing skewness of the model.
tt
pv=2*pnorm(tt)  # Compute p-value 
pv
plot(m6)
### two-pass analysis
yt=intc-mean(intc)
m1=arima(yt^2,order=c(1,0,1))
m1
mean(intc)
fit=yt^2-m1$residuals
v3=volatility(m6)  # m6 is GARCH(1,1) with skew-t innovations.
cor(v3,sqrt(fit))
####



source("D:/teaching/Tsayreference/chart4/Igarch.R")
m7=Igarch(intc)
names(m7)
###
y=intc*100   # Intel stock returns in percentages
source("D:/teaching/Tsayreference/chart4/garchM.R")  # Compile the script
m8=garchM(y)



##
m9=garchFit(~1+aparch(1,1),data=y,trace=F)
summary(m9)
m10=garchFit(~1+aparch(1,1),data=y,delta=2,include.delta=F,trace=F)
summary(m10)
plot(m10)
###

source("D:/teaching/Tsayreference/chart4/Ngarch.R")
m11=Ngarch(y)
res=m11$residuals
vol=m11$volatility
resi=res/vol
Box.test(resi,lag=10,type='Ljung')
Box.test(resi^2,lag=10,type='Ljung')
####################################
