da=read.table("D:/teaching/TSA2020UG/Rcode2019/q-ko-earns8309.txt")
head(da)
eps=log(da$value)
koeps=ts(eps,frequency=4, start=c(1983,1))
plot(koeps,type='l')
boxplot(koeps,frequency=4)
acf(koeps,lag=20)
par(mfcol=c(2,2))
diff_eps=diff(koeps)
sdiff_eps=diff(koeps,4)
dsdiff_eps=diff(sdiff_eps)
acf(koeps,lag=20)
acf(diff_eps,lag=20)
acf(sdiff_eps,lag=20)
acf(dsdiff_eps,lag=20)
par(mfcol=c(3,1))
plot(diff_eps,xlab='year',ylab='diff',type='l')
plot(sdiff_eps,xlab='year',ylab='sdiff',type='l')
plot(dsdiff_eps,xlab='year',ylab='dsdiff',type='l')
%%%%%%%%%%%%%%%%%%
  Estimation
%%%%%%%%%%%%%%%%
m1=arima(koeps,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=4))
tsdiag(m1,gof=20)
Box.test(m1$residuals,lag=12,type='Ljung')
pp=1-pchisq(13.30,10)
length(koeps)
y=koeps[1:100]
m1=arima(y,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=4))

%%%%%%% Prediction %%%%%%%%%%%%%%%
pm1=predict(m1,7)
names(pm1)
pred=pm1$pred
se=pm1$se

%%%%%%%% Anti-log transformation
ko=da$value
fore=exp(pred+se^2/2)
v1=exp(2*pred+se^2)*(exp(se^2)-1)
s1=sqrt(v1)
eps=ko[80:107]
tdx=(c(1:28)+3)/4+2002
upp=c(ko[100], fore+2*s1)
low=c(ko[100], fore-2*s1)
plot(tdx,eps,xlab='year',ylab='earnings',type='l',ylim=c(0.35,1.3))
points(tdx[22:28],fore,pch='*')
lines(tdx[21:28],upp,lty=2)
lines(tdx[21:28],low,lty=2)
points(tdx[22:28],ko[101:107],pch='o',cex=0.7)
