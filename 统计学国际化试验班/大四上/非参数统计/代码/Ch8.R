###
##### R code for Chapter 8: Bootstrap Methods
###

library(boot)

###
##### Example 8.1.1. Latch failure time, calculate bootstrap MSE, bias and variance
###
x = c(7, 11, 15, 16, 20, 22, 24, 25, 29, 33, 34, 37, 41, 42, 49, 57, 66, 71, 84, 90)

### view the data set  (optional)

hist(x,prob=T)
lines(density(x), col="red")

### Define my own function of median estimator
mymedian.func = function(x, d)
{
	#d is the index of resampled observations
	median(x[d])
	}

### bootstrap median estimates based on R bootstrap samples
median.boot = boot(x, mymedian.func, R=1000)


### outputs of boot()
### t0: The observed value of statistic applied to data. That is, \hat\theta
### t: A matrix with R rows each of which is a bootstrap replicate of statistic. That is, \hat\theta^*


### median estimate based on the observed sample
theta.obs = median(x) #equivalent to mymedian.func(x, d=1:length(x))
theta.obs
median.boot$t0

### bootstrap bias
bias1=mean(median.boot$t) - median.boot$t0

### bootstrap variance
var1=var(median.boot$t)

### bootstrap MSE
mean((median.boot$t-median.boot$t0)^2)   # close to bias1^2+var1

### the boostrap bias and se are summarized in the outputs of boot()
median.boot

###
##### You may also obtain bootstrap estimators directly without using boot() function
###
B = 1000
theta.boot = rep(0, B)
n = length(x)
for(b in 1:B)
{
    # resampled subjects indices
    ind = sample(1:n, replace=TRUE)
    # bootstrap sample
    x.boot = x[ind]
    theta.boot[b] = median(x.boot)
}
### bootstrap bias
mean(theta.boot) - theta.obs

### bootstrap variance
var(theta.boot)

### bootstrap MSE
mean((theta.boot-theta.obs)^2)

### the slight difference from the outputs of boot() is due to different seeds



###
##### Example 8.2.1 Construct 95% CI of the mean
###
x = c(7, 11, 15, 16, 20, 22, 24, 25, 29, 33, 34, 37, 41, 42, 49, 57, 66, 71, 84, 90)

### Define my own function of median estimator
mymean.func = function(x, d)
{
	#d is the index of resampled observations
	mean(x[d])
	}

mean.boot = boot(x, mymean.func, R=1000)
bias=mean(mean.boot$t)-mean.boot$t0
sd=sqrt(var(mean.boot$t))



### bootstrap CI of the mean estimator

boot.ci(mean.boot)




mean.boot$t0-bias-1.96*sd
mean.boot$t0-bias+1.96*sd

quantile(mean.boot$t,c(0.025, 0.975))
2*mean.boot$t0-quantile(mean.boot$t,c(0.975, 0.025)) 


# t-pivotal 
x0=mean(x)
n=length(x)

mymean2.func=function(x,d){
sqrt(n)*(mean(x[d])-x0)/sd(x[d])
}


boot.t=boot(x, mymean2.func, R=1000)
qr=quantile(boot.t$t, c(0.975, 0.025))
x0-qr[1]*sd(x)/sqrt(n)
x0-qr[2]*sd(x)/sqrt(n)


### nonp bootstrap
n=200
x=rnorm(n)
y=2*x+0.5*rnorm(n)

grade2=data.frame(V1=x,V2=y)

mybeta<-function(formula,data,d){  
new<-data[d,]  
fit<-lm(formula,data=new)  
return(fit$coef[2])  
}  
result<-boot(grade2,mybeta,R=1000,formula=V2~V1)  
boot.ci(result)  



# parametric bootstrap 

m1=lm(V2~V1, data=grade2)

attach(grade2)
fulldata=data.frame(V1=V1, V2=V2, yhat=fitted(m1), res=resid(m1))
detach(grade2)

mybeta2<-function(formula,data){  
fit<-lm(formula,data=data)  
return(fit$coef[2])  
}  


my.rg2 <- function(data, mle) {
attach(data)
ind=sample(1:mle)
resn=res[ind]
yn=yhat+resn
out=data.frame(V1=V1, V2=yn)
detach(data)
out
}


result2<-boot(data=fulldata, mybeta2, R=1000,formula=V2~V1, sim = "parametric",
                      ran.gen = my.rg2, mle=n)  

# library(bootstrap)
#  mybeta3 <- function(x,xdata){ coef(lm(xdata[x,1]~xdata[x,2]))[2] }
#        results <- bcanon(1:n,1000,mybeta3,grade2)
     




# Another 

hatb=coef(m1)

my.rg <- function(data, hatb) {
attach(data)
yhat=hatb[1]+hatb[2]*V1
res=V2-yhat
n=length(res)
ind=sample(1:n)
resn=res[ind]
yn=yhat+resn
out=data.frame(V1=V1, V2=yn)
detach(data)
out
}

result3<-boot(data=grade2, mybeta2, R=1000,formula=V2~V1, sim = "parametric",
                      ran.gen = my.rg, mle = hatb)  

### 8.4.1 



x = c(7, 11, 15, 16, 20, 22, 24, 25, 29, 33, 34, 37, 41, 42, 49, 57, 66, 71, 84, 90)






###
####### Bootstrap test of the mean
###
mu0 = 30

### Idea1 (bad)
t1 = mean(x)-mu0
mean.boot = boot(x, mymean.func, R=1000)
t1.star = mean.boot$t-mu0
pval1 = mean(t1.star>=t1)
pval1

hist(mean.boot$t);
abline(v=mu0,col="blue")
abline(v=mean(x),col="red")

### Idea2 (better but not ideal)
t2 = mean(x)-mu0
mean.boot = boot(x, mymean.func, R=1000)
t2.star = mean.boot$t-mean(x)
pval2 = mean(t2.star>=t2)
pval2

### Idea3 (best)
t3 = (mean(x)-mu0)/(sd(x)/sqrt(n))

myt.func = function(x, d)
{
	n = length(x)
	mean(x[d])/(sd(x[d])/sqrt(n))
	}
	
out = boot(x-mean(x), myt.func, R=1000)
t3.star = out$t
pval3 = mean(t3.star>=t3)
pval3

####
###### Regression: HW and Final Exam
####
grade <- read.csv("st745.csv")
x <- grade[,1] # Middle term exam score
y <- grade[,2] # Final exam score


### Least squares regression results (inference based on normal assumption)
summary(lm(y~x))

### bootstrap
### define my function of obtaining slope estimate
lm.func <- function(data, d)
{
    lm(data[d,1]~data[d,-1])$coef[2]
}
slope.boot = boot(cbind(y,x), lm.func , R=1000)
slope.boot
boot.ci(slope.boot)

bhat = slope.boot$t0
bhat.boot = slope.boot$t
apply(bhat.boot, 2, sd)

## Test H0: beta1 = 0.5
## approach 2
t1 = bhat-0.5
t1.star = bhat.boot-bhat
pval1 = mean(abs(t1.star)>=abs(t1))
pval1

hist(t1.star); abline(v=t1,col=2)

## approach 3: based on the pivotal
t2 = (bhat-0.5)/0.2662

t2.star = NULL
n = length(y)
for(b in 1:1000)
{
    d = sample(1:n, n, replace=TRUE)
    coef = summary(lm(y[d]~x[d]))$coef
    se = coef[2,2]
    slope = coef[2,1]
    t2.star = c(t2.star, (slope - bhat)/se)
}

pval2 = mean(abs(t2.star)>=abs(t2))
pval2

hist(t2.star); abline(v=t2,col=2)
