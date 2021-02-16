source("functions-Ch5.R")

		#####################################
		#  Permutation test for rho/beta1   #
		#####################################
x=c(0,96,65,58,56)
y=c(0,166,130,118,130)


### direct function for calculating the sample correlation: r
(r = cor(x,y))

### least squares coefficient estimates
(lm(y~x))

#### calculate r and beta1hat manually:
(sumx = sum(x))
(sumy = sum(y))
(sumxy = sum(x*y))
(sumxx = sum(x^2))
(sumyy = sum(y^2))

(n = length(x))
xbar = sumx/n
ybar = sumy/n
Sxy = 1/(n-1)*(sumxy - n*xbar*ybar)
Sx2 = 1/(n-1)*(sumxx - n*xbar^2)
Sy2 = 1/(n-1)*(sumyy - n*ybar^2)
Sx = sqrt(Sx2)
Sy = sqrt(Sy2)

(robs = Sxy/(Sx*Sy))
(bhat1.obs = Sxy/(Sx2))

### Test H0:r=0 versu Ha:r>0

##
##### t-test assuming normality
##
tt = sqrt((n-2)/(1-robs^2))*robs
1-pt(tt,n-2)

##
##### large sample approximation
##
(Zr = sqrt(n-1)*r)
1-pnorm(Zr)

permr <- perm.approx.r(x,y,1000)
mean(permr >= robs)
mean(abs(permr) >= abs(robs))

		#####################################
		#  Spearman Correlation (ties)      #
		#####################################
##
##### Reading ability data set
##
x=1:10  #rankings of 10 children on a reading test
y = c(3,2,1,4,5,6,8,7,10,9) #teacher's ranking on the reading ability

(x = rank(x))  ###no need for this data as they are already ranked
(y = rank(y))

## Spearman correlation
(rs.obs = cor(x, y))

## permutation test for the Spearman correlation
permr <- perm.approx.r(x, y, 1000)
mean(permr >= rs.obs)
mean(abs(permr) >= abs(rs.obs))



##### data set: scores of ten projects at a Science Fair (Table 5.2.2 of Higgins)
x = c(8,8,7,8,5,6,6,9,8,7) #Judge A
y = c(7,8,8,5,6,4,5,8,6,9) #Judge B
x
y
(x = rank(x))
(y = rank(y))

## Spearman correlation
(rs.obs = cor(x, y))

## permutation test for the Spearman correlation
permr <- perm.approx.r(x, y, 1000)
mean(permr >= rs.obs)
mean(abs(permr) >= abs(rs.obs))
		#######################
		#  Kendall's tau      #
		#######################
#A subset of the ST372 grades:
x=c(0,96,65,58,56)
y=c(0,166,130,118,130)

tauobs <- getTau(x,y)
tauobs

permtau <- perm.approx.tau(x,y,1000)
mean(permtau >= tauobs)


			###############################
			#  permutation chi-sq test    #
			###############################
##
##### example 1: Satisfaction v.s. Treatment
##
satisfy = c(1,1,2,2,2,3,3)
gender = c(1,1,1,1,2,2,2) #1: female; 2:male

(partytable <- table(gender, satisfy))  #convert from data to table

tabletodata(partytable)  #convert from table to data

f <- partytable
fi. <- matrix(apply(f,1,sum),ncol=1)
f.j <- matrix(apply(f,2,sum),nrow=1)
e <- (fi. %*% f.j)/sum(fi.)
X2p <- sum((f-e)^2/e)

chisq.test(satisfy, gender)
# or chisq.test(gender, satisfy)
# or chisq.test(table(satisfy, gender))

sweep(f,1,fi.,"/")  # table of row percentages
sweep(f,2,f.j,"/")  # table of column percentages

permX2 <- permapproxX2(satisfy,gender,1000)
mean(permX2 >= chisq.test(satisfy,gender)$statistic)

## compare with
chisq.test(satisfy, gender)$p.value

##
###### example2: Gender versus Party
##
party <- c(rep(1,279),rep(2,73),rep(3,225),
           rep(1,165),rep(2,47),rep(3,191))#1:Dem, 2: Ind, 3:Rep
gender <- c(rep(1,279+73+225),rep(2,165+47+191))  #1:female, 2:male

(partytable <- table(gender, party))

chisq.test(party, gender)

f <- partytable
fi. <- matrix(apply(f,1,sum),ncol=1)
f.j <- matrix(apply(f,2,sum),nrow=1)
e <- (fi. %*% f.j)/sum(fi.)
X2p <- sum((f-e)^2/e)

sweep(f,1,fi.,"/")  # table of row percentages
sweep(f,2,f.j,"/")  # table of column percentages

## Permutation test for contingency table

permX2 <- permapproxX2(party, gender,1000)
mean(permX2 >= chisq.test(party, gender)$statistic)

## compare with
chisq.test(party, gender)$p.value


##
###### example3: Presidential preference
##
president <- c(rep(1,10),rep(2,3))
region <- c(rep(1,6),rep(2,4),rep(1,1),rep(2,2))

(partytable <- table(president, region))

chisq.test(president, region)

f <- partytable
fi. <- matrix(apply(f,1,sum),ncol=1)
f.j <- matrix(apply(f,2,sum),nrow=1)
e <- (fi. %*% f.j)/sum(fi.)
X2p <- sum((f-e)^2/e)

sweep(f,1,fi.,"/")  # table of row percentages
sweep(f,2,f.j,"/")  # table of column percentages

## Permutation test for contingency table

permX2 <- permapproxX2(president, region,1000)
mean(permX2 >= chisq.test(president, region)$statistic)

## compare with
chisq.test(president, region)$p.value


##
###### Fisher's exact test
##

## Example 1: Tea tasting
TeaTasting <-
matrix(c(3, 1, 1, 3),
       nrow = 2,
       dimnames = list(Guess = c("Milk", "Tea"),
                       Truth = c("Milk", "Tea")))
                       
fisher.test(TeaTasting, alternative = "greater")
## => p=0.2429, association could not be established

fisher.test(TeaTasting)

## Example 2: presidental preference

presidental <- matrix(c(6,4,1,2),nrow=2, dimnames = list(state = c("West", "East"),
                       candidate = c("Bush", "Kerry")))

fisher.test(presidental, alternative = "greater")
## => p=0.4371, association could not be established

## hypergeometric distribution
m = 7
n = 6
k = 10
dhyper(6:10, m, n, k)

sum(dhyper(6:10, m, n, k))










######## Additional example: Kendall's tau
#package Kendall
M = c(70, 69, 65, 64, 66, 65, 64, 66, 60, 70, 66)
D = c(67,64,62,64,69,70,65,66,63,74,60)
x=M
y=D
tauobs <- getTau(x,y)
tauobs

 y[5]=66
 y[8]=69
 getTau(x,y)

library(Kendall)
Kendall(x,y)
getTau(x,y)
