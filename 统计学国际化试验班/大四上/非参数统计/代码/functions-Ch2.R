### functions used in chapter 2



# the combinations function from package gregmisc
combinations <- function (n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE) 
{
    if (mode(n) != "numeric" || length(n) != 1 || n < 1 || (n%%1) != 
        0) 
        stop("bad value of n")
    if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r%%1) != 
        0) 
        stop("bad value of r")
    if (!is.atomic(v) || length(v) < n) 
        stop("v is either non-atomic or too short")
    if ((r > n) & repeats.allowed == FALSE) 
        stop("r > n and repeats.allowed=FALSE")
    if (set) {
        v <- unique(sort(v))
        if (length(v) < n) 
            stop("too few different elements")
    }
    v0 <- vector(mode(v), 0)
    if (repeats.allowed) 
        sub <- function(n, r, v) {
            if (r == 0) 
                v0
            else if (r == 1) 
                matrix(v, n, 1)
            else if (n == 1) 
                matrix(v, 1, r)
            else rbind(cbind(v[1], Recall(n, r - 1, v)), Recall(n - 
                1, r, v[-1]))
        }
    else sub <- function(n, r, v) {
        if (r == 0) 
            v0
        else if (r == 1) 
            matrix(v, n, 1)
        else if (r == n) 
            matrix(v, 1, n)
        else rbind(cbind(v[1], Recall(n - 1, r - 1, v[-1])), 
            Recall(n - 1, r, v[-1]))
    }
    sub(n, r, v[1:n])
}

##################################################################################################
### A function for obtaining the permutation distribution of difference of means (E(X)-E(Y))   ###
### (can be adapted for comparing trimmed means, medians, variances, etc.)                     ###
##################################################################################################
 
perm.2sample = function(x, y, alternative = c("two.sided", "less", "greater"), stat=c("meandiff", "mediandiff", "trmdiff", "sumX"), trim=0)
{
     # exact permutation
    # test $H_0: delta=mu_x-mu_y=0$ versus $H_a: delta>0 (<0) (\neq 0)
    # stat="meandiff": difference of means
    # stat="mediandiff": difference of medians
    # stat="trmdiff": difference of trimmed means 
    # trim: the fraction (0 to 0.5) of observations to be trimmed from each end of the data set before the mean is computed. 
    #stat="sumX": sum of observations from the x group (equivalent with meandiff)

    m = length(x)
    n = length(y)
    N = m+n
    nperm = choose(N, m)
    if(nperm>5000){ warning("The number of permutations is over 5000, consider random permutation!")
        stop()
        }
    
    xy = c(x, y)
    idx = combinations(n=N, r=m)
    R = nrow(idx)
    permut = NULL
   for(i in 1:R)
    {
        permut = rbind(permut, c(xy[idx[i,]], xy[-idx[i,]]))
    }
    if(stat=="meandiff") trim=0
    if(stat %in% c("meandiff", "trmdiff"))
        D = apply(permut, 1, function(x) mean(x[1:m], trim)-mean(x[-(1:m)], trim))
    if(stat=="mediandiff")
        D = apply(permut, 1, function(x) median(x[1:m])-median(x[-(1:m)]))
        
    if(stat=="sumX")
        D = apply(permut, 1, function(x) sum(x[1:m]))
    
    if(alternative=="greater")
        pval = mean(D >= D[1])
    if(alternative=="less")
        pval = mean(D <= D[1])
     if(alternative=="two.sided") 
        pval = mean(abs(D) >= abs(D[1]))
            
     return(list(pval=pval, Dobs=D[1]))
    }



            ####################################################
            ##    Random Sampling two-sample Permutation      ##
            ####################################################
rand.perm = function(x, y, R, alternative = c("two.sided", "less", "greater"), stat=c("meandiff", "mediandiff", "trmdiff", "sumX"), trim=0)
{
    #stat="meandiff": mean difference between two groups
    #stat="mediandiff": median difference between two groups
    #stat="trmdiff: trimmed mean difference between two groups, trim=0.1 means 10% of the data will be trimmed.
    #stat="sumX": sum of observations from the x group
    
    m = length(x)
    n = length(y)
    N = m+n
    xy = c(x, y)
    permut = NULL
  for(i in 1:R)
    {
        idx = sample(1:N, replace=FALSE)
        permut = rbind(permut, c(xy[idx[1:m]], xy[idx[-(1:m)]]))
    }
    if(stat=="meandiff") trim=0
    if(stat %in% c("meandiff", "trmdiff"))
        {
            D = apply(permut, 1, function(x) -mean(x[-(1:m)], trim)+mean(x[1:m], trim))
            Dobs = mean(x, trim) - mean(y, trim)
            }
    if(stat=="sumX")
        {
            D = apply(permut, 1, function(x) sum(x[1:m]))
            Dobs = sum(x)
        }
    if(stat=="mediandiff")
        {
            D = apply(permut, 1, function(x) -median(x[-(1:m)])+median(x[1:m]))
            Dobs = -median(y) + median(x)
        }
    if(alternative=="greater")
        pval = mean(D >= Dobs)
    if(alternative=="less")
        pval = mean(D <= Dobs)
     if(alternative=="two.sided") 
        pval = mean(abs(D) >= abs(Dobs))
    
    hist(D, main="Hist of D under the null")
    abline(v=Dobs, col=2)
    
     return(list(pval=pval, Dobs=Dobs))
    }


####
######## A function for obtaining the permutation distribution of sum(x[1:n1]) 
###
perm.dist.sum <- function(x, n1)
{
    #the permutation of sum of x[1:n1]
    #n1: the number of observations in group 1
    #x: combined data
    #output: the sums of all possible permutations
    n = length(x)
    nperm = choose(n, n1)
    if(nperm>5000){ warning("The number of permutations is over 5000, consider random permutation!")
        stop()
        }
    temp <- combinations(length(x), n1)
    results <- rep(NA, nrow(temp))
    for(i in 1:nrow(temp))
        results[i] <- sum(x[temp[i,]])
    results
}
####
######## A function for obtaining the random permutation distribution of sum(x[1:n1]) 
###
rand.perm.dist.sum <- function(x, n1, R)
{
    #the permutation of sum of x[1:n1]
    #n1: the number of observations in group 1
    #x: combined data
    #R: number of permutations
    #output: the sums of all possible permutations
    N = length(x)
    results <- rep(NA, R)
    for(i in 1:R)
    {
         idx = sample(1:N, n1, replace=FALSE)
        results[i] <- sum(x[idx])
     }
    results
}


###
######## Function for calculating the siegelranks
###
siegelrank <- function(x,ansari=F)
{
    # if ansari=T, return the Ansari-Bradley's ranks
    if (ansari==T) ansval <- .5
    n <- length(x)
    mat <- cbind(x,1:n)
    mat <- mat[order(mat[,1]),1:2]
    i <- 1
    rx <- rep(NA,n)
    while (i <= n)
    {
        if (i %% 4 == 1) rx[(i %/% 4)*2 + 1] <- i 
        if (i %% 4 == 2) rx[n - (i %/% 4)*2] <- i 
        if (i %% 4 == 3) rx[n - (i %/% 4)*2 - 1] <- i 
        if (i %% 4 == 0) rx[(i %/% 4)*2 ] <- i 
        i <- i+1
    }
    if (ansari==T) rx <- (rx + rx[n:1])/2
    mat <- cbind(mat,rx)
    mat <- mat[order(mat[,2]),1:3]
    mat[,3]
}


########### Function for carrying out siegel-tukey test of equality of scales
########### The function can adjust for tied observations
siegel.tukey=function(x, y,
    id.col=FALSE, adjust.median=F, rnd=-1,alternative="two.sided",mu=0,paired=FALSE,
    exact=FALSE,correct=TRUE,conf.int=FALSE,conf.level=0.95)
{
##see details at    http://www.r-bloggers.com/siegel-tukey-a-non-parametric-test-for-equality-in-variability-r-code/
#   Usage:
#
#siegel.tukey(x,y,id.col=FALSE,adjust.median=FALSE,rnd=8, ...)
#Arguments:
#x: a vector of data
#y: Data of the second group (if id.col=FALSE) or group indicator (if id.col=TRUE). In the latter case, y MUST take 1 or 2 to indicate observations of group 1 and 2, respectively, and x must contain the data for both groups.
#id.col: If FALSE (default), then x and y are the data columns for group 1 and 2, respectively. If TRUE, the y is the group indicator.
#adjust.median: Should between-group differences in medians be leveled before performing the test? In certain cases, the Siegel-Tukey test is susceptible to median differences and may indicate significant differences in variability that, in reality, stem from differences in medians.
#rnd: Should the data be rounded and, if so, to which decimal? The default (-1) uses the data as is. Otherwise, rnd must be a non-negative integer. Typically, this option is not needed. However, occasionally, differences in the precision with which certain functions return values cause the merging of two data frames to fail within the siegel.tukey function. Only then rounding is necessary. This operation should not be performed if it affects the ranks of observations.
#??? arguments passed on to the Wilcoxon test. See ?wilcox.test
#Value: Among other output, the function returns rank sums for the two?groups, the associated Wilcoxon's W, 
# and the p-value for a Wilcoxon test on?tie-adjusted Siegel-Tukey ranks (i.e., it performs and returns 
#Siegel-Tukey test). If significant, the group with the smaller rank sum has greater variability.
#References: Sidney Siegel and John Wilder Tukey (1960) “A nonparametric sum of ranks procedure for relative spread in unpaired samples.??? Journal of the American Statistical Association. See also, David J. Sheskin (2004) ”Handbook of parametric and nonparametric statistical procedures.??? 3rd edition. Chapman and Hall/CRC. Boca Raton, FL.
#Notes: The Siegel-Tukey test has relatively low power and may, under certain conditions, indicate significance due to differences in medians rather than differences in variabilities (consider using the argument adjust.median).
#Output (in this order)
#1. Group medians
#2. Wilcoxon-test for between-group differences in median (after the median adjustment if specified)
#3. Unique values of x and their tie-adjusted Siegel-Tukey ranks
#4. Xs of group 1 and their tie-adjusted Siegel-Tukey ranks
#5. Xs of group 2 and their tie-adjusted Siegel-Tukey ranks
#6. Siegel-Tukey test (Wilcoxon test on tie-adjusted Siegel-Tukey ranks)

#id.col=FALSE; adjust.median=F; rnd=-1;alternative="two.sided";mu=0,paired=FALSE;
#   exact=FALSE;correct=TRUE;conf.int=FALSE;conf.level=0.95
    

 if(id.col==FALSE){
   data=data.frame(c(x,y),rep(c(1,2),c(length(x),length(y))))
   } 
   else {
    data=data.frame(x,y)
   }
 names(data)=c("x","y")
 data=data[order(data$x),]
 if(rnd> -1){
    data$x = round(data$x,rnd)}
 if(adjust.median==T){
    data$x[data$y==1]=data$x[data$y==1]-(median(data$x[data$y==1])-median(data$x[data$y==2]))/2
    data$x[data$y==2]=data$x[data$y==2]-(median(data$x[data$y==2])-median(data$x[data$y==1]))/2
 }
 cat("Median of group 1 = ",median(data$x[data$y==1]),"\n")
 cat("Median of group 2 = ",median(data$x[data$y==2]),"\n")
 cat("Test of median differences","\n")
 print(wilcox.test(data$x[data$y==1], data$x[data$y==2]))
 
a=rep(seq(ceiling(length(data$x)/4)),each=2)
 b=rep(c(0,1),ceiling(length(data$x)/4))
 rk.up=c(1,(a*4+b))[1:ceiling(length(data$x)/2)]
 rk.down=rev(c(a*4+b-2)[1:floor(length(data$x)/2)])
cat("Performing Siegel-Tukey rank transformation...","\n")
 rks=c(rk.up,rk.down) #Siegel ranks without accounting for ties
 
 unqs=unique(sort(data$x))
 
# correct Siegel ranks to account for ties
# by averaging the ranks of tied observations.
 corr.rks=tapply(rks,data$x,mean)
 
 cbind(unqs,corr.rks)
 rks.data=data.frame(unqs,corr.rks)
 names(rks.data)=c("unique values of x","tie-adjusted Siegel-Tukey rank")
 print(rks.data,row.names=F)
 names(rks.data)=c("unqs","corr.rks")
 data=merge(data,rks.data,by.x="x",by.y="unqs")
 rk1=data$corr.rks[data$y==1]
 rk2=data$corr.rks[data$y==2]
 cat("\n","Tie-adjusted Siegel-Tukey ranks of group 1","\n")
 group1=data.frame(data$x[data$y==1],rk1)
 names(group1)=c("x","rank")
 print(group1,row.names=F)
 cat("\n","Tie-adjusted Siegel-Tukey ranks of group 2","\n")
 group2=data.frame(data$x[data$y==2],rk2)
 names(group2)=c("x","rank")
 print(group2,row.names=F)
 cat("\n")
 cat("Siegel-Tukey test","\n")
 cat("Siegel-Tukey rank transformation performed.","Tie adjusted ranks computed.","\n")
 if(adjust.median==T) {cat("Medians adjusted to equality.","\n")} else {cat("Medians not adjusted.","\n")}
 cat("Rank sum of group 1 =", sum(rk1),"    Rank sum of group 2 =",sum(rk2),"\n")

# perform Wilcoxons rank sum test based on ranks rk1, rk2 
print(wilcox.test(rk1,rk2,alternative=alternative,mu=mu,paired=paired,exact=exact,correct=correct,conf.int=conf.int,conf.level=conf.level))
# the returned  W is the Mann Whiteney's test statistic, so W_X should be W+n(n+1)/2
}


            ########################
            #  Tests on Deviances  #
            ########################
######################  functions begin here for RMD test ##############
### permutation distribution of rmd2
perm.dist.rmd2 <- function(x, n1)
{
    #2-sided alternative
    #n1: the number of observations in group 1
    #x: combined data
    #output: the sums of all possible permutations the RMD test statistic (for two-sided alternatives)
    temp <- combinations(length(x), n1)
    results <- rep(NA, nrow(temp))
    for(i in 1:nrow(temp))
    {
        idx = temp[i,]
        results[i] <- max( mean(abs(x[idx])), mean(abs(x[-idx])))/
                      min( mean(abs(x[idx])), mean(abs(x[-idx])))
        }
        return(results)
}

### approximate permutation of rmd2
perm.approx.rmd2 <- function(x, n1, R=1000)
{
    #2-sided alternative
    #n1: the number of observations in group 1
    #x: combined data
    #output: the sums of R permutations of the RMD test statistic (for two-sided alternatives)
    results <- rep(NA, R)
    N = length(x)
    for(i in 1:R)
    {
        idx = sample(1:N, n1)
        results[i] <- max( mean(abs(x[idx])), mean(abs(x[-idx])))/
                      min( mean(abs(x[idx])), mean(abs(x[-idx])))
        }
    return(results)
}

### permutation distribution of rmd
perm.dist.rmd <- function(x, n1)
{
    #1-sided alternative
    #n1: the number of observations in group 1
    #x: combined data
    #output: the sums of all possible permutations the RMD test statistic (for two-sided alternatives)
    temp <- combinations(length(x), n1)
    results <- rep(NA, nrow(temp))
    for(i in 1:nrow(temp))
    {
        idx = temp[i,]
        results[i] <- mean(abs(x[idx])) / mean(abs(x[-idx]))
        }
        return(results)
}

### approximate permutation of rmd
perm.approx.rmd <- function(x, n1, R=1000)
{
    #1-sided alternative
    #n1: the number of observations in group 1
    #x: combined data
    #output: the sums of R permutations of the RMD test statistic (for two-sided alternatives)
    results <- rep(NA, R)
    N = length(x)
    for(i in 1:R)
    {
        idx = sample(1:N, n1)
        results[i] <- mean(abs(x[idx])) / mean(abs(x[-idx]))
        }
    return(results)
}

################### functions end here for RMD test ################


            ############################
            #  Kolmogrov-Smirnov test  #
            ############################

### permutation distribution of Kolmogrov-Smirnov test
perm.dist.ks <- function(x, n1)
{
    #n1: the number of observations in group 1
    #x: combined data
    #output: the sums of all possible permutations
    n = length(x)
    nperm = choose(n, n1)
    if(nperm>5000){ warning("The number of permutations is over 5000, consider random permutation!")
        stop()
        }
    temp <- combinations(length(x), n1)
    results <- rep(NA, nrow(temp))
    for(i in 1:nrow(temp))
    {
        idx = temp[i,]
        results[i] <- ks.test(x[idx], x[-idx])$stat
        }
        return(results)
}


### approximate permutation
perm.approx.ks <- function(x, n1, R=1000)
{
 #n1: the number of observations in group 1
    #x: combined data
    #output: the sums of all possible permutations
    results <- rep(NA, R)
    N = length(x)
    for(i in 1:R)
    {
        idx = sample(1:N, n1)
        results[i] <- ks.test(x[idx], x[-idx])$stat
        }
    return(results)
}
