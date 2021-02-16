

            ####################################################
            ##  F-test and Permutation F-test ##
            ####################################################
#Functions Needed
getmeans <- function(x,grps)
{
    #estimate the group-specific means
    junk <- table(grps)
    meanvec <- rep(NA,length(junk))
    for (i in 1:length(junk)) meanvec[i] <- mean(x[grps==names(junk)[i]])
    meanvec
}

getF <- function(x,grps)
{
    # estimate the F test statistic for one-way layout
    junk <- table(grps)
    SST <- sum(junk * (getmeans(x,grps))^2) -length(x)*(mean(x))^2
    SST / (length(junk)-1) / ( (sum((x-mean(x))^2) - SST)/(length(x)-length(junk)))
}

###### approximate permutation distribution of F statistic for one-way layout
perm.approx.F <- function(x,grps,R)
{
    #x: observed response vector
    #grps: group number
    #R: number of permutations
    results <- rep(NA,R)
    #for (i in 1:R) results[i] <- summary(aov(x[sample(1:(length(x)),length(x))] ~ factor(grps)))[[1]][1,4]
    for (i in 1:R) results[i] <- getF(x[sample(1:(length(x)),length(x))],grps)
      results
}


            #################################################
            #  Multiple Comparison Permutation Tests  #
            #################################################

perm.approx.diff <- function(x,n1,n2,R)
{
### approximate permutation of pairwise group means
### n1: sample size of group 1
### n2: sample size of group 2
### x: combined data
### R: number of random permutations
### assumption: the observations in x are exchangeable under H0
    results <- rep(NA,R)
    for (i in 1:R)
    {
        temp <- x[sample(1:(length(x)),n1+n2)]
        results[i] <- mean(temp[1:n1]) - mean(temp[(n1+1):(n1+n2)])
    }
    results
}

Bonf.adj = function(x, grps, k, alpha, R=1000, test=c("perm.meandiff", "ttest"))
{
    #Bonferroni correction
    #alpha: overall significance level
    #R: number of permutation
    #test="perm.meandiff": pvalues are calculated through permutation test based on mean difference from two groups
    #test="ttest": pvalues are calculated through pairwise two-sample t-test (requires normal distirbution assumption)

    #output:
    #adj.alpha: adjsuted significance level for each pairwise comparison
    # bonf.pvals: p-val from pairwise comparison
    # sig: TRUE means the two pairs are significantly different

    nn <- table(factor(grps))
    ## obtain the permutation p-value for comparing each pair (i, j)
    bonf.pvals = matrix(NA, k, k)
    if(test=="perm.meandiff")
    {
    for (i in 2:k) {
        for (j in 1:(i-1)) {
            diff.perms = perm.approx.diff(x[grps %in% c(i, j)], nn[i], nn[j], R=R)
            diff.obs = mean(x[grps==i]) - mean(x[grps==j])
            bonf.pvals[i,j] <- mean(abs(diff.perms) >= abs(diff.obs))
          }}
    }

    if(test=="ttest")
    {
        for (i in 2:k) {
        for (j in 1:(i-1)) {
            bonf.pvals[i,j] = t.test(x[grps==i], x[grps==j])$p.value
          }}
    }
    # For Bonferroni, compare each of the p-values against the adjusted level: alpha/(k(k-1)/2)
    (adj.alpha = alpha / (k*(k-1)/2))
    sig = (bonf.pvals <= adj.alpha)
    out = list(adj.alpha=adj.alpha, bonf.pvals=bonf.pvals, sig=sig)
    return(out)
}

Fisher.LSD = function(x, grps, k, alpha=0.05, R=1000)
{
    #Fisher's Protected LSD
    ##
    ### overall test (permutation F-test)
    Fobs = getF(x, grps)
    perm.F = perm.approx.F(x, grps, R=R)
    perm.pval = mean(perm.F >= Fobs)
   nn <- table(factor(grps))

    if(perm.pval<=alpha)
    {
        cat("the overall test is significant, proceed with pairwise comparison", "\n")

        unadj.pvals <- matrix(NA,k,k)
        for (i in 2:k) {
            for (j in 1:(i-1)){
                diff.perms = perm.approx.diff(x[grps %in% c(i, j)], nn[i], nn[j], R=R)
                diff.obs = mean(x[grps==i]) - mean(x[grps==j])
                unadj.pvals[i,j] <- mean(abs(diff.perms) >= abs(diff.obs))
          }}
        sig = (unadj.pvals <= alpha)
        out = list(pval.overall= perm.pval, pval.pairwise= unadj.pvals, sig=sig)
    }
    else if(perm.pval>alpha)
    {
        cat("the overall test is not significant, do not proceed with pairwise comparison", "\n")
        out = list(pval.overall= perm.pval)
        }
    return(out)
}



getmaxTij <- function(x, grps, MSE)
{
    # estimate the maximum of Tij (pairwise mean diff) of a given data x
    trtmeans <- getmeans(x,grps)
   nn <- table(factor(grps))
    k <- length(trtmeans)
    Tijs <- matrix(NA,k,k)
    for (i in 2:k) {
    for (j in 1:(i-1)){
         Tijs[i,j] <- abs(trtmeans[i] - trtmeans[j])/sqrt(MSE/2 * (1/nn[i] + 1/nn[j]))
     }}
    max(Tijs,na.rm=T)
}

perm.approx.maxTij <- function(x,grps,MSE,R)
{
    ### obtain the null permutation distribution of maxTij
    results <- rep(NA,R)
    for (i in 1:R) results[i] <- getmaxTij(x[sample(1:(length(x)),length(x))],grps,MSE)
    results
}

Tukey.HSD = function(x, grps, k, alpha=0.05, R=1000)
{
    #Tukey's HSD
    #summary(aov(x ~ factor(grps)))
    nn <- table(factor(grps))
    trtmeans <- getmeans(x,grps)

    (MSE <- summary(aov(x ~ factor(grps)))[[1]][2,3])
    ### observed Tij
    Tijs <- matrix(NA,k,k)
    for (i in 2:k){
        for (j in 1:(i-1)){
             Tijs[i,j] <- abs(trtmeans[i] - trtmeans[j])/sqrt(MSE/2 * (1/nn[i] + 1/nn[j]))
             }}

    ### observed maxTij
    #getmaxTij(x,grps,MSE)

    ### permutation maxTij
    perm.maxTij <- perm.approx.maxTij(x,grps,MSE,R)

    pvalsTij <- matrix(NA,k,k)
    for (i in 2:k){
    for (j in 1:(i-1)){
        pvalsTij[i,j] <- mean(perm.maxTij >= Tijs[i,j])
        }}

    ### compare the pairwise pvalue with alpha
    sig = (pvalsTij <= alpha)

    out = list(sig=sig, pvalsTij= pvalsTij)
    return(out)
}



Bonf.adj.rank = function(x, grps, k, alpha=0.05, R=1000)
{
    #Bonferroni correction based on ranks
    #alpha: overall significance level
    #R: number of permutation
    #pvalues are calculated through permutation test based on mean difference from two groups

    #output:
    #adj.alpha: adjsuted significance level for each pairwise comparison
 # bonf.pvals: p-val from pairwise comparison
    # sig: TRUE means the two pairs are significantly different

    ## obtain the permutation p-value for comparing each pair (i, j)
    nn <- table(factor(grps))  #size of each treatment group
    bonf.pvals = matrix(NA, k, k)
    for (i in 2:k) {
        for (j in 1:(i-1)) {
            diff.perms = perm.approx.diff(rank(x[grps %in% c(i, j)]), nn[i], nn[j], R=R)
            x1 = x[grps==i]
            x2 = x[grps==j]
            tmp.rank = rank(c(x1, x2))
            diff.obs = mean(tmp.rank[1:nn[i]]) - mean(tmp.rank[-(1:nn[i])])
            bonf.pvals[i,j] <- mean(abs(diff.perms) >= abs(diff.obs))
          }}
    # For Bonferroni, compare each of the p-values against the adjusted level: alpha/(k(k-1)/2)
    (adj.alpha = alpha / (k*(k-1)/2))
    sig = (bonf.pvals <= adj.alpha)

    out = list(adj.alpha=adj.alpha, bonf.pvals=bonf.pvals, sig=sig)
    return(out)
}


Fisher.LSD.rank = function(x, grps, k, alpha=0.05, R=1000)
{
    #Tukey's HSD based on ranks
    ### overall test (kruskal-Wallis test)
    pval.overall = kruskal.test(x, grps)$p.val
    nn <- table(factor(grps))

    if(pval.overall<=alpha)
    {
        cat("the overall test is significant, proceed with pairwise comparison", "\n")

    unadj.pvals <- matrix(NA,k,k)
    for (i in 2:k) {
        for (j in 1:(i-1)){
            diff.perms = perm.approx.diff(rank(x[grps %in% c(i, j)]), nn[i], nn[j], R=R)
            x1 = x[grps==i]
            x2 = x[grps==j]
            tmp.rank = rank(c(x1, x2))
            diff.obs = mean(tmp.rank[1:nn[i]]) - mean(tmp.rank[-(1:nn[i])])
            unadj.pvals[i,j] <- mean(abs(diff.perms) >= abs(diff.obs))
          }}
    ### compare pairwise $p$-values with alpha
        sig = (unadj.pvals <= alpha)
        out = list(pval.overall= pval.overall, pval.pairwise= unadj.pvals, sig=sig)
    }
    else if(pval.overall>alpha)
    {
        cat("the overall test is not significant, do not proceed with pairwise comparison", "\n")
        out = list(pval.overall= pval.overall)
        }
    return(out)
}

Tukey.HSD.rank = function(x, grps, k, alpha=0.05, R=1000)
{
    #Tukey's HSD based on ranks
    nn <- table(factor(grps))
    trtmeans <- getmeans(rank(x),grps)

    ##obtain pairwise mean difference of ranks (ranging from 1:(ni+nj))
    meandiffs <- matrix(NA,k,k)
    for (i in 2:k) {
        for (j in 1:(i-1)) {
        idx1 = which(grps==i)
        idx2 = which(grps==j)
        idx = c(idx1, idx2)
        temp <- rank(x[c(idx1, idx2)])
        rank.i = temp[1:length(idx1)]
        rank.j = temp[-(1:length(idx1))]
        meandiffs[i,j] <- abs(mean(rank.i) - mean(rank.j))
    }}


    MSE <- summary(aov(rank(x) ~ factor(grps)))[[1]][2,3]

    ### observed Tij
    Tijs <- matrix(NA,k,k)
    for (i in 2:k){
        for (j in 1:(i-1)){
             Tijs[i,j] <- abs(trtmeans[i] - trtmeans[j])/sqrt(MSE * (1/nn[i] + 1/nn[j]))
            }}


    ### observed maxTij
    #getmaxTij(rank(x), grps, MSE)

    ### permutation maxTij
    perm.maxTij <- perm.approx.maxTij(rank(x), grps, MSE, R)

    pvalsTij <- matrix(NA,k,k)
    for (i in 2:k){
        for (j in 1:(i-1)){
            pvalsTij[i,j] <- mean(perm.maxTij >= Tijs[i,j])
        }}

    ### compare the pairwise pvalue with alpha
    sig = (pvalsTij <= alpha)

    out = list(sig=sig, pvalsTij= pvalsTij, meandiffs= meandiffs)
    return(out)
}



            ################################
            #   Ordered Alternatives  #
            ################################
##
##### Jonckheere-Terpstra test for ordered alternative
##
JT.ranksum = function(x, grps, k, R=1000)
{
    ##Jonckheere-Terpstra test for Ha: mu1 <= mu2 <= ... <= muk
    ### Tij: Wilcoxn rank-sum test statistic
    Tijs = matrix(NA, k, k)
    for(i in 1:(k-1)){
        xi = x[grps==i]
        for(j in (i+1):k){
            xj = x[grps==j]
            Tijs[i,j] = sum(rank(c(xi, xj))[-(1:length(xi))])
            }
        }
    JTobs <- sum(Tijs,na.rm=T)

    getJT.ranksum <- function(x, grps, k)
    {
        Tijs = matrix(NA, k, k)
        for(i in 1:(k-1)){
        xi = x[grps==i]
        for(j in (i+1):k){
            xj = x[grps==j]
            Tijs[i,j] = sum(rank(c(xi, xj))[-(1:length(xi))])
            }
        }
        sum(Tijs,na.rm=T)
    }

    perm.approx.JT.ranksum <- function(x, grps, k, R)
    {
        results <- rep(NA,R)
        for (i in 1:R)
            results[i] <- getJT.ranksum(x[sample(1:(length(x)),length(x))], grps, k)
        results
    }

    permJT <- perm.approx.JT.ranksum(x, grps,k, R)

    pvalJT <- mean(permJT >= JTobs)

    out = list(pvalJT=pvalJT, JTobs= JTobs)
    return(out)
}


JT.MW = function(x, grps, k, R=1000)
{
    ##Jonckheere-Terpstra test for Ha: mu1 <= mu2 <= ... <= muk
    ### Tij: Mann-Whiteney statistic
    Tijs = matrix(NA, k, k)
    for(i in 1:(k-1)){
        xi = x[grps==i]
        for(j in (i+1):k){
            xj = x[grps==j]
            nj = length(xj)
            Tijs[i,j] = sum(rank(c(xi, xj))[-(1:length(xi))]) - nj*(nj+1)/2
            }
        }
    JTobs <- sum(Tijs,na.rm=T)

    getJT.MW <- function(x, grps, k)
    {
        Tijs = matrix(NA, k, k)
        for(i in 1:(k-1)){
        xi = x[grps==i]
        for(j in (i+1):k){
            xj = x[grps==j]
            nj = length(xj)
            Tijs[i,j] = sum(rank(c(xi, xj))[-(1:length(xi))])- nj*(nj+1)/2
            }
        }
        sum(Tijs,na.rm=T)
    }

    perm.approx.JT.MW <- function(x, grps, k, R)
    {
        results <- rep(NA,R)
        for (i in 1:R)
            results[i] <- getJT.MW(x[sample(1:(length(x)),length(x))], grps, k)
        results
    }

    permJT <- perm.approx.JT.MW(x, grps,k, R)
    pvalJT <- mean(permJT >= JTobs)

    out = list(pvalJT=pvalJT, JTobs= JTobs)
    return(out)
}
