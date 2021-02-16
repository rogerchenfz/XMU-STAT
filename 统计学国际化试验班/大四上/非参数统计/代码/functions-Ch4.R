

            #########################################
            ##  Paired-Comparison Permutation Test ##
            ##          permuting the signs of Di   ##
           #########################################
getbinvecs <- function(n)
{
    ### obtain the 2^n assignments of 0/1
    results <- matrix(NA,2^n,n)
    for (i in 1:(2^n)) 
    {   
        remain <- i-1
        for (j in 1:n)
       {
          results[i,j] <- remain %/% (2^(n-j)) #%/%: integer division
          remain <- remain %% (2^(n-j)) #%%: remainder
       }
   }
   results  
}

perm.dbar <- function(d)
{
    # exact permutation distribution of dbar
    # return the Dbar* based on 2^n permutations
    n <- length(d)
    junk <- getbinvecs(n)
    signs <- junk - !junk  ### permuted signs (+/-)
    apply( t(abs(d) * t(signs)) , 1, mean) 
}

perm.approx.dbar <- function(d,R)
{
    # an approximation of perm.dbar for large n
    # R: number of permutations
    n <- length(d)
    junk <- matrix(NA,R,n)
    for (i in 1:R) junk[i,] <- rbinom(n,1,.5)
    signs <- junk - !junk
    apply( t(abs(d) * t(signs)) , 1, mean) 
}



###calculate the signed rank
signrank <- function(x){rank(abs(x))* sign(x) }




            #################################
            ##  Permutation Test for RCBD ##
            #################################
getmeans <- function(x,grps)
{
    ## estimate the group-specific means
    junk <- table(grps)
    meanvec <- rep(NA,length(junk))
    for (i in 1:length(junk)) meanvec[i] <- mean(x[grps==names(junk)[i]])
    meanvec
}

### calculate observed SSTM
getSSTM <- function(x,grps)
{
    ## estimate SSTM=sum of squared group-means 
    sum( (getmeans(x,grps))^2 )
    }

### obtain the permutation distribution of SSTM under H0
###             i.e. obtain SSTM*'s with R permutations
perm.approx.RCBD <- function(x,grps,blocks,R)
{
    # obtain the aprpoximate permutation distribution of SSTM under the null hypothesis
    mat <- cbind(x,grps,blocks)
    mat <- mat[order(mat[,3]),] 
    k <- length(table(grps))
    b <- length(table(blocks))
    results <- rep(NA,R)
    for (i in 1:R)
    {
        junk <- rep(NA,b*k)
        for (j in 1:b) junk[k*(j-1)+(1:k)] <- sample(1:k,k)
        results[i] <- getSSTM(mat[,1],junk)
    }
    results
}

            #################################
            ##   Friedman Test for RCBD    ##
            #################################


rankinblock <- function(x, blocks)  ## assumes RCBD
{
    ### obtain the ranks of observations within each block
   mat <- cbind(x,blocks,1:length(x))
   mat <- mat[order(mat[,2]),]
    b <- length(table(blocks))
    n <- length(x)
   for (j in 1:b) mat[n/b*(j-1)+(1:(n/b)),1] <- rank( mat[n/b*(j-1)+(1:(n/b)),1] )
   mat <- mat[order(mat[,3]),]
   mat[,1]
}



            ############################################
            ##  Page's test for ordered alternatives  ##
            ############################################

getsums <- function(x,grps)
{
    junk <- table(grps)
    sumvec <- rep(NA,length(junk))
    for (i in 1:length(junk)) sumvec[i] <- sum(x[grps==names(junk)[i]])
    sumvec
}

getPage <- function(x,grps,blocks)
{
    sum( (1:(length(table(grps))))*getsums(rankinblock(x,blocks),grps) )
    }

perm.approx.Page <- function(x,grps,blocks,R)
{
    mat <- cbind(x,grps,blocks)
    mat <- mat[order(mat[,3]),]  # sort matrix rows by block
    k <- length(table(grps))
    b <- length(table(blocks))
    results <- rep(NA,R)
    for (i in 1:R)
    {
        junk <- rep(NA,b*k)
        for (j in 1:b) junk[k*(j-1)+(1:k)] <- sample(1:k,k)
        results[i] <- getPage(mat[,1],junk,blocks)
    }
    results
}




