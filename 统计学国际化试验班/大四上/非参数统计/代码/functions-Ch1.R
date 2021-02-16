### Chapter 1

### function for constructing the finite sample (1-alpha) confidence interval for the median theta

# A distribution-free confidence interval on the median assuming only iid is based on the order statistics and the binomial distribution. See Lehmann
# (Nonparametrics: Statistical Methods Based on Ranks, Holden-Day, 1975, p. 182-183). 

# The following R function can be used to set up such a 100(1-alpha)% distribution-free confidence interval on the median:

conf.med<-function(x, alpha)
{
  v <- sort(x, na.last = NA)
  n <- length(x)
  if(n > 0) {
    m <- median(x)
    l <- qbinom(alpha/2, n, 0.5)
    if(l > 0) 
    {
      u = n-l+1
      r <- c(m, v[l], v[u])
    }
    else r <- c(m, NA, NA)
  }
  else r <- c(NA, NA, NA)
  r <- as.data.frame(list(median = r[1], lower = r[2], upper = r[3]))
  class(r) <- "table"
  r
}

# The following R function can be used to set up such a 100(1-alpha)% distribution-free confidence interval on the p-th quantile:

conf.quantile<-function(x, alpha, p)
{
  # 1-alpha confidence interval for theta_p: the pth quantile, e.g. p=0.75 corresponding to the 75th percentile 
  v <- sort(x, na.last = NA)
  n <- length(x)
  if(n > 0) {
    m <- quantile(x, p)
    l <- qbinom(alpha/2, n, p)
    u = qbinom(1-alpha/2, n, p)+1
    if(l > 0) 
    {
      r <- c(m, v[l], v[u])
    }
    else r <- c(m, NA, NA)
  }
  else r <- c(NA, NA, NA)
  r <- as.data.frame(list(p_quantile = r[1], lower = r[2], upper = r[3]))
  class(r) <- "table"
  r
}


## An example
# x=c(20,  22,  31,  34,  37,  43,  45,  48,  51,  54,  62,  83,  84,  87,  93, 117, 119, 150, 200, 216)
# p=0.5, conf.med and conf.quantile are identical
# alpha=0.01; p=0.5
# conf.med(x, alpha)
# conf.quantile(x, alpha, p)

# p=0.75
# alpha=0.01; p=0.75
# conf.quantile(x, alpha, p)

### another example (1.2.4 on the slides)
# x=c(79,74,88,80,80,66,65,86,84,80,78,72,71,74,86,96,77,81,76,80,76,75,78,87,87,74,85,84,76,77,76,74,85,74,76,77,76,74,81,76)
# conf.quantile(x, alpha=0.05, p=0.75)
