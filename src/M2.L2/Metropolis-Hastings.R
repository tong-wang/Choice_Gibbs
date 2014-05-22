####
# Implementation of the Metropolis-Hastings Algorithm
#    --- Continuous case:
#            (1) MH.mvnorm(): Multivariate Normal proposal distribution
#    --- Discrete case:
#            #(1) discreteMH.norm(): discretized univariate Normal proposal
#            (2) discreteMH.mvnorm(): discretized multivariate Normal proposal
#            #(3) discreteMH.nbinom(): univariate Negative Binomial proposal (for non-negative discrete random variables)
#            #(4) discreteMH.truncated.norm(): discretized univariate truncated Normal proposal (for non-negative discrete random variables)
####


# implementation of Metropolis-Hastings algorithm with multivariate Normal proposal
# tested for one- and high-dimensional distribution
MH.mvnorm <- function (logpost, start, sigma=diag(length(start)), scale=1, nrun=1000, thin=1, ...) 
{
    dim <- length(start)
    MC <- array(0, c(nrun %/% thin, dim))
    b1 <- start
    ll.b1 <- logpost(start, ...)
    a <- chol(sigma)


    accept <- 0
    for (i in 1:nrun) {
        b2 <- as.vector(b1 + scale * t(a) %*% rnorm(dim))
        
        ll.b2 <- logpost(b2, ...)
        ll.ratio <- exp(ll.b2 - ll.b1)

        if (!is.na(ll.ratio)) {
            if (runif(1) <= ll.ratio) {
                ll.b1 <- ll.b2
                b1 <- b2
                accept <- accept + 1
            }
        }
        
        if (i %% thin == 0) {
            MC[i %/% thin, ] <- b1
        }
    }
    accept <- accept/nrun
    
    list(MC = MC, accept = accept)
}




# implementation of discrete Metropolis-Hastings algorithm with discretized symmetric multivariate Normal proposal
# tested for one- or high-dimensional discrete distribution
discreteMH.mvnorm <-function (logpost, start, sigma=diag(length(start)), scale=1, nrun=1000, thin=1, ...) 
{
    dim <- length(start)
    MC <- array(0, c(nrun %/% thin, dim))
    b1 <- start
    ll.b1 <- logpost(start, ...)
    a <- chol(sigma)

    
    accept <- 0
    for (i in 1:nrun) {
        #proposed next point is discrete Multivariate Normal using the current position as its mean, variance is controled by $sigma$ and $scale$
        b2 <- as.vector(round(b1 + scale * t(a) %*% rnorm(dim)))
        
        ll.b2 <- logpost(b2, ...) 
        ll.ratio <- exp(ll.b2 - ll.b1)
        
        if (!is.na(ll.ratio)) {
            if (runif(1) <= ll.ratio) {
                ll.b1 <- ll.b2
                b1 <- b2
                accept <- accept + 1
            }
        }
        
        if (i %% thin == 0) {
            MC[i %/% thin, ] <- b1
        }
    }
    accept <- accept/nrun
    
    list(MC = MC, accept = accept)
}
