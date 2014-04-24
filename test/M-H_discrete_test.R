require("mvtnorm")

# implementation of discrete Metropolis-Hastings algorithm with discretized symmetric univariate Normal proposal
# tested for one- and high-dimensional discrete distribution
discreteMH.norm <-function (logpost, start, scale, nrun, ...) 
{
    dim = length(start)
    MC = array(0, c(nrun, dim))
    b1 = start
    ll.b1 = logpost(start, ...)
    
    
    accept = 0
    for (i in 1:nrun) {
        #proposed next point is discrete Multivariate Normal using the current position as its mean, variance is controled by $scale$
        b2 <- round(rnorm(dim, mean=b1, sd=scale))
        
        ll.b2 = logpost(b2, ...) 
        ll.ratio = exp(ll.b2 - ll.b1)
        
        if (!is.na(ll.ratio)) {
            if (runif(1) <= ll.ratio) {
                ll.b1 = ll.b2
                b1 = b2
                accept = accept + 1
            }
        }
        MC[i, ] = b1
        
    }
    accept = accept/nrun
    
    list(MC = MC, accept = accept)
}

# implementation of discrete Metropolis-Hastings algorithm with discretized symmetric multivariate Normal proposal
# tested for one- or high-dimensional discrete distribution
# not working with COV input
discreteMH.mvnorm <-function (logpost, start, scale, nrun, ...) 
{
    dim = length(start)
    MC = array(0, c(nrun, dim))
    b1 = start
    ll.b1 = logpost(start, ...)

    
    accept = 0
    for (i in 1:nrun) {
        #proposed next point is discrete Multivariate Normal using the current position as its mean, variance is controled by $scale$
        b2 <- round(rmvnorm(1, mean=b1, sigma=diag(dim)*scale))
        
        ll.b2 = logpost(b2, ...) 
        ll.ratio = exp(ll.b2 - ll.b1)
        
        if (!is.na(ll.ratio)) {
            if (runif(1) <= ll.ratio) {
                ll.b1 = ll.b2
                b1 = b2
                accept = accept + 1
            }
        }
        MC[i, ] = b1

    }
    accept = accept/nrun
    
    list(MC = MC, accept = accept)
}



# implementation of discrete Metropolis-Hastings algorithm with (discrete) Negative Binomial proposal
# tested for one- and high-dimensional discrete distribution
discreteMH.nbinom <-function (logpost, start, scale, nrun, ...) 
{
    dim = length(start)
    MC = array(0, c(nrun, dim))
    b1 = start
    ll.b1 = logpost(start, ...)
    
    
    accept = 0
    for (i in 1:nrun) {
        #proposed next point is Negative Binomial using the current position as its mean, variance is controled by $scale$
        b2 <- rnbinom(dim, mu=b1, size=scale)
        
        ll.b2 = logpost(b2, ...) 
        #since proposal is asymmetric, need to adjust the jumping probability accordingly        
        ll.ratio = exp(ll.b2 - ll.b1 + sum(dnbinom(b1, mu=b2, size=scale, log=TRUE) - dnbinom(b2, mu=b1, size=scale, log=TRUE)) )
        
        if (!is.na(ll.ratio)) {
            if (runif(1) <= ll.ratio) {
                ll.b1 = ll.b2
                b1 = b2
                accept = accept + 1
            }
        }
        MC[i, ] = b1

    }
    accept = accept/nrun
    
    list(MC = MC, accept = accept)
}


# implementation of discrete Metropolis-Hastings algorithm with discretized zero-truncated Normal proposal
# tested for one-dimensional discrete distribution
discreteMH.truncated.norm <-function (logpost, start, scale, nrun, ...) 
{
    require("truncnorm")
    
    dim = length(start)
    MC = array(0, c(nrun, dim))
    b1 = start
    ll.b1 = logpost(start, ...)
    
    
    accept = 0
    for (i in 1:nrun) {
        # proposed next point is discrete Multivariate Normal using the current position as its mean, variance is controled by $scale$
        # any negative value proposed will be discarded
        b2 <- round(rtruncnorm(dim, mean=b1, sd=scale, a=-0.5))
        
        ll.b2 = logpost(b2, ...) 
        ll.ratio = exp(ll.b2 - ll.b1 + sum(log(ptruncnorm(b1+0.5, mean=b2, sd=scale, a=-0.5) - ptruncnorm(b1-0.5, mean=b2, sd=scale, a=-0.5)) - log(ptruncnorm(b2+0.5, mean=b1, sd=scale, a=-0.5) - ptruncnorm(b2-0.5, mean=b1, sd=scale, a=-0.5)) ) )
        
        if (!is.na(ll.ratio)) {
            if (runif(1) <= ll.ratio) {
                ll.b1 = ll.b2
                b1 = b2
                accept = accept + 1
            }
        }
        MC[i, ] = b1
        
    }
    accept = accept/nrun
    
    list(MC = MC, accept = accept)
}




#sample
nrun <- 1000000

## test univariate discrete distributions
s <- discreteMH.truncated.norm(logpost=dpois, start=1, scale=10, nrun=nrun, lambda=20, log=TRUE)
#s <- discreteMH.truncated.norm(logpost=dbinom, start=1, scale=10, nrun=nrun, size=100, prob=0.2, log=TRUE)
#s <- discreteMH.truncated.norm(logpost=dnbinom, start=1, scale=10, nrun=nrun, size=100, mu=20, log=TRUE)

#s <- discreteMH.nbinom(logpost=dpois, start=1, scale=5, nrun=nrun, lambda=20, log=TRUE)
#s <- discreteMH.nbinom(logpost=dbinom, start=1, scale=5, nrun=nrun, size=100, prob=0.2, log=TRUE)
#s <- discreteMH.nbinom(logpost=dnbinom, start=1, scale=5, nrun=nrun, size=100, mu=20, log=TRUE)


s$accept
plot(s$MC, type="l")

s.truncated <- s$MC[(0.2*nrun):nrun]
mean(s.truncated)
hist(s.truncated, breaks=100, freq=FALSE, main="histogram")

density.true <- dpois(0:49, 20)
#density.true <- dbinom(0:49, 100, 0.2)
#density.true <- dnbinom(0:49, 100, mu=20)
density.MHsim <- hist(s.truncated, breaks=0:50, plot=F, right=F)$density
plot(0:49, (density.MHsim - density.true)/density.true )


sim <- rpois(length(s.truncated), 20)
#sim <- rbinom(length(s.truncated), 100, 0.2)
#sim <- rnbinom(length(s.truncated), 100, mu=20)
density.rsim <- hist(sim, breaks=0:50, plot=F, right=F)$density
plot(0:49, (density.rsim - density.true)/density.true  )


plot(0:49, (density.MHsim - density.rsim)/density.true )





## test multivariate discrete distributions
ll.multinom <- function(x, size, prob, log)
{
    allx <- c(x, size-sum(x))
    if (any(allx < 0)) 
        return(if (log) -Inf else 0)
    else
        return(dmultinom(allx, prob=prob, log=log))
}

#s2 <- discreteMH.norm(logpost=ll.multinom, start=c(1,1), scale=c(10,5), nrun=nrun, size=100, prob=c(0.5, 0.1, 0.4), log=TRUE)
s2 <- discreteMH.truncated.norm(logpost=ll.multinom, start=c(1,1), scale=c(10,5), nrun=nrun, size=100, prob=c(0.5, 0.1, 0.4), log=TRUE)
#s2 <- discreteMH.nbinom(logpost=ll.multinom, start=c(1,1), scale=c(100,25), nrun=nrun, size=100, prob=c(0.5, 0.1, 0.4), log=TRUE)

s2$accept
plot(s2$MC[,1], type="l")
plot(s2$MC[,2], type="l")

s2.truncated <- s2$MC[(0.2*nrun):nrun, ]
colMeans(s2.truncated)
hist(s2.truncated[,1], breaks=100, freq=FALSE, main="histogram")
hist(s2.truncated[,2], breaks=100, freq=FALSE, main="histogram")


density1.true <- dbinom(0:100, size=100, prob=0.5)
density1.MHsim <- hist(s2.truncated[,1], breaks=0:101, plot=F, right=F)$density
plot(0:100, (density1.MHsim - density1.true)/density1.true )

density2.true <- dbinom(0:40, size=100, prob=0.1)
density2.MHsim <- hist(s2.truncated[,2], breaks=0:41, plot=F, right=F)$density
plot(0:40, (density2.MHsim - density2.true)/density2.true )
