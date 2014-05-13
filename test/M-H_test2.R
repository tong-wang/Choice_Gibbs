require("LearnBayes")
require("MCMCpack")

MH.mvnorm <- function (logpost, start, sigma, scale, nrun, ...) 
{
    dim = length(start)
    MC = array(0, c(nrun, dim))
    a = chol(sigma)
    b1 = start
    ll.b1 = logpost(start, ...)
    
    
    accept = 0
    for (i in 1:nrun) {
        b2 = as.vector(b1 + scale * a %*% rnorm(dim))
        
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




#sample
nrun <- 100000

## test univariate discrete distributions
s1 <- MH.mvnorm(logpost=dgamma, sigma=1, scale=10, start=1, nrun=nrun, shape=10, rate=1, log=TRUE)
s1$accept
plot(s1$MC, type="l")
s1.truncated <- s1$MC[20000:100000]
mean(s1.truncated)

s2 <- rwmetrop(logpost=dgamma, proposal=list(var=1,scale=10), start=1, m=nrun, shape=10, rate=1, log=TRUE)
s2$accept
plot(s2$par, type="l")
s2.truncated <- s2$par[20000:100000]
mean(s2.truncated)


s3 <- MCMCmetrop1R(dgamma, V=diag(1), tune=3, theta.init=1, mcmc=nrun, verbose=500, shape=10, rate=1, log=TRUE)
plot(s3)

s3 <- MCMCmetrop1R(fun=dnorm, V=diag(1), tune=0.0001, theta.init=1, mcmc=nrun, verbose=500, logfun=TRUE, log=TRUE)





density.true <- dgamma(0:49, shape=10, rate=1)
#density.true <- dbinom(0:49, 100, 0.2)
#density.true <- dnbinom(0:49, 100, mu=20)
density.1 <- hist(s1.truncated, breaks=0:50, plot=F, right=F)$density
plot(0:49, (density.1 - density.true)/density.true )
density.2 <- hist(s2.truncated, breaks=0:50, plot=F, right=F)$density
plot(0:49, (density.2 - density.true)/density.true )


