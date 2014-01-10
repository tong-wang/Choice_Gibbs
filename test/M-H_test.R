require("LearnBayes")

#log-posterior of theta
logpost.theta <- function(theta, data) {
    if (theta <= 0) return(NaN)
    else {
        d <- data$D
        p <- data$P
        
        loglikelihood <- ifelse(d==1, pexp(q=p, rate=theta, lower.tail=FALSE, log.p=TRUE), pexp(q=p, rate=theta, log.p=TRUE))
        logprior <- dgamma(x=theta, shape=theta.alpha, rate=theta.beta, log=TRUE)
        
        return(sum(loglikelihood) + logprior)
    }
}

  

#generate data
theta.true <- 2
N <- 1000
data <- data.frame(matrix(nrow=N, ncol=0))

data$P <- runif(n=N, min=0, max=2)
for (i in 1:N) {
  data$D[i] <- rbinom(1, 1, 1-pexp(data$P[i], rate=theta.true) )
}


#prior
theta.alpha <- 1
theta.beta <- 1

#sample
varcov <- 0.01
proposal <- list(var=varcov, scale=2)

start <- 1
m <- 10000
s=rwmetrop(logpost.theta,proposal,start,m,data)

s$accept
plot(s$par, type="l")

s.truncated <- s$par[1000:10000]
mean(s.truncated)
hist(s.truncated, 50, freq=FALSE, main="histogram")



require("MCMCpack")
post.samp <- MCMCmetrop1R(logpost.theta, theta.init=1,
                          data=data,
                          thin=1, mcmc=10000, burnin=1000,
                          tune=2,
                          verbose=500, logfun=TRUE)
raftery.diag(post.samp)
plot(post.samp)
summary(post.samp)
