####
# Scenario 1. No-purchase is not observed (Binary Choice case)
# !!! This case does not seem to converge !!!
####

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src")

require("MCMCpack")
require("mvtnorm")



### Known parameters
M <- 2 # number of alternatives (alternative 2 is dummy for no-purchase)
L <- 1 # number of covariates
K <- 360 # number of periods

#X is the attributes of the alternatives; in each period, [Xij] is an (M-1)*L matrix, i=1...M-1, j=1...L.
#by row: [X11 X12]
X_Mean <- 5

#for each of the K periods, generate an X matrix
#dim of X_Mat is K, (M-1)*L
X_Mat <- rnorm(K, mean=X_Mean, sd=1)


## true values of the parameters to be estimated
# beta is the MNL coefficient
beta <- 0.1; # L-dimensional
# lambda is the poisson demand rate in each perid
lambda <- 100


### simulate data
# simulate number of demand per period, ~Poisson(lambda)
N <- rpois(K, lambda) 

#in each period, score of choice 1 and 2 (col) is exp(X*beta)
#dim of score is M by K
score <- rbind(exp(X_Mat * beta), rep(1, K))

#choice probabilities in each period (M by K)
choice.prob <- apply(score, 2, function(x) x/sum(x)) # M*N matrix
#simulate actual choice in each period (M by K)
choice.mat <- matrix(0, nrow=M, ncol=K)
for (k in 1:K) {
    choice.mat[, k] <- rmultinom(1, N[k], choice.prob[, k])    
}    

rowMeans(choice.mat)

observation1 <- choice.mat[1:M-1,]



#log-posterior of beta, to be called by M-H algorithm
logpost.beta <- function(beta, data) {
    
    if (any(beta<0))
        return(-Inf)
    else {
        score <- rbind(exp(X_Mat * beta), rep(1, K))
        choice.prob <- apply(score, 2, function(x) x/sum(x)) # M*N matrix
        
        logLikelihood <- data*log(choice.prob)

        logprior <- dnorm(log(beta), mean=beta.mu, sd=beta.sd, log=TRUE)
        
        return(sum(logLikelihood) + logprior)
    }
}


logpost.d0 <- function(d0, data, lambda, beta, k) {
    
    if (any(d0<0))
        return(-Inf)
    else {
        score <- c(exp(X_Mat[k] * beta), 1)
        choice.prob <- score/sum(score) # M*N matrix
        
        dd <- c(data,d0)
        nn <- sum(dd)
    
        return(dmultinom(x=dd, prob=choice.prob, log=TRUE) + dpois(nn, lambda, log=TRUE))
    }
}




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




sample = function(data, parameters, nrun=1000) {

    
    # initialize arrays to save samples
    lambdas = array(0, dim=c(nrun, 1))
    d0s = array(0, dim=c(nrun, K))
    betas = array(0, dim=c(nrun, L))
    
    # initial parameters
    lambda1 <- parameters$lambda
    d01 <- parameters$d0
    beta1 <- parameters$beta
    
    
    # start sampling loop
    for(i in 1:nrun) {

        #update posterior of lambda by conjugacy
        alpha2 <- lambda.alpha + sum(data) + sum(d01)
        beta2 <- lambda.beta + K
        
        #simulate lambda2
        lambda2 <- rgamma(1, shape=alpha2, rate=beta2)
        
        
        
        #simulate beta2 by Metropolis-Hastings
        beta2 <- MCMCmetrop1R(logpost.beta, theta.init=beta1,
                         data=rbind(data,d01),
                         thin=1, mcmc=1, burnin=10, tune=0.01,
                         verbose=0,  V=matrix(c(1),1,1))[1,]


        #simulate d0 by discrete Metropolis-Hastings
        d02 <- d01
        d0.accept <- rep(0, K)
        for (j in 1:K) {
            sim <- discreteMH.norm(logpost.d0, start=d01[j], scale=15, nrun=10, 
                            data=data[j], lambda=lambda2, beta=beta2, k=j)
            d02[j] <- sim$MC[10]
            d0.accept[j] <- sim$accept
        }
        cat("discreteMH acceptance rate was ", mean(d0.accept), "\n\n")
        
        
        

        # save the samples obtained in the current iteration
        lambdas[i,] = lambda2
        d0s[i,] = d02
        betas[i,] = beta2
        
        cat("Run:", i, "\tlambda:", lambda2, "\tbeta2:", beta2, "\td0[1]:", d02[1], "\n", sep=" ")
        
        lambda1 <- lambda2
        d01 <- d02
        beta1 <- beta2
    }
    
    # results
    return(list(lambdas=lambdas, d0s=d0s, betas=betas))

} # end function 'sample'




### initialize input before sampling
# beta prior ~ logN(beta.mu, beta.sd)
beta.mu <- -2
beta.sd <- 10


# lambda ~ Gamma(lambda.alpha, lambda.beta)
lambda.alpha <- 0.005
lambda.beta <- 0.0001


## initial sampling input
param0 <- list(beta=rep(0, L), lambda=lambda.alpha/lambda.beta, d0=rep(10,K))
nrun <- 5000
burnin <- 0.5



### sample
z1 <- sample(data=observation1, parameters=param0, nrun=nrun)


save(z1, observation1, file="MNL_Scenario1.binary.RData")




### Visualize results
start <- burnin*nrun+1

#plot lambda
samples.lambda <- z1$lambdas
plot(samples.lambda, type="l")


samples.lambda.truncated <- samples.lambda[start:nrun,]
quantile(samples.lambda.truncated, c(.025,.5,.975))
mean(samples.lambda.truncated)
hist(samples.lambda.truncated)



#plot beta
samples.beta <- z1$betas
plot(samples.beta, type="l")


samples.beta.truncated <- samples.beta[start:nrun,]
quantile(samples.beta.truncated, c(.025,.5,.975))
mean(samples.beta.truncated)
hist(samples.beta.truncated)



#plot d0[1]
samples.d0 <- data.frame(z1$d0s)
plot(samples.d0$X1, type="l")
plot(samples.d0$X2, type="l")
plot(samples.d0$X3, type="l")
plot(samples.d0$X4, type="l")
plot(samples.d0$X5, type="l")


samples.d0.truncated <- samples.d0[start:nrun,]
quantile(samples.d0.truncated$X1, c(.025,.5,.975))
quantile(samples.d0.truncated$X2, c(.025,.5,.975))
quantile(samples.d0.truncated$X3, c(.025,.5,.975))
quantile(samples.d0.truncated$X4, c(.025,.5,.975))
quantile(samples.d0.truncated$X5, c(.025,.5,.975))
colMeans(samples.d0.truncated)
hist(samples.d0.truncated$X1)
hist(samples.d0.truncated$X2)
hist(samples.d0.truncated$X3)
hist(samples.d0.truncated$X4)
hist(samples.d0.truncated$X5)



