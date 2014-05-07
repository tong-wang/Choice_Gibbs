####
# Scenario 1. No-purchase is not observed 
#
####

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src")

require("MCMCpack")
require("mvtnorm")



## Load simulated choice data (NEED TO RUN MNL_InitData.R TO GENERATE THE DATA FIRST)
load(file="MNL_InitData.RData")



observation1 <- choice.mat[1:M-1,]


#log-posterior of beta, to be called by M-H algorithm
logpost.beta <- function(beta, data) {
    
    if (any(beta<0))
        return(-Inf)
    else {
        score <- apply(X_Mat, 1, function (x) exp(matrix(c(x,0,0), nrow=M, ncol=L, byrow=TRUE) %*% beta))
        choice.prob <- apply(score, 2, function(x) x/sum(x)) # M*N matrix
        
        logLikelihood <- data*log(choice.prob)
        
        logprior <- dmvnorm(log(beta), mean=beta.mu, sigma=beta.sg, log=TRUE)
        
        return(sum(logLikelihood) + logprior)
    }
}



logpost.d0 <- function(d0, data, lambda, beta, k) {
    
    if (any(d0<0))
        return(-Inf)
    else {
        score <-  exp(matrix(c(X_Mat[k,],0,0), nrow=M, ncol=L, byrow=TRUE) %*% beta)
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
                            thin=1, mcmc=1, burnin=10, tune=0.005,
                            verbose=0,  V=matrix(c(1,0,0,1),2,2))[1,]
                            #optim.lower=1e-6, optim.method="L-BFGS-B")[1,]


        #simulate d0 by discrete Metropolis-Hastings
        d02 <- d01
        d0.accept <- rep(0, K)
        for (j in 1:K) {
            sim <- discreteMH.norm(logpost.d0, start=d01[j], scale=10, nrun=10, 
                              data=data[,j], lambda=lambda2, beta=beta2, k=j)
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
# beta prior ~ logN(beta.mu, beta.sg)
beta.mu <- c(-2, -2)
beta.sg <- matrix(c(10, 0, 0, 10), 2, 2)

# lambda ~ Gamma(lambda.alpha, lambda.beta)
lambda.alpha <- 0.005
lambda.beta <- 0.0001


## initial sampling input
param0 <- list(beta=rep(0.1, L), lambda=lambda.alpha/lambda.beta, d0=rep(10,K))
nrun <- 5000
burnin <- 0.5



### sample
z1 <- sample(data=observation1, parameters=param0, nrun=nrun)


save(z1, observation1, file="MNL_Scenario1.RData")





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
samples.beta <- data.frame(z1$betas)
plot(samples.beta$X1, type="l")
plot(samples.beta$X2, type="l")

plot(c(1,nrun), c(-0.05, 0.2), type="n", xlab="Samples", ylab="a", xaxt="n", yaxt="n")
lines(samples.beta$X1, col="GREEN")
lines(samples.beta$X2, col="BLUE")
axis(side=1)
axis(side=2)

samples.beta.truncated <- samples.beta[start:nrun,]
quantile(samples.beta.truncated$X1, c(.025,.5,.975))
quantile(samples.beta.truncated$X2, c(.025,.5,.975))
colMeans(samples.beta.truncated)
hist(samples.beta.truncated$X1)
hist(samples.beta.truncated$X2)

### save plots
require(ggplot2)
pdf('MNL_Scenario1.lambda.pdf', width = 8, height = 8)
ggplot(data=data.frame(samples.lambda.truncated)) + geom_density(aes(x=samples.lambda.truncated), color="black") + scale_x_continuous(limits=c(90, 110))
dev.off()
pdf('MNL_Scenario1.beta1.pdf', width = 8, height = 8)
ggplot(data=samples.beta.truncated) + geom_density(aes(x=X1), color="black") + scale_x_continuous(limits=c(0, 0.15))
dev.off()
pdf('MNL_Scenario1.beta2.pdf', width = 8, height = 8)
ggplot(data=samples.beta.truncated) + geom_density(aes(x=X2), color="black") + scale_x_continuous(limits=c(0, 0.1))
dev.off()



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



