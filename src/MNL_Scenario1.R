####
# Scenario 1. No-purchase is not observed 
#
####

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src")

require("MCMCpack")
require("mvtnorm")



### Known parameters
M <- 3 # number of alternatives (alternative 3 is dummy for no-purchase)
L <- 2 # number of covariates
K <- 360 # number of periods

#X is the attributes of the alternatives; in each period, [Xij] is an (M-1)*L matrix, i=1...M-1, j=1...L.
#by row: [X11 X12; X21 X22]
X_Mean <- c(4, 3, 6, 1)

#for each of the K periods, generate an X matrix
#dim of X_Mat is K, (M-1)*L
X_Mat <- rmvnorm(K, mean=X_Mean)


## true values of the parameters to be estimated
# beta is the MNL coefficient
beta <- c(0.06, 0.03); # L-dimensional
# lambda is the poisson demand rate in each perid
lambda <- 100


### simulate data
# simulate number of demand per period, ~Poisson(lambda)
N <- rpois(K, lambda) 

#in each period, score of choice 1 and 2 (col) is exp(X*beta)
#dim of score is M by K
score <- apply(X_Mat, 1, function (x) exp(matrix(c(x,0,0), nrow=M, ncol=L, byrow=TRUE) %*% beta))
#choice probabilities in each period (M by K)
choice.prob <- apply(score, 2, function(x) x/sum(x)) # M*N matrix
#simulate actual choice in each period (M by K)
choice.mat <- matrix(0, nrow=3, ncol=K)
for (k in 1:K) {
    choice.mat[, k] <- rmultinom(1, N[k], choice.prob[, k])    
}    

rowMeans(choice.mat)

observation1 <- choice.mat[1:2,]



#log-posterior of beta, to be called by M-H algorithm
logpost.beta <- function(beta, data) {
    
    score <- apply(X_Mat, 1, function (x) exp(matrix(c(x,0,0), nrow=M, ncol=L, byrow=TRUE) %*% beta))
    choice.prob <- apply(score, 2, function(x) x/sum(x)) # M*N matrix
    

    logLikelihood <- data*log(choice.prob)
    
    return(sum(logLikelihood))
}


logpost.d0 <- function(d0, data, lambda, beta, k) {
    
    score <-  exp(matrix(c(X_Mat[k,],0,0), nrow=M, ncol=L, byrow=TRUE) %*% beta)
    choice.prob <- score/sum(score) # M*N matrix
    
    dd <- c(data,d0)
    nn <- sum(dd)

    #choice.ll <- rep(0,K)
    #for (k in 1:K) {
    #    choice.ll[k] <- dmultinom(x=dd[, k], prob=choice.prob[, k], log=TRUE)    
    #}    
    return(dmultinom(x=dd, prob=choice.prob, log=TRUE) + dpois(nn, lambda, log=TRUE))
}

#implementation of discrete Metropolic-Hastings with Negative Binomial proposal
discreteMH <-function (logpost, proposal, start, m, ...) 
{
    pb = length(start)
    Mpar = array(0, c(m, pb))
    b = matrix(t(start))
    lb = logpost(start, ...)
    
    size = proposal$size
    
    accept = 0
    for (i in 1:m) {
        #proposed next point is Negative Binomial using the current position as its mean, variance is controled by $size$
        bc <- rnbinom(pb, size=size, mu=b)
        
        lbc = logpost(t(bc), ...) 
        #since proposal is asymmetric, need to adjust the jumping probability accordingly
        prob = exp(lbc + dnbinom(b, size=size, mu=bc, log=TRUE) - lb - dnbinom(bc, size=size, mu=b, log=TRUE) )
        
        if (is.na(prob) == FALSE) {
            if (runif(1) < prob) {
                lb = lbc
                b = bc
                accept = accept + 1
            }
        }
        Mpar[i, ] = b
    }
    accept = accept/m
    stuff = list(par = Mpar, accept = accept)
    return(stuff)
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
        lambda2 <- rpois(1, alpha2/beta2)
        
        
        
        #simulate beta2 by Metropolis-Hastings
        beta2 <- MCMCmetrop1R(logpost.beta, theta.init=beta1,
                         data=rbind(data,d01),
                         thin=1, mcmc=1, burnin=1000, tune=2,
                         verbose=0, logfun=TRUE)[1,]


        #simulate d0 by discrete Metropolis-Hastings
        d02 <- rep(0, K)
        for (j in 1:K) {
            d02[j] <- discreteMH(logpost.d0, proposal=list(size=10), start=d01[j], m=1000, 
                             data=data[,j], lambda=lambda2, beta=beta2, k=j)$par[1000]
        }
        
        

        # save the samples obtained in the current iteration
        lambdas[i,] = lambda2
        d0s[i,] = d02
        betas[i,] = beta2
        
        cat("Run:", i, "\tlambda:", lambda2, "\tbeta2:", beta2, "\n", sep=" ")
        
        lambda1 <- lambda2
        d01 <- d02
        beta1 <- beta2
    }
    
    # results
    return(list(lambdas=lambdas, d0s=d0s, betas=betas))

} # end function 'sample'




### initialize input before sampling
# b has a diffuse prior

# lambda ~ Gamma(lambda.alpha, lambda.beta)
lambda.alpha <- 0.5
lambda.beta <- 0.01


## initial sampling input
param0 <- list(beta=rep(0, L), lambda=lambda.alpha/lambda.beta, d0=rep(10,K))
nrun <- 10000



### sample
z <- sample(data=observation1, parameters=param0, nrun=nrun)



#stopCluster(cl)



### Visualize results
burnin <- 0.1*nrun

#plot lambda
samples.lambda <- z$lambdas
plot(samples.lambda, type="l")


samples.lambda.truncated <- samples.lambda[burnin:nrun,]
quantile(samples.lambda.truncated, c(.025,.5,.975))
mean(samples.lambda.truncated)
hist(samples.lambda.truncated)



#plot beta
samples.beta <- data.frame(z$betas)
plot(samples.beta$X1, type="l")
plot(samples.beta$X2, type="l")

plot(c(1,nrun), c(-0.05, 0.2), type="n", xlab="Samples", ylab="a", xaxt="n", yaxt="n")
lines(samples.beta$X1, col="GREEN")
lines(samples.beta$X2, col="BLUE")
axis(side=1)
axis(side=2)

samples.beta.truncated <- samples.beta[burnin:nrun,]
quantile(samples.beta.truncated$X1, c(.025,.5,.975))
quantile(samples.beta.truncated$X2, c(.025,.5,.975))
colMeans(samples.beta.truncated)
hist(samples.beta.truncated$X1)
hist(samples.beta.truncated$X2)



#plot d0[1]
samples.d0 <- data.frame(z$d0s)
plot(samples.d0$X1, type="l")
plot(samples.d0$X2, type="l")
plot(samples.d0$X3, type="l")
plot(samples.d0$X4, type="l")
plot(samples.d0$X5, type="l")


samples.d0.truncated <- samples.d0[burnin:nrun,]
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



