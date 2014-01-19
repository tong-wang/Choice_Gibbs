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
K <- 60 # number of periods

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


logpost.d0 <- function(d0, data, lambda, beta) {
    
    score <- apply(X_Mat, 1, function (x) exp(matrix(c(x,0,0), nrow=M, ncol=L, byrow=TRUE) %*% beta))
    choice.prob <- apply(score, 2, function(x) x/sum(x)) # M*N matrix
    
    dd <- rbind(data,d0)
    nn <- as.integer(colSums(dd))

    choice.ll <- rep(0,K)
    for (k in 1:K) {
        choice.ll[k] <- dmultinom(x=dd[, k], prob=choice.prob[, k], log=TRUE)    
    }    

    logLikelihood <- sum(choice.ll) + sum(dpois(nn, lambda, log=TRUE))
    
    return(sum(logLikelihood))
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
                         thin=1, mcmc=1, burnin=0, tune=2,
                         verbose=0, logfun=TRUE)[1,]


        #simulate d0 by Metropolis-Hastings
        ######NEED DISCRETE M-H ALGORITHM HERE!!!!#######
        d02 <- MCMCmetrop1R(logpost.d0, theta.init=d01,
                              data=data, lambda=lambda2, beta=beta2,
                              thin=1, mcmc=1, burnin=0, tune=2,
                              verbose=0, logfun=TRUE)[1,]
        
        
        

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
    return(list(bs=bs, Ws=Ws, betas=betas))

} # end function 'sample'




### initialize input before sampling
# b has a diffuse prior

# lambda ~ Gamma(lambda.alpha, lambda.beta)
lambda.alpha <- 0.5
lambda.beta <- 0.01


## initial sampling input
param0 <- list(beta=rep(0, L), lambda=lambda.alpha/lambda.beta, d0=rep(0,K))
nrun <- 1000



### sample
z <- sample(data=observation1, parameters=param0, nrun=nrun)



stopCluster(cl)


### test results
b2 <- c(-1.835909, 1.347521)
W2 <-  matrix(c(0.1995626, 0.2219391, 0.2219391, 0.2571687), nrow=L, ncol=L)
beta2 <- t(rmvnorm(N, mean=b2, sigma=W2)) # L*N matrix
score2 <- exp(XMAT %*% beta2) # M*N matrix
choice.prob2 <- apply(score2, 2, function(x) x/sum(x)) # M*N matrix
choice.mat2 <- apply(choice.prob2, 2, function(x) rmultinom(1, 1, x)) # M*N matrix
rowMeans(choice.mat2)
rowMeans(choice.mat)




### Visualize results
burnin <- 0.1*nrun

#plot b
#samples.b <- data.frame(matrix(unlist(z$bs), ncol=ncol(beta), byrow=T))
samples.b <- data.frame(z$bs)
plot(samples.b$X1, type="l")
plot(samples.b$X2, type="l")

plot(c(1,nrun), c(-2, 4), type="n", xlab="Samples", ylab="a", xaxt="n", yaxt="n")
lines(samples.b$X1, col="GREEN")
lines(samples.b$X2, col="RED")
axis(side=1)
axis(side=2)

samples.b.truncated <- samples.b[burnin:nrun,]
quantile(samples.b.truncated$X1, c(.025,.5,.975))
quantile(samples.b.truncated$X2, c(.025,.5,.975))
colMeans(samples.b.truncated)
hist(samples.b.truncated$X1)
hist(samples.b.truncated$X2)


#plot W
#samples.W <- data.frame(matrix(unlist(z$Ws), ncol=4, byrow=T))
samples.W <- data.frame(matrix(z$Ws, ncol=4))
plot(samples.W$X1, type="l")
plot(samples.W$X4, type="l")
plot(samples.W$X2, type="l")

plot(c(1,nrun), c(0, 5), type="n", xlab="Samples", ylab="a", xaxt="n", yaxt="n")
lines(samples.W$X1, col="GREEN")
lines(samples.W$X4, col="RED")
lines(samples.W$X2, col="BLUE")
axis(side=1)
axis(side=2)

samples.W.truncated <- samples.W[burnin:nrun,]
quantile(samples.W.truncated$X1, c(.025,.5,.975))
quantile(samples.W.truncated$X2, c(.025,.5,.975))
quantile(samples.W.truncated$X4, c(.025,.5,.975))
colMeans(samples.W.truncated)
hist(samples.W.truncated$X1)
hist(samples.W.truncated$X4)
hist(samples.W.truncated$X2)




#plot beta
samples.beta4 <- data.frame(z$betas[,,4])
plot(samples.beta4$X1, type="l")
plot(samples.beta4$X2, type="l")

plot(c(1,nrun), c(-5, 5), type="n", xlab="Samples", ylab="a", xaxt="n", yaxt="n")
lines(samples.beta4$X1, col="GREEN")
lines(samples.beta4$X2, col="BLUE")
axis(side=1)
axis(side=2)

samples.beta4.truncated <- samples.beta4[burnin:nrun,]
quantile(samples.beta4.truncated$X1, c(.025,.5,.975))
quantile(samples.beta4.truncated$X2, c(.025,.5,.975))
colMeans(samples.beta4.truncated)
hist(samples.beta4.truncated$X1)
hist(samples.beta4.truncated$X2)

