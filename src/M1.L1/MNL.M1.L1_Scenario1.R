####
# Scenario 1. No-purchase is not observed 
#
####

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src/M1.L1")

require("mvtnorm")
source(file="Metropolis-Hastings.R")



scenarioName <- "MNL.M1.L1_Scenario1"

## Load simulated choice data (NEED TO RUN MNL_InitData.binary.L1.R TO GENERATE THE DATA FIRST)
load(file="MNL.M1.L1_InitData.RData")

# final observation consists of the Demand
observation1 <- list(demand=Demand)



# log-posterior of betaT, to be called by M-H algorithm
logpost.betaT <- function(betaT, data) {
    
    beta <- c(betaT[1], betaT[2:(L+M)]*betaT[1])
    beta.coef <- beta[1:L]
    beta.const <- beta[(L+1):(L+M)]
    
    score <- rbind(t(exp(beta.const + X_Mat %*% beta.coef)), 1)
    choice.prob <- apply(score, 2, function(x) x/sum(x))
    
    logLikelihood <- data*log(choice.prob)
    
    logprior <- dmvnorm(beta, mean=beta.mu, sigma=beta.sg, log=TRUE)
    
    return(sum(logLikelihood) + logprior)
    
}


logpost.nopurchase <- function(nopurchase, demand, lambda, beta, k) {
    
    if (any(nopurchase<0))
        return(-Inf)
    else {
        beta.coef <- beta[1:L]
        beta.const <- beta[(L+1):(L+M)]
        
        score <- rbind(exp(beta.const + X_Mat[k,] %*% beta.coef), 1)
        choice.prob <- score/sum(score)
        
        dd <- c(demand, nopurchase)
        nn <- sum(dd)
        
        
        return(dmultinom(x=dd, prob=choice.prob, log=TRUE) + dpois(nn, lambda, log=TRUE))
        
    }
}


sample = function(data, parameters, nrun=1000) {

    demand <- data$demand
    
    # initialize arrays to save samples
    lambdas <- array(0, dim=c(nrun, 1))
    betas <- array(0, dim=c(nrun, L+M))
    nopurchases <- array(0, dim=c(nrun, K))
    
    # initial parameters
    #lambda1 <- parameters$lambda
    beta1 <- parameters$beta
    nopurchase1 <- parameters$nopurchase
    
    
    # start sampling loop
    for(i in 1:nrun) {

        #update posterior of lambda by conjugacy
        alpha2 <- lambda.alpha + sum(demand) + sum(nopurchase1)
        beta2 <- lambda.beta + K
        
        #simulate lambda2
        lambda2 <- rgamma(1, shape=alpha2, rate=beta2)
        
        
        
        #simulate beta2 by Metropolis-Hastings
        betaT1 <- c(beta1[1], beta1[2:(L+M)]/beta1[1])
        MH <- MH.mvnorm(logpost.betaT, start=betaT1, scale=c(0.2, 0.02), nrun=10, data=rbind(demand, nopurchase1))
        betaT2 <- MH$MC[10,]
        beta2 <- c(betaT2[1], betaT2[2:(L+M)]*betaT2[1])
        cat("MH acceptance rate: ", MH$accept, "\n")
        
        
        
        #simulate nopurchase by discrete Metropolis-Hastings
        nopurchase2 <- nopurchase1
        dMH.accept <- rep(0, K)
        for (j in 1:K) {
            dMH <- discreteMH.mvnorm(logpost.nopurchase, start=nopurchase1[j], scale=10, nrun=10, 
                              demand=demand[j], lambda=lambda2, beta=beta2, k=j)
            nopurchase2[j] <- dMH$MC[10]
            dMH.accept[j] <- dMH$accept
        }
        cat("discreteMH acceptance rate was ", mean(dMH.accept), "\n\n")
        
        
        
        # save the samples obtained in the current iteration
        lambdas[i,] <- lambda2
        betas[i,] <- beta2
        nopurchases[i,] <- nopurchase2
        
        cat("Run:", i, "\tlambda:", lambda2, "\tbeta2:", beta2, "\tnopurchase[1]:", nopurchase2[1], "\n", sep=" ")
        
        lambda1 <- lambda2
        beta1 <- beta2
        nopurchase1 <- nopurchase2
    }
    
    # results
    return(list(lambdas=lambdas, betas=betas, nopurchases=nopurchases))

} # end function 'sample'




### initialize input before sampling
# beta prior ~ N(beta.mu, beta.sg)
beta.mu <- rep(0, L+M)
beta.sg <- 100*diag(L+M)

# lambda ~ Gamma(lambda.alpha, lambda.beta)
lambda.alpha <- 0.005
lambda.beta <- 0.0001


## initial sampling input
param0 <- list(beta=c(-1, 1), nopurchase=rep(10, K))
nrun <- 5000
burnin <- 0.5
start <- burnin*nrun+1



### sample
z1 <- sample(data=observation1, parameters=param0, nrun=nrun)

save(z1, observation1, file=paste0(scenarioName, ".RData"))




### Visualize results

#plot lambda
samples.lambda <- z1$lambdas
plot(samples.lambda, type="l")

samples.lambda.truncated <- samples.lambda[start:nrun,]
mean(samples.lambda.truncated)
quantile(samples.lambda.truncated, c(.025,.5,.975))
hist(samples.lambda.truncated)


#plot beta
samples.beta <- data.frame(z1$betas)
plot(samples.beta$X1, type="l")
plot(samples.beta$X2, type="l")

samples.beta.truncated <- samples.beta[start:nrun,]
colMeans(samples.beta.truncated)
quantile(samples.beta.truncated$X1, c(.025,.5,.975))
quantile(samples.beta.truncated$X2, c(.025,.5,.975))
hist(samples.beta.truncated$X1)
hist(samples.beta.truncated$X2)


### save plots
require(ggplot2)

pdf(paste0(scenarioName, ".lambda.pdf"), width = 8, height = 8)
ggplot(data=data.frame(samples.lambda.truncated)) + geom_density(aes(x=samples.lambda.truncated), color="black")
dev.off()

pdf(paste0(scenarioName, ".beta.coef.pdf"), width = 8, height = 8)
ggplot(data=samples.beta.truncated) + geom_density(aes(x=X1), color="black")
dev.off()

pdf(paste0(scenarioName, ".beta.const.pdf"), width = 8, height = 8)
ggplot(data=samples.beta.truncated) + geom_density(aes(x=X2), color="black")
dev.off()


#plot nopurchase[1]
samples.nopurchase <- data.frame(z1$nopurchases)
plot(samples.nopurchase$X1, type="l")
plot(samples.nopurchase$X2, type="l")
plot(samples.nopurchase$X3, type="l")

samples.nopurchase.truncated <- samples.nopurchase[start:nrun,]
colMeans(samples.nopurchase.truncated)
quantile(samples.nopurchase.truncated$X1, c(.025,.5,.975))
quantile(samples.nopurchase.truncated$X2, c(.025,.5,.975))
quantile(samples.nopurchase.truncated$X3, c(.025,.5,.975))
hist(samples.nopurchase.truncated$X1)
hist(samples.nopurchase.truncated$X2)
hist(samples.nopurchase.truncated$X3)



