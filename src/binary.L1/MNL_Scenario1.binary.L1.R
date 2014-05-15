####
# Scenario 1. No-purchase is not observed 
#
####

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src/binary.L1")

require("mvtnorm")
source(file="Metropolis-Hastings.R")



## Load simulated choice data (NEED TO RUN MNL_InitData.binary.R TO GENERATE THE DATA FIRST)
load(file="MNL_InitData.binary.L1.RData")



observation1 <- list(demand=Demand)


#log-posterior of beta, to be called by M-H algorithm
logpost.beta <- function(beta, data) {
    
    score <- rbind(t(exp(X_Mat %*% beta)), 1)
    choice.prob <- apply(score, 2, function(x) x/sum(x)) # M*N matrix
    
    logLikelihood <- data*log(choice.prob)
    
    logprior <- dmvnorm(beta, mean=beta.mu, sigma=beta.sg, log=TRUE)
    
    return(sum(logLikelihood) + logprior)
    
}



logpost.d0 <- function(d0, data, lambda, beta, k) {
    
    if (any(d0<0))
        return(-Inf)
    else {
        score <- rbind(t(exp(X_Mat[k,] %*% beta)), 1)
        choice.prob <- score/sum(score) # M*N matrix
        
        dd <- c(data,d0)
        nn <- sum(dd)
        
        
        return(dmultinom(x=dd, prob=choice.prob, log=TRUE) + dpois(nn, lambda, log=TRUE))
        
    }
}



sample = function(data, parameters, nrun=1000) {

    demand = data$demand
    
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
        alpha2 <- lambda.alpha + sum(demand) + sum(d01)
        beta2 <- lambda.beta + K
        
        #simulate lambda2
        lambda2 <- rgamma(1, shape=alpha2, rate=beta2)
        
        
        
        #simulate beta2 by Metropolis-Hastings
        MH <- MH.mvnorm(logpost.beta, sigma=diag(L), scale=0.1, start=beta1, nrun = 10, data=rbind(demand,d01))
        beta2 <- MH$MC[10,]
        cat("MH acceptance rate: ", MH$accept, "\n")
        

        #simulate d0 by discrete Metropolis-Hastings
        d02 <- d01
        d0.accept <- rep(0, K)
        for (j in 1:K) {
            sim <- discreteMH.norm(logpost.d0, start=d01[j], scale=10, nrun=10, 
                              data=demand[j], lambda=lambda2, beta=beta2, k=j)
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
# beta prior ~ N(beta.mu, beta.sg)
beta.mu <- rep(0, L)
beta.sg <- 100*diag(L)

# lambda ~ Gamma(lambda.alpha, lambda.beta)
lambda.alpha <- 0.005
lambda.beta <- 0.0001


## initial sampling input
param0 <- list(beta=rep(0.1, L), lambda=lambda.alpha/lambda.beta, d0=rep(10,K))
nrun <- 5000
burnin <- 0.5



### sample
z1 <- sample(data=observation1, parameters=param0, nrun=nrun)


save(z1, observation1, file="MNL_Scenario1.binary.L1.RData")





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


### save plots
require(ggplot2)
pdf('MNL_Scenario1.binary.L1.lambda.pdf', width = 8, height = 8)
ggplot(data=data.frame(samples.lambda.truncated)) + geom_density(aes(x=samples.lambda.truncated), color="black")
dev.off()
pdf('MNL_Scenario1.binary.L1.beta.pdf', width = 8, height = 8)
ggplot(data=data.frame(samples.beta.truncated)) + geom_density(aes(x=samples.beta.truncated), color="black")
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



