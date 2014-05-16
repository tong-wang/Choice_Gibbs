####
# Scenario 0c. All choices are observed 
#    --- but demand is censored by inventory
####

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src/binary.L1")

require("mvtnorm")
source(file="Metropolis-Hastings.R")



## Load simulated choice data (NEED TO RUN MNL_InitData.binary.R TO GENERATE THE DATA FIRST)
load(file="MNL_InitData.binary.L1.RData")

# final observation consists of Sales, Stockout status, and Rest
observation0c <- list(sales=Sales, stockout=Stockout, rest=NoPurchase+LostSales)



#log-posterior of beta, to be called by M-H algorithm
logpost.beta <- function(beta, data) {
    
    score <- rbind(t(exp(X_Mat %*% beta)), 1)
    choice.prob <- apply(score, 2, function(x) x/sum(x)) # M*N matrix
    
    logLikelihood <- data*log(choice.prob)
    
    logprior <- dmvnorm(beta, mean=beta.mu, sigma=beta.sg, log=TRUE)
    
    return(sum(logLikelihood) + logprior)
    
}


logpost.nopurchase <- function(nopurchase, sales, rest, lambda, beta, k) {

    if (any(nopurchase<0) | any(nopurchase>rest))
        return(-Inf)
    else {
        score <- rbind(t(exp(X_Mat[k,] %*% beta)), 1)
        choice.prob <- score/sum(score) # M*N matrix
        
        dd <- c(sales + rest - nopurchase, nopurchase)
        nn <- sum(dd)
        
        
        return(dmultinom(x=dd, prob=choice.prob, log=TRUE) + dpois(nn, lambda, log=TRUE))
        
    }
}



sample = function(data, parameters, nrun=1000) {

    sales <- data$sales
    stockout <- data$stockout
    rest <- data$rest
    
    # initialize arrays to save samples
    lambdas <- array(0, dim=c(nrun, 1))
    betas <- array(0, dim=c(nrun, L))
    nopurchases <- array(0, dim=c(nrun, K))
    
    # initial parameters
    lambda1 <- parameters$lambda
    beta1 <- parameters$beta
    nopurchase1 <- rest
    
    
    # start sampling loop
    for(i in 1:nrun) {
        
        demand <- sales + rest - nopurchase1

        #update posterior of lambda by conjugacy
        alpha2 <- lambda.alpha + sum(demand) + sum(nopurchase1)
        beta2 <- lambda.beta + K
        
        #simulate lambda2
        lambda2 <- rgamma(1, shape=alpha2, rate=beta2)
        
        
        
        #simulate beta2 by Metropolis-Hastings
        MH <- MH.mvnorm(logpost.beta, start=beta1, scale=0.1, nrun=10, data=rbind(demand, nopurchase1))
        beta2 <- MH$MC[10,]
        cat("MH acceptance rate: ", MH$accept, "\n")
        
        
        
        #simulate nopurchase by discrete Metropolis-Hastings
        nopurchase2 <- nopurchase1
        dMH.accept <- rep(NA, K)
        for (j in which(stockout)) {
            dMH <- discreteMH.mvnorm(logpost.nopurchase, start=nopurchase1[j], scale=10, nrun=10, 
                          sales=sales[j], rest=rest[j], lambda=lambda2, beta=beta2, k=j)
            nopurchase2[j] <- dMH$MC[10]
            dMH.accept[j] <- dMH$accept
        }
        cat("discreteMH acceptance rate was ", mean(dMH.accept, na.rm=TRUE), "\n\n")
        
        
        
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
beta.mu <- rep(0, L)
beta.sg <- 100*diag(L)

# lambda ~ Gamma(lambda.alpha, lambda.beta)
lambda.alpha <- 0.005
lambda.beta <- 0.0001


## initial sampling input
param0 <- list(lambda=lambda.alpha/lambda.beta, beta=rep(0.1,L))
nrun <- 5000
burnin <- 0.5
start <- burnin*nrun+1



### sample
z0c <- sample(data=observation0c, parameters=param0, nrun=nrun)


save(z0c, observation0c, file="MNL_Scenario0c.binary.L1.RData")





### Visualize results

#plot lambda
samples.lambda <- z0c$lambdas
plot(samples.lambda, type="l")

samples.lambda.truncated <- samples.lambda[start:nrun,]
quantile(samples.lambda.truncated, c(.025,.5,.975))
mean(samples.lambda.truncated)
hist(samples.lambda.truncated)


#plot beta
samples.beta <- z0c$betas
plot(samples.beta, type="l")

samples.beta.truncated <- samples.beta[start:nrun,]
quantile(samples.beta.truncated, c(.025,.5,.975))
mean(samples.beta.truncated)
hist(samples.beta.truncated)


### save plots
require(ggplot2)
pdf('MNL_Scenario0c.binary.L1.lambda.pdf', width = 8, height = 8)
ggplot(data=data.frame(samples.lambda.truncated)) + geom_density(aes(x=samples.lambda.truncated), color="black")
dev.off()
pdf('MNL_Scenario0c.binary.L1.beta.pdf', width = 8, height = 8)
ggplot(data=data.frame(samples.beta.truncated)) + geom_density(aes(x=samples.beta.truncated), color="black")
dev.off()


#plot nopurchase[1]
samples.nopurchase <- data.frame(z0c$nopurchases)
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

