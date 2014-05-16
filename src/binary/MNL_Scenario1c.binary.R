####
# Scenario 1c. No-purchase is not observed 
#    --- but demand is censored by inventory
####

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src/binary")

require("mvtnorm")
source(file="Metropolis-Hastings.R")



## Load simulated choice data (NEED TO RUN MNL_InitData.binary.R TO GENERATE THE DATA FIRST)
load(file="MNL_InitData.binary.RData")

# final observation consists of the Sales and Stockout status
observation1c <- list(sales=Sales, stockout=Stockout)



#log-posterior of beta, to be called by M-H algorithm
logpost.beta <- function(beta, data) {
    
    score <- rbind(t(exp(X_Mat %*% beta)), 1)
    choice.prob <- apply(score, 2, function(x) x/sum(x)) # M*N matrix
    
    logLikelihood <- data*log(choice.prob)
    
    logprior <- dmvnorm(beta, mean=beta.mu, sigma=beta.sg, log=TRUE)
    
    return(sum(logLikelihood) + logprior)
    
}


logpost.nopurchase <- function(nopurchase, demand, lambda, beta, k) {
    
    if (any(nopurchase<0))
        return(-Inf)
    else {
        score <- rbind(t(exp(X_Mat[k,] %*% beta)), 1)
        choice.prob <- score/sum(score) # M*N matrix
        
        dd <- c(demand, nopurchase)
        nn <- sum(dd)
        
        
        return(dmultinom(x=dd, prob=choice.prob, log=TRUE) + dpois(nn, lambda, log=TRUE))
        
    }
}


logpost.demand_nopurchase <- function(x, sales, lambda, beta, k) {
    
    demand <- x[1]
    nopurchase <- x[2]
    
    if (any(nopurchase<0) | any(demand<sales))
        return(-Inf)
    else {
        score <- rbind(t(exp(X_Mat[k,] %*% beta)), 1)
        choice.prob <- score/sum(score) # M*N matrix
        
        dd <- c(demand, nopurchase)
        nn <- sum(dd)
        
        return(dmultinom(x=dd, prob=choice.prob, log=TRUE) + dpois(nn, lambda, log=TRUE))
    }
}



sample = function(data, parameters, nrun=1000) {

    sales <- data$sales
    stockout <- data$stockout
    
    # initialize arrays to save samples
    lambdas <- array(0, dim=c(nrun, 1))
    betas <- array(0, dim=c(nrun, L))
    nopurchases <- array(0, dim=c(nrun, K))
    demands <- array(0, dim=c(nrun, K))
    
    # initial parameters
    lambda1 <- parameters$lambda
    beta1 <- parameters$beta
    nopurchase1 <- parameters$nopurchase
    demand1 <- sales
    
    
    # start sampling loop
    for(i in 1:nrun) {
        
        #update posterior of lambda by conjugacy
        alpha2 <- lambda.alpha + sum(demand1) + sum(nopurchase1)
        beta2 <- lambda.beta + K
        
        #simulate lambda2
        lambda2 <- rgamma(1, shape=alpha2, rate=beta2)
        
        
        
        #simulate beta2 by Metropolis-Hastings
        MH <- MH.mvnorm(logpost.beta, start=beta1, scale=c(0.06, 0.03), nrun=10, data=rbind(demand1, nopurchase1))
        beta2 <- MH$MC[10,]
        cat("MH acceptance rate: ", MH$accept, "\n")
        
        
        
        #simulate nopurchase by discrete Metropolis-Hastings
        nopurchase2 <- nopurchase1
        demand2 <- demand1
        dMH.accept <- rep(NA, K)
        dMH2.accept <- rep(NA, K)
        for (j in 1:K) {
            if (!stockout[j]) {
                
                dMH <- discreteMH.mvnorm(logpost.nopurchase, start=nopurchase1[j], scale=10, nrun=10, 
                                         demand=demand1[j], lambda=lambda2, beta=beta2, k=j)
                nopurchase2[j] <- dMH$MC[10]
                dMH.accept[j] <- dMH$accept
                
            } else {
                
                dMH2 <- discreteMH.mvnorm(logpost.demand_nopurchase, start=c(demand1[j], nopurchase1[j]), scale=5, nrun=10, 
                                          sales=sales[j], lambda=lambda2, beta=beta2, k=j)
                demand2[j] <- dMH2$MC[10,1]
                nopurchase2[j] <- dMH2$MC[10,2]
                dMH2.accept[j] <- dMH2$accept
            }
            
        }
        cat("discreteMH acceptance rate was ", mean(dMH.accept, na.rm=TRUE), "\t", mean(dMH2.accept, na.rm=TRUE), "\n\n")
        
        
        
        # save the samples obtained in the current iteration
        lambdas[i,] <- lambda2
        betas[i,] <- beta2
        nopurchases[i,] <- nopurchase2
        demands[i,] <- demand2
        
        cat("Run:", i, "\tlambda:", lambda2, "\tbeta2:", beta2, "\tnopurchase[1]:", nopurchase2[1], "\td[1]:", demand2[1], "\n", sep=" ")
        
        lambda1 <- lambda2
        beta1 <- beta2
        nopurchase1 <- nopurchase2
        demand1 <- demand2
    }
    
    # results
    return(list(lambdas=lambdas, betas=betas, nopurchases=nopurchases, demands=demands))

} # end function 'sample'




### initialize input before sampling
# beta prior ~ N(beta.mu, beta.sg)
beta.mu <- rep(0, L)
beta.sg <- 100*diag(L)

# lambda ~ Gamma(lambda.alpha, lambda.beta)
lambda.alpha <- 0.005
lambda.beta <- 0.0001


## initial sampling input
param0 <- list(lambda=lambda.alpha/lambda.beta, beta=rep(0.1,L), nopurchase=rep(10,K))
nrun <- 5000
burnin <- 0.5
start <- burnin*nrun+1



### sample
z1c <- sample(data=observation1c, parameters=param0, nrun=nrun)


save(z1c, observation1c, file="MNL_Scenario1c.binary.RData")





### Visualize results

#plot lambda
samples.lambda <- z1c$lambdas
plot(samples.lambda, type="l")

samples.lambda.truncated <- samples.lambda[start:nrun,]
quantile(samples.lambda.truncated, c(.025,.5,.975))
mean(samples.lambda.truncated)
hist(samples.lambda.truncated)


#plot beta
samples.beta <- data.frame(z1c$betas)
plot(samples.beta$X1, type="l")
plot(samples.beta$X2, type="l")


samples.beta.truncated <- samples.beta[start:nrun,]
quantile(samples.beta.truncated$X1, c(.025,.5,.975))
quantile(samples.beta.truncated$X2, c(.025,.5,.975))
colMeans(samples.beta.truncated)
hist(samples.beta.truncated$X1)
hist(samples.beta.truncated$X2)


### save plots
require(ggplot2)
pdf('MNL_Scenario1c.binary.lambda.pdf', width = 8, height = 8)
ggplot(data=data.frame(samples.lambda.truncated)) + geom_density(aes(x=samples.lambda.truncated), color="black")
dev.off()
pdf('MNL_Scenario1c.binary.beta1.pdf', width = 8, height = 8)
ggplot(data=samples.beta.truncated) + geom_density(aes(x=X1), color="black")
dev.off()
pdf('MNL_Scenario1c.binary.beta2.pdf', width = 8, height = 8)
ggplot(data=samples.beta.truncated) + geom_density(aes(x=X2), color="black")
dev.off()


#plot nopurchase[1]
samples.nopurchase <- data.frame(z1c$nopurchases)
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


#plot demand[1]
samples.demand <- data.frame(z1c$demands)
plot(samples.demand$X1, type="l")
plot(samples.demand$X2, type="l")
plot(samples.demand$X3, type="l")

samples.demand.truncated <- samples.demand[start:nrun,]
colMeans(samples.demand.truncated)
quantile(samples.demand.truncated$X1, c(.025,.5,.975))
quantile(samples.demand.truncated$X2, c(.025,.5,.975))
quantile(samples.demand.truncated$X3, c(.025,.5,.975))
hist(samples.demand.truncated$X1)
hist(samples.demand.truncated$X2)
hist(samples.demand.truncated$X3)

