####
# Scenario 3c. A noisy observation of No-purchase is observed 
# --- Noisy observation (say, traffic flow) is modeled as T = exp^epsilon1 * N + exp^epsilon2, where N is total realized number of potential customers and epsilon1 and epsilon2 are noise terms ~ N(epsilon1.mean, epsilon1.sd) and N(epsilon2.mean, epsilon2.sd)
# --- Case 1: only multiplicative noise (only epsilon1, epsilon2=-Inf)
#    --- but demand is censored by inventory
####

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src/binary.L1")

require("mvtnorm")
source(file="Metropolis-Hastings.R")



## Load simulated choice data (NEED TO RUN MNL_InitData.binary.L1.R TO GENERATE THE DATA FIRST)
load(file="MNL_InitData.binary.L1.RData")

# final observation consists of the Sales, Stockout status and the Traffic flow
observation3c.m <- list(sales=Sales, stockout=Stockout, traffic=Traffic.m)



#log-posterior of beta, to be called by M-H algorithm
logpost.beta <- function(beta, data) {
    
    score <- rbind(t(exp(X_Mat %*% beta)), 1)
    choice.prob <- apply(score, 2, function(x) x/sum(x)) # M*N matrix
    
    logLikelihood <- data*log(choice.prob)
    
    logprior <- dmvnorm(beta, mean=beta.mu, sigma=beta.sg, log=TRUE)
    
    return(sum(logLikelihood) + logprior)
    
}


logpost.nopurchase <- function(nopurchase, demand, traffic, lambda, beta, k, eps1.mu, eps1.sd) {
    
    if (any(nopurchase<0))
        return(-Inf)
    else {
        score <- rbind(t(exp(X_Mat[k,] %*% beta)), 1)
        choice.prob <- score/sum(score) # M*N matrix
        
        dd <- c(demand, nopurchase)
        nn <- sum(dd)
        
        return(dmultinom(x=dd, prob=choice.prob, log=TRUE) + dpois(nn, lambda, log=TRUE) + dnorm(x=log(traffic/nn), mean=eps1.mu, sd=eps1.sd, log=TRUE))
    }
}


logpost.demand_nopurchase <- function(x, sales, traffic, lambda, beta, k, eps1.mu, eps1.sd) {
    
    demand <- x[1]
    nopurchase <- x[2]
    
    if (any(nopurchase<0) | any(demand<sales))
        return(-Inf)
    else {
        score <- rbind(t(exp(X_Mat[k,] %*% beta)), 1)
        choice.prob <- score/sum(score) # M*N matrix
        
        dd <- c(demand, nopurchase)
        nn <- sum(dd)
        
        return(dmultinom(x=dd, prob=choice.prob, log=TRUE) + dpois(nn, lambda, log=TRUE) + dnorm(x=log(traffic/nn), mean=eps1.mu, sd=eps1.sd, log=TRUE))
    }
}


sample = function(data, parameters, nrun=1000) {
    
    sales <- data$sales
    stockout <- data$stockout
    traffic <- data$traffic

    
    # initialize arrays to save samples
    lambdas <- array(0, dim=c(nrun, 1))
    betas <- array(0, dim=c(nrun, L))
    eps1.mus <- array(0, dim=c(nrun, 1))
    eps1.sds <- array(0, dim=c(nrun, 1))
    nopurchases <- array(0, dim=c(nrun, K))
    demands <- array(0, dim=c(nrun, K))
    
    # initial parameters
    lambda1 <- parameters$lambda
    beta1 <- parameters$beta
    #eps1.mu1 <- parameters$eps1.mu
    eps1.sd1 <- parameters$eps1.sd
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
        MH <- MH.mvnorm(logpost.beta, start=beta1, scale=0.1, nrun=10, data=rbind(demand1, nopurchase1))
        beta2 <- MH$MC[10,]
        cat("MH acceptance rate: ", MH$accept, "\n")



        #update and simulate eps1.mu and eps1.sd
        eps1 <- log(traffic) - log(demand1 + nopurchase1)
        #eps1.mu2 <- rnorm(1, mean=mean(eps1), sd=eps1.sd1/sqrt(K))
        eps1.mu2 <- epsilon1.mean    #eps1.mu is known, no updating
        eps1.pr2 <- rgamma(1, shape=eps1.pr.alpha0+K/2, rate=eps1.pr.beta0+sum((eps1-eps1.mu2)^2)/2)
        eps1.sd2 <- 1/sqrt(eps1.pr2)
        
        
        
        #simulate nopurchase by discrete Metropolis-Hastings
        nopurchase2 <- nopurchase1
        demand2 <- demand1
        dMH.accept <- rep(NA, K)
        dMH2.accept <- rep(NA, K)
        for (j in 1:K) {
            if (!stockout[j]) {
                
                dMH <- discreteMH.mvnorm(logpost.nopurchase, start=nopurchase1[j], scale=10, nrun=10, 
                                       demand=demand1[j], traffic=traffic[j], lambda=lambda2, beta=beta2, k=j, eps1.mu=eps1.mu2, eps1.sd=eps1.sd2)
                nopurchase2[j] <- dMH$MC[10]
                dMH.accept[j] <- dMH$accept
                
            } else {
                
                dMH2 <- discreteMH.mvnorm(logpost.demand_nopurchase, start=c(demand1[j], nopurchase1[j]), scale=5, nrun=10, 
                                          sales=sales[j], traffic=traffic[j], lambda=lambda2, beta=beta2, k=j, eps1.mu=eps1.mu2, eps1.sd=eps1.sd2)
                demand2[j] <- dMH2$MC[10,1]
                nopurchase2[j] <- dMH2$MC[10,2]
                dMH2.accept[j] <- dMH2$accept
            }
            
        }
        cat("discreteMH acceptance rate was ", mean(dMH.accept, na.rm=TRUE), "\t", mean(dMH2.accept, na.rm=TRUE), "\n\n")
        
        
        
        # save the samples obtained in the current iteration
        lambdas[i,] <- lambda2
        betas[i,] <- beta2
        eps1.mus[i,] <- eps1.mu2
        eps1.sds[i,] <- eps1.sd2
        nopurchases[i,] <- nopurchase2
        demands[i,] <- demand2
        
        cat("Run:", i, "\tlambda:", lambda2, "\tbeta2:", beta2, "\teps1.mu2:", eps1.mu2, "\teps1.sd2:", eps1.sd2, "\tnopurchase[1]:", nopurchase2[1], "\td[1]:", demand2[1], "\n", sep=" ")
        
        lambda1 <- lambda2
        beta1 <- beta2
        eps1.mu1 <- eps1.mu2
        eps1.sd1 <- eps1.sd2
        nopurchase1 <- nopurchase2
        demand1 <- demand2
    }
    
    # results
    return(list(lambdas=lambdas, betas=betas, eps1.mus=eps1.mus, eps1.sds=eps1.sds, nopurchases=nopurchases, demands=demands))

} # end function 'sample'




### initialize input before sampling
# beta prior ~ N(beta.mu, beta.sg)
beta.mu <- rep(0, L)
beta.sg <- 100*diag(L)

# lambda ~ Gamma(lambda.alpha, lambda.beta)
lambda.alpha <- 0.005
lambda.beta <- 0.0001

# prior of the precision of epsilon
eps1.pr.alpha0 <- 0.0000001 #0.0001
eps1.pr.beta0 <- 0.000000001 #0.0001

## initial sampling input
param0 <- list(lambda=lambda.alpha/lambda.beta, beta=rep(0.1, L), eps1.sd=eps1.pr.alpha0/eps1.pr.beta0, nopurchase=rep(10,K))
nrun <- 5000
burnin <- 0.5
start <- burnin*nrun+1



### sample
z3c.m <- sample(data=observation3c.m, parameters=param0, nrun=nrun)

save(z3c.m, observation3c.m, file="MNL_Scenario3c_M.binary.L1.m.RData")




### Visualize results

#plot lambda
samples.lambda <- z3c.m$lambdas
plot(samples.lambda, type="l")

samples.lambda.truncated <- samples.lambda[start:nrun,]
quantile(samples.lambda.truncated, c(.025,.5,.975))
mean(samples.lambda.truncated)
hist(samples.lambda.truncated)


#plot beta
samples.beta <- z3c.m$betas
plot(samples.beta, type="l")

samples.beta.truncated <- samples.beta[start:nrun,]
quantile(samples.beta.truncated, c(.025,.5,.975))
mean(samples.beta.truncated)
hist(samples.beta.truncated)


### save plots
require(ggplot2)
pdf('MNL_Scenario3c_M.binary.L1.lambda.m.pdf', width = 8, height = 8)
ggplot(data=data.frame(samples.lambda.truncated)) + geom_density(aes(x=samples.lambda.truncated), color="black")
dev.off()
pdf('MNL_Scenario3c_M.binary.L1.beta.m.pdf', width = 8, height = 8)
ggplot(data=data.frame(samples.beta.truncated)) + geom_density(aes(x=samples.beta.truncated), color="black")
dev.off()


#plot eps1.mu and eps1.sd
samples.eps1 <- data.frame(mu=z3c.m$eps1.mus, sd=z3c.m$eps1.sds)
plot(samples.eps1$mu, type="l")
plot(samples.eps1$sd, type="l")

samples.eps1.truncated <- samples.eps1[start:nrun,]
quantile(samples.eps1.truncated$mu, c(.025,.5,.975))
quantile(samples.eps1.truncated$sd, c(.025,.5,.975))
colMeans(samples.eps1.truncated)
hist(samples.eps1.truncated$mu)
hist(samples.eps1.truncated$sd)


#plot nopurchase[1]
samples.nopurchase <- data.frame(z3c.m$nopurchases)
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
samples.demand <- data.frame(z3c.m$demands)
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


