####
# Scenario 3. A noisy observation of No-purchase is observed 
# --- Noisy observation (say, traffic flow) is modeled as T = exp^epsilon1 * N, where N is total realized number of potential customers and epsilon1 and epsilon2 are noise terms ~ N(epsilon1.mean, epsilon1.sd) and N(epsilon2.mean, epsilon2.sd)
# --- Case 1: only multiplicative noise (only epsilon1, epsilon2=-Inf)
####

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src/binary.L1")

require("mvtnorm")
source(file="Metropolis-Hastings.R")



## Load simulated choice data (NEED TO RUN MNL_InitData.binary.R TO GENERATE THE DATA FIRST)
load(file="MNL_InitData.binary.L1.RData")



# parameters for the noise terms epsilon1
epsilon1.mean <- 3
epsilon1.sd <- 0.5 #unknown variance to be estimated


### simulate traffic observation
epsilon1 <- rnorm(K, mean=epsilon1.mean, sd=epsilon1.sd)
Traffic <- exp(epsilon1) * colSums(choice.mat) 

# final observation consists of the demand of alternative 1 and 2 and the traffic flow
observation3m <- list(demand=Demand, traffic=Traffic)



#log-posterior of beta, to be called by M-H algorithm
logpost.beta <- function(beta, data) {
    
    score <- rbind(t(exp(X_Mat %*% beta)), 1)
    choice.prob <- apply(score, 2, function(x) x/sum(x)) # M*N matrix
    
    logLikelihood <- data*log(choice.prob)
    
    logprior <- dmvnorm(beta, mean=beta.mu, sigma=beta.sg, log=TRUE)
    
    return(sum(logLikelihood) + logprior)
    
}


logpost.d0 <- function(d0, data, lambda, beta, k, eps1.mu, eps1.sd) {
    
    if (any(d0<0))
        return(-Inf)
    else {
        score <- rbind(t(exp(X_Mat[k,] %*% beta)), 1)
        choice.prob <- score/sum(score) # M*N matrix
        
        dd <- c(data$demand,d0)
        nn <- sum(dd)
        
        return(dmultinom(x=dd, prob=choice.prob, log=TRUE) + dpois(nn, lambda, log=TRUE) + dnorm(x=log(data$traffic/nn), mean=eps1.mu, sd=eps1.sd, log=TRUE))
    }
}



sample = function(data, parameters, nrun=1000) {
    
    demand <- data$demand
    traffic <- data$traffic

    
    # initialize arrays to save samples
    lambdas = array(0, dim=c(nrun, 1))
    d0s = array(0, dim=c(nrun, K))
    betas = array(0, dim=c(nrun, L))
    eps1.mus = array(0, dim=c(nrun, 1))
    eps1.sds = array(0, dim=c(nrun, 1))
    
    # initial parameters
    #lambda1 <- parameters$lambda
    d01 <- parameters$d0
    beta1 <- parameters$beta
    #eps1.mu1 <- parameters$eps1.mu
    eps1.sd1 <- parameters$eps1.sd
    
    
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



        #update and simulate eps1.mu and eps1.sd
        eps1 <- log(traffic) - log(demand+d01)
        #eps1.mu2 <- rnorm(1, mean=mean(eps1), sd=eps1.sd1/sqrt(K))
        eps1.mu2 <- epsilon1.mean    #eps1.mu is known, no updating
        eps1.pr2 <- rgamma(1, shape=eps1.pr.alpha0+K/2, rate=eps1.pr.beta0+sum((eps1-eps1.mu2)^2)/2)
        eps1.sd2 <- 1/sqrt(eps1.pr2)
        
        
        
        #simulate d0 by discrete Metropolis-Hastings
        d02 <- d01
        d0.accept <- rep(0, K)
        for (j in 1:K) {
            dataj <- list(demand=demand[j], traffic=traffic[j])
            sim <- discreteMH.norm(logpost.d0, start=d01[j], scale=10, nrun=10, 
                                   data=dataj, lambda=lambda2, beta=beta2, k=j, eps1.mu=eps1.mu2, eps1.sd=eps1.sd2)
            d02[j] <- sim$MC[10]
            d0.accept[j] <- sim$accept
        }
        cat("discreteMH acceptance rate was ", mean(d0.accept), "\n\n")
        

        # save the samples obtained in the current iteration
        lambdas[i,] = lambda2
        d0s[i,] = d02
        betas[i,] = beta2
        eps1.mus[i,] = eps1.mu2
        eps1.sds[i,] = eps1.sd2

        cat("Run:", i, "\tlambda:", lambda2, "\tbeta2:", beta2, "\teps1.mu2:", eps1.mu2, "\teps1.sd2:", eps1.sd2, "\td0[1]:", d02[1], "\n", sep=" ")
        
        lambda1 <- lambda2
        d01 <- d02
        beta1 <- beta2
        eps1.mu1 <- eps1.mu2
        eps1.sd1 <- eps1.sd2
    }
    
    # results
    return(list(lambdas=lambdas, d0s=d0s, betas=betas, eps1.mus=eps1.mus, eps1.sds=eps1.sds))

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
param0 <- list(beta=rep(0.1, L), d0=rep(10,K), eps1.sd=eps1.pr.alpha0/eps1.pr.beta0)
nrun <- 5000
burnin <- 0.5



### sample
z3m <- sample(data=observation3m, parameters=param0, nrun=nrun)

save(z3m, observation3m, file="MNL_Scenario3_M.binary.L1.RData")




### Visualize results
start <- burnin*nrun+1

#plot lambda
samples.lambda <- z3m$lambdas
plot(samples.lambda, type="l")


samples.lambda.truncated <- samples.lambda[start:nrun,]
quantile(samples.lambda.truncated, c(.025,.5,.975))
mean(samples.lambda.truncated)
hist(samples.lambda.truncated)



#plot beta
samples.beta <- z3m$betas
plot(samples.beta, type="l")



samples.beta.truncated <- samples.beta[start:nrun,]
quantile(samples.beta.truncated, c(.025,.5,.975))
mean(samples.beta.truncated)
hist(samples.beta.truncated)


### save plots
require(ggplot2)
pdf('MNL_Scenario3.binary.L1.lambda.pdf', width = 8, height = 8)
ggplot(data=data.frame(samples.lambda.truncated)) + geom_density(aes(x=samples.lambda.truncated), color="black")
dev.off()
pdf('MNL_Scenario3.binary.L1.beta.pdf', width = 8, height = 8)
ggplot(data=data.frame(samples.beta.truncated)) + geom_density(aes(x=samples.beta.truncated), color="black")
dev.off()



#plot d0[1]
samples.d0 <- data.frame(z3m$d0s)
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



#plot eps1.mu and eps1.sd
samples.eps1 <- data.frame(mu=z3m$eps1.mus, sd=z3m$eps1.sds)
plot(samples.eps1$mu, type="l")
plot(samples.eps1$sd, type="l")

samples.eps1.truncated <- samples.eps1[start:nrun,]
quantile(samples.eps1.truncated$mu, c(.025,.5,.975))
quantile(samples.eps1.truncated$sd, c(.025,.5,.975))
colMeans(samples.eps1.truncated)
hist(samples.eps1.truncated$mu)
hist(samples.eps1.truncated$sd)

