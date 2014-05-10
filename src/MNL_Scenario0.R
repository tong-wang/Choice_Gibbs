####
# Scenario 0. All choices are observed 
#
####

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src")

require("MCMCpack")
require("mvtnorm")


## Load simulated choice data (NEED TO RUN MNL_InitData.R TO GENERATE THE DATA FIRST)
load(file="MNL_InitData.RData")



observation0 <- choice.mat



### Estimating beta by M-H sampling
# log-posterior of beta, to be called by M-H algorithm
logpost.beta <- function(beta, data) {
    
    if (any(beta<0))
        return(-Inf)
    else {
        score <- rbind(apply(X_Mat, 1, function (x) exp(matrix(x, nrow=M-1, ncol=L, byrow=TRUE) %*% beta)), 1)
        choice.prob <- apply(score, 2, function(x) x/sum(x)) # M*N matrix
        
        logLikelihood <- data*log(choice.prob)
        
        logprior <- dmvnorm(log(beta), mean=beta.mu, sigma=beta.sg, log=TRUE)
        
        return(sum(logLikelihood) + logprior)
    }
}



## initialize input before sampling
# beta prior ~ logN(beta.mu, beta.sg)
beta.mu <- rep(-2, L)
beta.sg <- diag(10, nrow=L, ncol=L)

nrun <- 50000
burnin <- 0.5


#direct sampling
z0 <- MCMCmetrop1R(logpost.beta, theta.init=rep(0.1, L),
                  data=observation0,
                  thin=1, mcmc=nrun*(1-burnin), burnin=nrun*burnin, tune=0.01,
                  verbose=500,  V=diag(1,L,L))


summary(z0)
plot(z0)
beta.estm <- colMeans(z0)


### Estimate lambda by conjugate prior

# lambda ~ Gamma(lambda.alpha, lambda.beta)
lambda.alpha <- 0.005
lambda.beta <- 0.0001

#update posterior of lambda by conjugacy
alpha2 <- lambda.alpha + sum(observation0)
beta2 <- lambda.beta + K

lambda.estm <- alpha2/beta2

score.estm <- exp(matrix(X_Mean, M-1, L) %*% beta.estm)
score.estm <- c(score.estm, 1)
choice.prob.estm <- score.estm / sum(score.estm)
demand.estm <- lambda.estm * choice.prob.estm




z0 <- list(lambdas=rgamma(nrun*(1-burnin), shape=alpha2, rate=beta2), betas=z0[,])

save(z0, observation0, file="MNL_Scenario0.RData")



### save plots
require(ggplot2)

pdf('MNL_Scenario0.lambda.pdf', width = 8, height = 8)
ggplot(data=as.data.frame(z0$lambdas)) + geom_density(aes(x=z0$lambdas), color="black") + scale_x_continuous(limits=c(90, 110))
dev.off()

pdf('MNL_Scenario0.beta1.pdf', width = 8, height = 8)
ggplot(data=as.data.frame(z0$betas)) + geom_density(aes(x=V1), color="black") + scale_x_continuous(limits=c(0.1, 0.5))
dev.off()

pdf('MNL_Scenario0.beta2.pdf', width = 8, height = 8)
ggplot(data=as.data.frame(z0$betas)) + geom_density(aes(x=V2), color="black") + scale_x_continuous(limits=c(0, 0.3))
dev.off()






