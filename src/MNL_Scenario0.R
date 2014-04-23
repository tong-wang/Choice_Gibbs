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
        score <- apply(X_Mat, 1, function (x) exp(matrix(c(x,0,0), nrow=M, ncol=L, byrow=TRUE) %*% beta))
        choice.prob <- apply(score, 2, function(x) x/sum(x)) # M*N matrix
        
        logLikelihood <- data*log(choice.prob)
        
        logprior <- dmvnorm(log(beta), mean=beta.mu, sigma=beta.sg, log=TRUE)
        
        return(sum(logLikelihood) + logprior)
    }
}



# beta prior ~ logN(beta.mu, beta.sg)
beta.mu <- c(-2.5, -2.5)
beta.sg <- matrix(c(0.5, 0, 0, 0.5), 2, 2)



#direct sampling
z <- MCMCmetrop1R(logpost.beta, theta.init=rep(0.1, L),
                  data=observation0,
                  thin=1, mcmc=10000, burnin=2000, tune=0.01,
                  verbose=500,  V=matrix(c(1,0,0,1),2,2))


summary(z)
plot(z)
beta.estm <- colMeans(z)

### Estimate lambda by conjugate prior
## initialize input before sampling
# lambda ~ Gamma(lambda.alpha, lambda.beta)
lambda.alpha <- 0.5
lambda.beta <- 0.01

#update posterior of lambda by conjugacy
alpha2 <- lambda.alpha + sum(observation0)
beta2 <- lambda.beta + K

lambda.estm <- alpha2/beta2

score.estm <- exp(matrix(X_Mean, M-1,L) %*% beta.estm)
score.estm <- c(score.estm, 1)
choice.prob.estm <- score.estm / sum(score.estm)
demand.estm <- lambda.estm * choice.prob.estm




### save plots
require(ggplot2)

pdf('MNL_Scenario0.beta1.pdf', width = 8, height = 8)
ggplot(data=as.data.frame(z)) + geom_density(aes(x=V1), color="black")
dev.off()
pdf('MNL_Scenario0.beta2.pdf', width = 8, height = 8)
ggplot(data=as.data.frame(z)) + geom_density(aes(x=V2), color="black")
dev.off()






