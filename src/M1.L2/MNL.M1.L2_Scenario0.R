####
# Scenario 0. All choices are observed 
#
####

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs/src/M1.L2")

require("mvtnorm")
source(file="Metropolis-Hastings.R")



scenarioName <- "MNL.M1.L2_Scenario0"

## Load simulated choice data (NEED TO RUN MNL_InitData.binary.L2.R TO GENERATE THE DATA FIRST)
load(file="MNL.M1.L2_InitData.RData")

# final observation consists of Demand, and NoPurchase
observation0 <- list(demand=Demand, nopurchase=NoPurchase)



#log-posterior of betaT, to be called by M-H algorithm
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



## initialize input before sampling
# beta prior ~ N(beta.mu, beta.sg)
beta.mu <- rep(0, L+M)
beta.sg <- 100*diag(L+M)


nrun <- 5000
burnin <- 0.5
start <- burnin*nrun+1


### direct sampling beta
MH <- MH.mvnorm(logpost.betaT, start=c(-1, -1, -1), scale=c(0.1, 0.01, 0.02), nrun=nrun*10, thin=10, data=rbind(observation0$demand, observation0$nopurchase))
cat("MH acceptance rate: ", MH$accept, "\n")
betaTs <- MH$MC
betas <- cbind(betaTs[,1], betaTs[,2:(L+M)]*betaTs[,1])


### Estimate lambda by conjugate prior

# lambda ~ Gamma(lambda.alpha, lambda.beta)
lambda.alpha <- 0.005
lambda.beta <- 0.0001

#update posterior of lambda by conjugacy
lambda.alpha2 <- lambda.alpha + sum(observation0$demand, observation0$nopurchase)
lambda.beta2 <- lambda.beta + K



### save results
z0 <- list(lambdas=rgamma(nrun, shape=lambda.alpha2, rate=lambda.beta2), betas=betas)

save(z0, observation0, file=paste0(scenarioName, ".RData"))




### Visualize results

#plot lambda
samples.lambda <- z0$lambdas
plot(samples.lambda, type="l")

samples.lambda.truncated <- samples.lambda[start:nrun]
mean(samples.lambda.truncated)
quantile(samples.lambda.truncated, c(.025,.5,.975))
hist(samples.lambda.truncated)


#plot beta
samples.beta <- data.frame(z0$betas)
plot(samples.beta$X1, type="l")
plot(samples.beta$X2, type="l")
plot(samples.beta$X3, type="l")

samples.beta.truncated <- samples.beta[start:nrun,]
colMeans(samples.beta.truncated)
quantile(samples.beta.truncated$X1, c(.025,.5,.975))
quantile(samples.beta.truncated$X2, c(.025,.5,.975))
quantile(samples.beta.truncated$X3, c(.025,.5,.975))
hist(samples.beta.truncated$X1)
hist(samples.beta.truncated$X2)
hist(samples.beta.truncated$X3)


### save plots
require(ggplot2)

pdf(paste0(scenarioName, ".lambda.pdf"), width = 8, height = 8)
ggplot(data=data.frame(samples.lambda.truncated)) + geom_density(aes(x=samples.lambda.truncated), color="black")
dev.off()

pdf(paste0(scenarioName, ".beta.coef.pdf"), width = 8, height = 8)
ggplot(data=samples.beta.truncated) + geom_density(aes(x=X1), color="black")
ggplot(data=samples.beta.truncated) + geom_density(aes(x=X2), color="black")
dev.off()

pdf(paste0(scenarioName, ".beta.const.pdf"), width = 8, height = 8)
ggplot(data=samples.beta.truncated) + geom_density(aes(x=X3), color="black")
dev.off()
